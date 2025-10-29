import random
import requests
from requests.adapters import HTTPAdapter
from urllib3.util.retry import Retry
from tqdm import tqdm
import os
import time
import logging
from typing import Optional, Any, Dict
import socket
import urllib3

#Hepler client class to interact with NCBI Datasets API, biggest concern was around retries and rate limiting
#Continue to add to this class as needed

class NCBIDatasetClient:
    """
    Handles all NCBI API interactions with advanced retry logic and error handling.
    """
    
    def __init__(self, 
                 rate_limit: float = 0.1,
                 max_retries: int = 2,
                 backoff_factor: float = 0.5,
                 retry_statuses: tuple = (408, 429, 500, 502, 503, 504),
                 timeout: tuple = (10, 30)):
        """
        Initialize the NCBI API client with retry configuration.
        
        Args:
            rate_limit: Minimum time between requests in seconds
            max_retries: Maximum number of retry attempts
            backoff_factor: Exponential backoff factor between retries
            retry_statuses: HTTP status codes to retry on
            timeout: (connect timeout, read timeout) in seconds
        """
        self.logger = logging.getLogger(__name__)
        self.rate_limit = rate_limit
        self.last_request_time = 0
        # Timeout is a tuple of (connect, read) timeouts
        self.timeout = timeout
        
        # Configure retry strategy
        retry_strategy = Retry(
            total=max_retries,
            backoff_factor=backoff_factor,
            status_forcelist=retry_statuses,
            allowed_methods=["HEAD", "GET", "POST", "PUT", "DELETE", "OPTIONS", "TRACE"],
            raise_on_status=False,
            respect_retry_after_header=True
        )
        
        self.session = requests.Session()
        adapter = HTTPAdapter(max_retries=retry_strategy)
        self.session.mount("https://", adapter)
        self.session.mount("http://", adapter)
    
    def _wait_for_rate_limit(self):
        """Enforce rate limiting between API calls"""
        elapsed = time.time() - self.last_request_time
        if elapsed < self.rate_limit:
            time.sleep(self.rate_limit - elapsed)
        self.last_request_time = time.time()
    
    def _handle_network_errors(self, url: str, error: Exception, attempt: int, max_attempts: int) -> None:
        """Log network-related errors with appropriate context"""
        if isinstance(error, (requests.exceptions.ConnectionError, urllib3.exceptions.NewConnectionError)):
            self.logger.warning(
                f"Network connection error on attempt {attempt}/{max_attempts} "
                f"for URL {url}. Will retry after backoff. Error: {str(error)}"
            )
        elif isinstance(error, requests.exceptions.Timeout):
            self.logger.warning(
                f"Request timed out on attempt {attempt}/{max_attempts} "
                f"for URL {url}. Will retry after backoff. Error: {str(error)}"
            )
        elif isinstance(error, socket.gaierror):
            self.logger.warning(
                f"DNS resolution failed on attempt {attempt}/{max_attempts} "
                f"for URL {url}. Will retry after backoff. Error: {str(error)}"
            )
        else:
            self.logger.warning(
                f"Unexpected error on attempt {attempt}/{max_attempts} "
                f"for URL {url}. Will retry after backoff. Error: {str(error)}"
            )
    
    def make_request(self, 
                    url: str, 
                    method: str = 'GET', 
                    max_attempts: int = 5,
                    headers: Optional[Dict[str, str]] = None,
                    **kwargs) -> requests.Response:
        """
        Make an API request with advanced error handling and retries.
        
        Args:
            url: The URL to request
            method: HTTP method to use
            max_attempts: Maximum number of attempts for this specific request
            **kwargs: Additional arguments to pass to requests
        
        Returns:
            requests.Response object
        
        Raises:
            requests.exceptions.RequestException: If all retry attempts fail
        """
        self._wait_for_rate_limit()
        
        # Ensure timeout is set
        kwargs['timeout'] = kwargs.get('timeout', self.timeout)
        
        request_headers = {}
        if headers:
            request_headers.update(headers)
        
        last_exception = None
        for attempt in range(max_attempts):
            try:
                response = self.session.request(method, url, headers=request_headers, **kwargs)
                
                # If the request was successful, return the response
                if response.ok:
                    return response
                
                # Handle specific status codes commonly observed in datasets api
                if response.status_code == 429:  # Too Many Requests
                    retry_after = int(response.headers.get('Retry-After', 60))
                    self.logger.warning(f"Rate limit hit. Waiting {retry_after} seconds...")
                    time.sleep(retry_after)
                    continue
                    
                if response.status_code >= 500:
                    self.logger.warning(
                        f"Server error {response.status_code} on attempt {attempt + 1}/{max_attempts}. "
                        f"URL: {url}"
                    )
                    time.sleep((2 ** attempt) + random.uniform(0, 1))
                    continue
                
                # If we get here, it's an unexpected error
                response.raise_for_status()
                
            except (requests.exceptions.RequestException,
                   urllib3.exceptions.NewConnectionError,
                   socket.gaierror) as e:
                last_exception = e
                self._handle_network_errors(url, e, attempt + 1, max_attempts)
                
                if attempt < max_attempts - 1:
                    # Don't want thundering herd problem
                    backoff = (2 ** attempt) + random.uniform(0, 1)
                    self.logger.info(f"Backing off for {backoff:.2f} seconds before retry...")
                    time.sleep(backoff)
                    continue
                break
        
        # This means all the retries have failed
        error_msg = f"Failed after {max_attempts} attempts for URL: {url}"
        if last_exception:
            error_msg += f". Last error: {str(last_exception)}"
        raise requests.exceptions.RequestException(error_msg)
    
    def get_json(self, url: str, **kwargs) -> Dict[str, Any]:
        """
        Convenience method to make a request and return JSON data.
        
        Args:
            url: The URL to request
            **kwargs: Additional arguments to pass to make_request
            
        Returns:
            Parsed JSON response
            
        Raises:
            requests.exceptions.RequestException: If request fails
            json.JSONDecodeError: If response is not valid JSON
        """
        response = self.make_request(url, **kwargs)
        return response.json()
    
    def download_file(self, 
                     url: str, 
                     output_path: str, 
                     chunk_size: int = 9000,
                     **kwargs) -> None:
        """
        Download a file with progress tracking and verification.
        
        Args:
            url: URL to download from
            output_path: Where to save the file
            chunk_size: Size of chunks to download
            **kwargs: Additional arguments to pass to make_request
        """
        response = self.make_request(url, stream=True, **kwargs)
        total_size = int(response.headers.get('content-length', 0))
        
        with open(output_path, 'wb') as f, tqdm(
            desc=f"Downloading {output_path}",
            total=total_size,
            unit='B',
            unit_scale=True,
            unit_divisor=1024,
        ) as pbar:
            for chunk in response.iter_content(chunk_size=chunk_size):
                size = f.write(chunk)
                pbar.update(size)
        
        # Want to make sure the file is not corrupted
        actual_size = os.path.getsize(output_path)
        if total_size > 0 and actual_size < total_size * 0.8:  # Allow 20% margin, this is observational 
            raise ValueError(
                f"Downloaded file size ({actual_size} bytes) is significantly "
                f"smaller than expected ({total_size} bytes)"
            )