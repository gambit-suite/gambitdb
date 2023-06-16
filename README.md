# GambitDB

This is a tool to generate a Gambit database. Its primary input is a spreadsheet from GTDB such as this [release](https://data.gtdb.ecogenomic.org/releases/release214/214.0/bac120_metadata_r214.tar.gz).

## Installation
This software can be installed using pip:
```
pip install .
```
Additional dependancies are:
* Gambit
* ncbi-genome-download

## Usage

To run the analysis, execute the `run_gambitdb.sh` script:

```bash
./run_gambitdb.sh /path/to/working_directory /path/to/spreadsheet.tsv num_cores
```

This script will create a Gambit database from a GTDB spreadsheet. It parses the spreadsheet, downloads data with ncbi-genome-download

## Results

The results of the analysis will be output to '/path/to/working_directory/final'

## Contributing

Contributions to this project are welcome. To contribute, please fork the repository and submit a pull request.

## License

This project is licensed under the GNU GPL 3 License - see the [LICENSE](LICENSE) file for details.