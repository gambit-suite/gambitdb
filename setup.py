import glob
from setuptools import setup

# Minimal setup.py to handle scripts only
# All other configuration is in pyproject.toml
setup(
    scripts=[s for s in glob.glob('scripts/*') if not s.endswith('.txt')]
)
