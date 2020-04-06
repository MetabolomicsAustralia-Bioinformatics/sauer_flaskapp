# Web app for Sauer Enrichment

Built in `flask`, which admittedly is a slightly overblown framework (Python `dash` would have been fine). Sauer enrichment logic copy-pasted from the [SauerEnrichmentAnalysis repo](https://github.com/MetabolomicsAustralia-Bioinformatics/SauerEnrichmentAnalysis), with some minor (mostly cosmetic) modifications. 

## Local Setup

Sample input data files and reference output files available in `/sample_files`. 

```
# install virtualenv package
pip3 install virtualenv

# Create virtual env
python3 -m venv venv

# Activate virtual env
source venv/bin/activate

# Install all required packages
pip install -r requirements.txt

# Run the app on localhost (port 5000 by default)
python3 app.py

# Deactivate virtual env when done
deactivate
```

### Making Requirements

Automatically collect packages installed in the virtual environment to create `requirements.txt` with `gunicorn`:

```
pip3 install gunicorn
pip3 freeze > requirements.txt
```
