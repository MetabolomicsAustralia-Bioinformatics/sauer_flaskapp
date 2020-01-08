# Web app for Sauer Enrichment

Built in `flask`, which admittedly is a slightly overblown framework (Python `dash` would have been fine). 

## Local Setup

Sample input data files and reference output files available in `/sample_files`. 

```
# install virtualenv package
pip3 install virtualenv

# Create virtual env
python3 -m venv sauer-env

# Activate virtual env
source sauer-env/bin/activate

# Deactivate virtual env when done
deactivate
```

### Requirements

Requirements in `reauirements.txt`. 
Dev note: automatically collect packages installed in the virtual environment to create `requirements.txt` with `gunicorn`:

```
pip3 install gunicorn
pip3 freeze > requirements.txt
```