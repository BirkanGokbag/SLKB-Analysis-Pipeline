## Installation

### Installing from PyPI (Recommended - Upon Release)

```
pip install SLKB
```

### Creating from source

After creating your desired environment, run the following under package folder

```
python setup.py sdist
pip install dist/SLKB-1.0.0.tar.gz --user
```

### Requirements

The project requires several packages to work appropriately. These packages and their versions can be found within requirements.txt or alternatively down below.

```
pandas==1.5.3
numpy==1.21.0
SQLAlchemy==2.0.12
scipy==1.7.3
ipykernel==6.9.1
ipython==8.2.0
ipython-genutils==0.2.0
ipywidgets==7.6.5
jupyterlab==3.3.2
mysql-connector-python==8.0.29
```