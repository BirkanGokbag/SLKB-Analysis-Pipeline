## Installation

## Creating from source

After creating your desired environment, run the following under package folder

```
python setup.py sdist
pip install dist/DGCyTOF-1.0.0.tar.gz --user
```

## Requirements

The project requires several packages to work appropriately. These packages and their versions can be found within requirements.txt or alternatively down below.

```
scipy==1.5.2
numpy==1.19.1
seaborn==0.10.1
matplotlib==3.3.0
pandas==1.1.0
umap_learn==0.4.6
hdbscan==0.8.26
scikit_learn==0.23.2
scikit_plot==0.3.7
scikit_bio==0.5.6
umap-learn==0.4.6
torch==1.6.0
torchvision==0.7.0
tensorboard==2.3.0
```