#!/usr/bin/env python
# -*- coding: utf-8 -*-

import io
import os
from setuptools import setup

DEPENDENCIES = ['pandas==1.5.3',
'numpy==1.21.0',
'SQLAlchemy==2.0.12',
'scipy==1.7.3',
'ipykernel==6.9.1',
'ipython==8.2.0',
'ipython-genutils==0.2.0',
'ipywidgets==7.6.5',
'jupyterlab==3.3.2',
'mysql-connector-python==8.0.29']

CURDIR = os.path.abspath(os.path.dirname(__file__))

with io.open(os.path.join(CURDIR, "README.md"), "r", encoding="utf-8") as f:
    README = f.read()

setup(
    name="SLKB",
    version="1.0.11",
    author="Birkan Gökbağ",
    author_email="birkan.gokbag@gmail.com",
    description="SLKB: Synthetic lethality knowledge base for gene combination double knockout experiments",
    long_description=README,
    url="https://github.com/BirkanGokbag/SLKB-Analysis-Pipeline",
    package_dir={'SLKB': 'SLKB'},
    packages=['SLKB'],
    include_package_data=True,
    keywords=[],
    scripts=[],
    zip_safe=False,
    install_requires=DEPENDENCIES,
    license="License :: OSI Approved :: GPL 3.0",
    classifiers=[
        "Programming Language :: Python",
        "Operating System :: OS Independent",
    ],
)