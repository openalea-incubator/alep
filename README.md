# Alinea.Alep

**Authors** : XXX

**Institutes** : INRIA / CIRAD

**Status** : Python package

**License** : Cecill-C

**URL** : https://github.com/openalea-incubator/alep.git

## Description

Alinea.Alep is a XXX

## Content

The Alinea.Alep package contains :

## Installation

### Requirements

* python = 3.8
* setuptools = 49.6
* numpy
* pandas
* scipy
* matplotlib
* openalea.deploy
* openalea.core
* openalea.mtg
* openalea.plantgl
* alinea.adel

### Users installation

```cmd
conda create -n alep -c conda-forge python=3.8 setuptools=49.6 numpy pandas scipy matplotlib

```

Activate conda environement

```cmd

conda activate alep
```

Install packages from conda

```cmd
conda install -c fredboudon -c conda-forge openalea.plantgl
```

Install packages from github repositories

* deploy

```cmd
git clone https://github.com/openalea/deploy.git
cd deploy
python setup.py install
cd..
```

* core

```cmd
git clone https://github.com/openalea/core.git
cd core
python setup.py install
cd..
```

* alinea.mtg

```cmd
git clone https://github.com/openalea/mtg.git
cd mtg
python setup.py install
cd..
```

* alinea.adel

```cmd
git clone https://github.com/openalea-incubator/adel.git
cd adel
python setup.py install
cd..
```

* alinea.astk

```cmd
git clone https://github.com/openalea-incubator/astk.git
cd astk
python setup.py install
cd..
```

* alinea.alep

```cmd
git clone https://github.com/openalea-incubator/alep.git
cd alep
python setup.py install
```

## Quick start


## Documentation

http://openalea.gforge.inria.fr/doc/alinea/alep/doc/_build/html/contents.html

