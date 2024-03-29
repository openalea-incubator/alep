# -*- coding: utf-8 -*-
__revision__ = "$Id: $"

import sys
import os

# CPL 
#from setuptools import setup, find_packages
from setuptools import setup, find_namespace_packages

# Remove metainfo and set it directly in the setup
#from openalea.deploy.metainfo import read_metainfo

# Reads the metainfo file
#metadata = read_metainfo('metainfo.ini', verbose=True)
#for key,value in metadata.items():
#    exec("%s = '%s'" % (key, value))
name='openalea.alep'

_version = {}
with open("src/alinea/alep/version.py") as fp:
    exec(fp.read(), _version)
version = _version["__version__"]

short_descr = "OpenAlea.Alep is a generic pathosystem model simulating 3D plant / microclimate and foliar pathogen interactions."
version=version
description=short_descr
long_description=description
authors="Guillaume Garin, Marc Labadie, Christian Fournier, Corinne Robert, Christophe Pradal",
authors_email="Christian.Fournier@inrae.fr, christophe dot pradal _at_ cirad fr ",
url="https://github.com/openalea-incubator/alep"
license="CeCILL-C"
keywords = 'openalea, FSPM, plant pathosystem',	



# Packages list, namespace and root directory of packages

packages = find_namespace_packages(where='src', include=['alinea.*'])
package_dir={'': 'src'}


# dependencies to other eggs
setup_requires = ['openalea.deploy']

setup(
    name=name,
    version=version,
    description=description,
    long_description=long_description,
    author=authors,
    author_email=authors_email,
    url=url,
    license=license,
    keywords = keywords,	

    # package installation
    packages= packages,	
    package_dir= package_dir,

    # Namespace packages creation by deploy
    #namespace_packages = ['alinea'],
    #create_namespaces = False,
    zip_safe= False,

    # Dependencies
    setup_requires = setup_requires,


    # Eventually include data in your package
    # (flowing is to include all versioned files other than .py)
    include_package_data = True,
    # (you can provide an exclusion dictionary named exclude_package_data to remove parasites).
    # alternatively to global inclusion, list the file to include   
    #package_data = {'' : ['*.pyd', '*.so'],},

    # postinstall_scripts = ['',],

    # Declare scripts and wralea as entry_points (extensions) of your package 
    entry_points = { 
        'wralea' : ['alep = alinea.alep_wralea'],
        'alep.disease' :  ['template = alinea.alep.fungal_objects:Fungus',
                           'septoria= alinea.alep.septoria:SeptoriaFungus',
                           'powdery_mildew = alinea.alep.powdery_mildew:PowderyMildewFungus',
                           'brown_rust = alinea.alep.brown_rust:BrownRustFungus',],
        },
    )


