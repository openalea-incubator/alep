
# This file has been generated at Wed Aug 01 09:27:07 2012

from openalea.core import *


__name__ = 'Alinea.SimCycle.pesticide'

__editable__ = True
__description__ = 'This module implements the model of Milne et al (2007) of the efficay of a (mix of) pesticide compounds'
__license__ = 'CeCILL-C'
__url__ = 'http://openalea.gforge.inria.fr'
__alias__ = []
__version__ = '0.1.0'
__authors__ = ''
__institutes__ = None
__icon__ = ''


__all__ = ['pesticide_nodes_decay', 'pesticide_nodes_erode_products', '_113433296', 'pesticide_nodes_global_efficacy', 'pesticide_nodes_compounds_from_csv']

default_compounds = [{'compound': 'Epoxiconazole', 'Ap': 0.71, 'Kp': 6.0, 'Ae': 0.5, 'type_code': 0, 'Ke': 7.01, 'decay_rate': 0.069, 'dose_max_ha': 125}]
default_doses = {'Epoxiconazole': 0.0125}
default_erosion_rates = {'Epoxiconazole' : 0.5}

pesticide_nodes_decay = Factory(name='decay',
                authors=' (wralea authors)',
                description='Update doses after decay over time',
                category='model',
                nodemodule='pesticide_nodes',
                nodeclass='decay',
                inputs=[{'interface': IDict, 'name': 'doses', 'value': default_doses, 'desc': ''}, 
                        {'interface': IFloat, 'name': 'days', 'value': 1.0, 'desc': ''}, 
                        {'interface': ISequence, 'name': 'compound parameters', 'value': default_compounds, 'desc': ''}],
                outputs=[{'interface': IDict, 'name': 'updated doses', 'desc': ''}],
                widgetmodule=None,
                widgetclass=None,
               )

pesticide_nodes_erode_products = Factory(name='erode products',
                authors=' (wralea authors)',
                description='Compute loss of activity of compounds due to evolution of resistance of strains',
                category='model',
                nodemodule='pesticide_nodes',
                nodeclass='erode_products',
                inputs=[{'interface': IDict, 'name': 'erosion_rates', 'value': default_erosion_rates, 'desc': ''}, 
                        {'interface': ISequence, 'name': 'compound parameters', 'value': default_compounds, 'desc': ''}],
                outputs=[{'interface': ISequence, 'name': 'updated parameters', 'desc': ''}],
                widgetmodule=None,
                widgetclass=None,
               )


_113433296 = DataFactory(name='products_parameters.csv',
                    description='Sample parameter file for tutorial of pesticide module',
                    editors=None,
                    includes=None,
                    )



pesticide_nodes_global_efficacy = Factory(name='global efficacy',
                authors=' (wralea authors)',
                description='compute global efficacy of a mixture of coumpounds',
                category='model',
                nodemodule='pesticide_nodes',
                nodeclass='global_efficacy',
                inputs=[{'interface': IDict, 'name': 'doses', 'value': default_doses, 'desc': ''}, 
                        {'interface': ISequence, 'name': 'compounds parameters', 'value': default_compounds, 'desc': ''}],
                outputs=[{'interface': IDict, 'name': 'efficacy of products', 'desc': ''}],
                widgetmodule=None,
                widgetclass=None,
               )




pesticide_nodes_compounds_from_csv = Factory(name='compounds from csv',
                authors=' (wralea authors)',
                description='',
                category='data i/o',
                nodemodule='pesticide_nodes',
                nodeclass='compounds_from_csv',
                inputs=[{'interface': IFileStr(filter="*.csv", save=False), 'name': 'filename', 'value': None, 'desc': ''}, {'interface': IStr, 'name': 'separator', 'value': ';', 'desc': ''}],
                outputs=[{'interface': ISequence, 'name': 'compound list', 'desc': ''}],
                widgetmodule=None,
                widgetclass=None,
               )




