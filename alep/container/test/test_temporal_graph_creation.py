# -*- python -*-
#
#       OpenAlea.image
#
#       Copyright 2012 INRIA - CIRAD - INRA
#
#       File author(s):  Jonathan Legrand <jonathan.legrand@ens-lyon.fr>
#                        Frederic Boudon <frederic.boudon@cirad.fr>
#
#       Distributed under the Cecill-C License.
#       See accompanying file LICENSE.txt or copy at
#           http://www.cecill.info/licences/Licence_CeCILL-C_V1-en.html
#
#       OpenAlea WebSite: http://openalea.gforge.inria.fr
#
################################################################################
"""Test creation of TemporalPropertyGraph"""

from temporal_property_graph_input import create_random_TPG
from openalea.container import TemporalPropertyGraph
import networkx as nx


def test_export_TPG_2_networkx():
    """ Create a random graph and export it to networkx """
    g = create_random_TPG()
    nxg = g.to_networkx()

def test_import_TPG_from_networkx():
    """ Create a random graph and import it to networkx """
    g = create_random_TPG()
    nxg = g.to_networkx()
    gg = TemporalPropertyGraph().from_networkx(nxg)
    #~ assert g==gg

def test_TGP_display_by_networkx(display = False):
    """ Test the display of a TPG  by networkx """
    g = create_random_TPG()
    nxg = g.to_networkx()
    import matplotlib.pyplot as plt
    nx.draw(nxg)
    if display:
        plt.show()


