"""Manual definition of plugins (can replace plugin discovery when unavailable"""

from alinea.alep.septoria import SeptoriaFungus
from alinea.alep.powdery_mildew import PowderyMildewFungus
from alinea.alep.brown_rust import BrownRustFungus


def get_disease(model='septoria'):
    if model == 'septoria':
        return SeptoriaFungus()
    elif model == 'powdery_mildew':
        return PowderyMildewFungus()
    elif model == 'brown_rust':
        return BrownRustFungus()
    else:
        raise Exception('Unknown disease: ' + model + ', should be one of : septoria, brown_rust, powdery_mildew')
