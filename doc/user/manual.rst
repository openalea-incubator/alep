Package Layout
###############

The Alep model is implemented in the directory src/alinea/alep.
The OpenAlea components are in src/alinea/alep_wralea

Main files are:
    * cycle2.py:
        - define the functioning of a dispersal_unit and a lesion
        - examples are provided for wheat septoria and grapevine powdery mildew
    * protocol.py: 
        - Generic algorithms used to simulate an epidemic on the canopy
    * Septo3DDispersion.py
        - Example of model used to disperse spores of septoria over a wheat canopy

Deprecated files are:
    * cycle.py
    * DispersalUnit.py
