import gcslab
import numpy as np
from math import pi

grating = gcslab.GCSlabModel()

# Define useful parameters
wavelength = 1.55
angle = 8.0/180*pi

n_Si = 3.45
n_SiO2 = 1.44

h_Si = 0.22
h_etch = 0.11
h_BOX = 2.0
df = 0.8

# Build the grating layer-by-layer. For each layer, we specify the thickness
# and refractive index. At least one layer must also have a non-zero source
# amplitude assigned to it
grating.add_slab(np.inf, n_SiO2)
grating.add_slab(h_etch, n_Si*df + n_SiO2*(1-df), source_amplitude=1.0)
grating.add_slab(h_Si-h_etch, n_Si)
grating.add_slab(h_BOX, n_SiO2)
grating.add_slab(np.inf, n_Si)

# Visualize the layer stack
grating.visualize_stack() 

# Calculate the directionality
directionality = grating.directionality(wavelength, angle)
print(f'The directionality = {directionality}.')

