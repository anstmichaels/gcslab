# GCSlab
Quickly estimate the directionality of a grating coupler using a homogenized-slab model. 

## Prerequisites
Before using GCSlab, the following python libraries must be installed:
* numpy
* Matplotlib

Furthermore, currently gcslab only supports python3.

## Example
The following code snippet models a simple partially-etch silicon grating coupler. The
grating coupler is defined by stacking up slabs of different refractive index. The
etched region of the grating is assigned an average refractive index and a non-zero
source amplitude. In this region, waves are injected and the total fraction of power
leaving the top of the grating is estimated.

```
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

```

This example can be found in the `example_gc.py`.

# Citing this Code
This code is based on the doctoral thesis titled "A Hierarchical Approach to the
Design and Optimization of Photonics." If using this code in a publication,
please reference the afformentioned work.
