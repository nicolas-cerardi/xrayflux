# xrayflux

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

A package that wraps `pyatomdb`, performs spectrum redshifting, optionally applies instrumental response, and outputs the X-ray photon flux.
The main functions are:
 - `compute_flux`: returns the photon flux for a given temperature $T$, abundancy $Z$ and redshift $z$.
 - `make_fluxtable`: returns tabulated photon fluxes in the $T - z$ dimensions, for a fixed $Z$.

### Installation

```
git clone https://github.com/nicolas-cerardi/xrayflux.git
cd xrayflux
pip install .
```

### Minimal example

_using AtomDB 3.0.9 tables._

```python
from xrayflux.xraytables import compute_flux

T = 5 #in keV
z = 0.
flux = compute_flux(T=T, z=z)
print(flux, "ph cm^3 s-1")
>>> 2.672476983729378e-15 ph cm^3 s-1
```

For more see the [demo notebook](https://github.com/nicolas-cerardi/xrayflux/blob/main/notebooks/demo_xrayflux.ipynb).