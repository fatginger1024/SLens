# SLens
## Introduction

Galaxy-galaxy strong lensing simualtion package. Assuming all lenses are axis-symmetric, whose mass density profiles can be modelled by a gnfw+Sersic model.
## Data
Lens galaxy catalogue:

[Crocce, M., Castander, F. J., Gaztañaga, E., Fosalba, P., & Carretero, J. 2015,
MNRAS, 453, 1513](https://arxiv.org/abs/1312.2013)

Source galaxy cataglogue:

[Laigle, C., McCracken, H. J., Ilbert, O., et al. 2016, ApJS, 224, 24](https://arxiv.org/abs/1604.02350)

## Basic models

- gNFW:

refer: Resolving the Central Density Profile of Dark Matter Halos Gravitational Lensing Statistics by Masamune Oguri, section E.1.5
- Sersic:

refer: [Sonnenfeld, A. & Cautun, M. 2021, A&A, 651, A18](https://arxiv.org/abs/2102.08973)
section 3.1
- concentration:

[Ludlow, A. D., Bose, S., Angulo, R. E., et al. 2016, MNRAS, 460, 1214](https://arxiv.org/abs/1601.02624) appendix c
- mass-size:

[Sonnenfeld, A., Wang, W., & Bahcall, N. 2019, A&A, 622, A30](https://arxiv.org/abs/1811.04934) section 3.1
- angular diamter distance:

[Hogg, David W. "Distance measures in cosmology." arXiv preprint astro-ph/9905116 (1999).](https://arxiv.org/abs/astro-ph/9905116) section 6


## Modules

![plot](./docs/mod1.jpg)
![plot](./docs/mod2.jpg)
![plot](./docs/mod3.jpg)
![plot](./docs/mod4.jpg)

## Workflow

![plot](./docs/workflow.png)

## Analysing composite (gnfw+Sersic) lens statistics

One can analyse the lens statistics using the <code>lens_statistics</code> module as provided in the tutorials.
The lens statistics (as functions of the dimensionless radial coordiate <img src="https://render.githubusercontent.com/render/math?math=x"> ) are:
- <img src="https://render.githubusercontent.com/render/math?math=\Sigma(x)">: the projected mass density;
- <img src="https://render.githubusercontent.com/render/math?math=\alpha(x)">: the deflection angle;
- <img src="https://render.githubusercontent.com/render/math?math=\beta(x)">: the lens equation, which describes the mapping from the lens plane <img src="https://render.githubusercontent.com/render/math?math=x"> to the source plane <img src="https://render.githubusercontent.com/render/math?math=\beta">;
- <img src="https://render.githubusercontent.com/render/math?math=\kappa(x)">: the convergence, which describes the magnification of the image by increasing its size;
- <img src="https://render.githubusercontent.com/render/math?math=\gamma(x)">: the shear, which describes the how much the shape of the image is changed tangentially.

<img src="./plots/Sigma.png" width="512"/>
<img src="./plots/alpha.png" width="512"/>
<img src="./plots/beta.png" width="512"/>
<img src="./plots/kappa.png" width="512"/>
<img src="./plots/gamma.png" width="512"/>
