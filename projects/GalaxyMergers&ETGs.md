---
layout: project
type: project
image: img/position_angle_plot_crop.png
title: "GalaxyMergers&ETGs"
date: 2024-08-15
published: false
labels:
  - Early-Type Galaxies
  - Python
  - Firefly
summary: "This study investigates the impact of galaxy mergers on the stellar initial mass function (IMF) of early-type galaxies (ETGs), using stacked spectra from the MaNGA-DR17 survey and SDSS Stripe 82 region"
---
We examine the effects of galaxy mergers on the stellar initial mass function (IMF) of early-type galaxies (ETGs) using stacked spectra of galaxies in the MaNGA-DR17 survey and the Stripe 82 region of the Sloan Digital Sky Survey. Galaxy mergers leave remnants of their gravitational interactions behind as tidal features, which are structures such as streams, shells and tails. A previous study revealed that ETGs with constant IMFs that have undergone recent mergers have lower metallicities, shallower metallicity gradients (at $M_* > (M_{\odot})^{10.6}$), more positive age gradients (at $M_* > (M_{\odot})^{10.6}$) and younger stellar populations (at $M_* < M_{\odot}^{11.1}) than ETGs that have no tidal features [1]. Other recent work has confirmed that allowing for gradients in the stellar IMF returns considerably larger gradients in the stellar mass-to-light ratio than if the IMF is held fixed [2]. IMF-driven stellar mass-to-light gradients can have extreme effects for how one estimates the stellar masses of ETGs. Using a variable IMF, this study aims to investigate the impact of galaxy mergers on the shape of the IMF of ETGs.

The spectral analysis code used in this project was originally developed in IDL by Professor Mariangela Bernardi. This code enables spectral stacking, smoothing, Lick index calculation, and stellar population model fitting to determine best-fitting parameters. I adapted and applied this code to investigate the impact of galaxy mergers on the IMF of ETGs.

[1] Yoon, Yongmin, Jongwan Ko, and Jae-Woo Kim. "Impact of Galaxy Mergers on Stellar Population Profiles of Early-type Galaxies." The Astrophysical Journal 946, no. 1 (2023): 41. doi: 10.3847/1538-4357/acbcc5
[2] Bernardi, M., R. K. Sheth, H. Dominguez-Sanchez, J.-L. Fischer, K.-H. Chae, M. Huertas-Company, and F. Shankar. "M*/L gradients driven by IMF variation: large impact on dynamical stellar mass estimates." Monthly Notices of the Royal Astronomical Society 477, no. 2 (2018): 2560-2571. doi: 10.1093/mnras/sty781.
