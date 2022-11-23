https://github.com/ClaireYLin/regularized-field-map-estimation

Matlab code and data for reproducing the results in the paper:
["Efficient Regularized Field Map Estimation in 3D MRI"](http://doi.org/10.1109/TCI.2020.3031082),
IEEE Trans. Computational Imaging, 6:1451-8, 2020,
by Claire Lin and Jeffrey A. Fessler, University of Michigan.

This code requires the Matlab version of the Michigan Image Reconstruction Toolbox (MIRT)
from [http://web.eecs.umich.edu/~fessler/code/index.html]
or [https://github.com/JeffFessler/mirt].
Please set up MIRT before running the examples.

The following scripts are in the example folder:
1. `example_simu_fieldmap.m`: Figs. 3,4
(Same procedure is used to reproduce Figs. 5,6,
except using ESPIRiT for sensemap.)
2. `example_simu_waterfat.m`: Figs. 7,8,9
(Simulation data generation requires
[ISMRM fat-water toolbox](https://www.ismrm.org/workshops/FatWater12/data.htm).)
3. `example_ankle_waterfat.m`: Figs. 10,11


### Errata

The original version of the code was missing a factor of `π`
in the calculation of RMSD in Hz.
See
[Issue #5](https://github.com/ClaireYLin/regularized-field-map-estimation/issues/5).
We have
[updated the code](https://github.com/ClaireYLin/regularized-field-map-estimation/pull/9)
to correct that calculation,
which leads to even lower errors
than were shown in the paper.
We also modified the code
to compute the RMSE
over the mask,
instead over the entire image.
Here is the updated Figure 4.

![fig4new](https://github.com/ClaireYLin/regularized-field-map-estimation/blob/main/fig/fig4new.png)


### Related packages

A Julia version of this code is available at
https://github.com/MagneticResonanceImaging/MRIfieldmaps.jl.
Version `v0.0.1` of that repo
is Julia code that basically reproduces Fig. 4.
Later versions offer improvements in speed and convenience.
