Matlab code and data for reproducing the results in the paper:
["Efficient Regularized Field Map Estimation in 3D MRI"](http://doi.org/10.1109/TCI.2020.3031082),
IEEE Trans. Computational Imaging</em>, 6:1451-8, 2020,
by Claire Lin and Jeffrey A. Fessler, University of Michigan.

This code requires the Matlab version of the Michigan Image Reconstruction Toolbox (MIRT)
from [http://web.eecs.umich.edu/~fessler/code/index.html]
or [https://github.com/JeffFessler/mirt].
Please set up MIRT before running the examples.

The following scripts are in the example folder:
1. example_simu_fieldmap.m: Figs. 3,4 
(Same procedure is used to reproduce Figs. 5,6, except using ESPIRiT for sensemap.)
2. example_simu_waterfat.m: Figs. 7,8,9
(Simulation data generation requires 
[ISMRM fat-water toolbox](https://www.ismrm.org/workshops/FatWater12/data.htm).)
3. example_ankle_waterfat.m: Figs. 10,11
