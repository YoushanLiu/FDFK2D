![Lanuguage](https://img.shields.io/badge/-Fortran-734f96?logo=fortran)
![Version](https://img.shields.io/static/v1?label=version&message=0.8&color=blue)
![Downloads](https://img.shields.io/github/downloads/YoushanLiu/FDFK2D/total)
![Commits](https://img.shields.io/github/commit-activity/m/YoushanLiu/FDFK2D)
![Size](https://img.shields.io/github/languages/code-size/YoushanLiu/FDFK2D)
![License](https://img.shields.io/github/license/YoushanLiu/FDFK2D)
# FDFK2D
Efficient two-dimensional teleseismic wavefiled hybrid simulation method for receiver function analysis. This software is a highly optimized Fortran package supporting parallel computing using the OpenMP, which can be used to synthetic P- and S-waves RFs.

# Background
Until now, simulation of high-frequency (< 8 s) seismic wavefield in global-scale is still challenging even using large-scale computer clusters. One feasible way to mitigate this issue is to adopt a so-called hybrid method. This method usually computes the wavefield from source to the boundaries of the target heterogeneous media using a fast and efficient 1D wavefield solver (such as AxiSEM, DSM, general Ray method, Normal mode, FK, etc.), which is referred to background wavefield. Then, the complex wavefield in the localized target heterogeneous media is computed using 2D/3D wavefield solver (such as finite difference method, pseduospectral method, finite element method, spectral element method, etc.) with the background wavefield as incident wavefield, which is referred to complete wavefield. In the boundaries of the target heterogenous model, a interface between background wavefield and complete wavefiled is applied. This method can simulate the interactions of incoming teleseismic wavefield with local heterogeneities but reducing the computational region to much smaller local domain, which can significantly reduce the computing costs of high-frequency teleseismic wavefield with less computational resources. 

![Schematic diagram of hybrid method](https://github.com/YoushanLiu/FDFK2D/blob/master/images/Hybrid%20method.png)

# About
Purpose: This software is a hybrid method coupling the FD and FK methods to simulate teleseismic wavefield.

FD: Seismic wavefield modeling using collocated grid finite difference method in two-dimensional elastic isotropic meida.

FK: Seismic wavefield modeling using the frequency-wavenumber or propagator matrix method in three-dimensional elastic isotropic media.

# Citation
If you use this software in your researches and prepare publications, please cite at least one of the following papers:

- Youshan Liu, Chenglong Wu, Tao Xu, Liang Zhao, Jiwen Teng. 2024. FDFK2D: Efficient two-dimensional teleseismic wavefield modeling
   for receiver function analysis using a hybrid method. Seismological Research Letters, accepted.
- Chenglong Wu, Tao Xu, Xiaobo Tian, Ross N. Mitchell, Jiyan Lin, Jianfeng Yang, Xin Wang, Zhanwu Lu. 2024. Underthrusting of Tarim
   lower crust beneath the Tibetan Plateau revealed by receiver function imaging. Geophysical Research Letters, 51, e2024GL108220.

# Prerequisites
Before you compile this software, you should first install the following packages.
- blas
- lapack
- fftw3

Usually, these packages can be installed by the Package Manager of your system.

# Installation and Usage
Before you installing this software, you can type `make lib` to copy the library files in the `lib` folder to `/usr/lib64` with root permission, then type `exit`.
Or, you can add the directory to the environment variable by setting `export LD_LIBRARY_PATH=/your_path:$LD_LIBRARY_PATH`.
Now, you can install this software as follows.
- `cd ./src`
- `make install`
   - It is better to use the user permission. It will install the executable file in ~/bin, thus you can run FDFK2D anywhere.
- `cd ../examples/test1`

Then, you can run the tests.
- `./run_this_example`

# Documentation
All input files are placed in a folder, such as `input`.
In this folder, there are six files.
- `inpar.data`         -> list of filenames
- `FD_model.dat`       -> Parameters for FD modeling.
- `Source.dat`         -> Parameters for an earthquake. You can use the TauP software to compute the ray parameter of an event.
- `Receiver.dat`       -> Parameters for stations.
- `FK_model_left.dat`  -> 1D media on left side of the 2D model.
- `FK_model_right.dat` -> 1D media on right side of the 2D model.
- `Vp_model.su`        -> Vp velocity model in SU format.
- `Vs_model.su`        -> Vs velocity model in SU format.

The parameters are self-explaining. All physical quantities are international units.

# Utils
There are several auxiliary programs in `utils`.
- gauss_smooth.f90 -> A Fortran program to smooth the velocity model. Usually, the geological structure in real world is relatively smooth. You can use this program to smooth your velocity models. It provides two options to smooth, i.e., gaussian smooth and triangular smooth, the latter may has better structure-preserving characteristics. You can just smooth part of the model. You can see the usage by typing `./gauss_smooth`.
- readsu.m         -> A MATLAB script to read velocity model in SU format (i.e., the SEGY format after removing 3200 bytes textual header and 400 bytes binary header).
- writesu.m        -> A MATLAB script to write velocity model in SU format.
- makeRFitdecon.m  -> A MATLAB script to compute RFs using time-domain iterative deconvolution written by Bailey, et al. (2010).
- plot_PRF.m       -> A MATLAB script to plot P-wave receiver function.
- plot_SRF.m       -> A MATLAB script to plot S-wave receiver function.
- plot_seism.m     -> A MATLAB script to plot seismograms.
- plot_snap.m      -> A MATLAB script to plot snapshots.
- plot_model.m     -> A MATLAB script to plot velocity models.
- create_animation -> A MATLAB script to convert snapshots to animation.

# Examples


## Animation of wavefield snapshots
![Animation](https://github.com/YoushanLiu/FDFK2D/blob/master/images/Altyn.gif)


## Wavefield snapshots
![Snapshots](https://github.com/YoushanLiu/FDFK2D/blob/master/images/snapshots.png)


## P-wave Receiver Function
![PRF](https://github.com/YoushanLiu/FDFK2D/blob/master/images/PRF.png)


# License
FDFK2D is a free software, you can redistribute it and/or modify it under the terms of the MIT License. A copy of this license is provided in LICENSE.
