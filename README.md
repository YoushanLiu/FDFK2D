# FDFK2D
Efficient two-dimensional teleseismic wavefiled hybrid simulation method for receiver function analysis. This software is written in Fortran with highly optimization and enable parallel computing using the OpenMP.

# Background
Until now, simulation of high-frerquency seismic wavefield in global-scale is still challenging even using large-scale computer clusters. One feasible way to mitigate this issue is to adopt a so-called hybrid method. This method usually computes the wavefield from source to the boundaries of the target heterogeneous media using a fast 1D wavefield sovler (such as AxiSEM, DSM, general Ray method, Normal mode, FK, etc.), which is referred to background wavefield. Then, the complex wavefield in the localized target heterogeneous media is computed using 2D/3D wavefield solver (such as finite difference method, pseduospectral method, finite element method, spectral element method, etc.) with the background wavefield as incident wavefield, which is referred to complete wavefield. In the boundaries of the target heterogenous model, a interface between background wavefield and complete wavefiled is applied. This method can simulate the interactions of incoming teleseismic wavefield with local heterogeneities but reducing these heterogeneities to much smaller local computational domain, which can significantly reduce the computing costs of high-frequency teleseismic wavefield with less computational resources. 

# About
Purpose: This software is a hybrid method coupling the FD and FK methods to simulate teleseismic wavefield.

FD: Seismic wavefield modeling using collocated grid finite difference method in two-dimensional elastic isotropic meida.

FK: Seismic wavefield modeling using the frequency-wavenumber or propagator matrix method in three-dimensional elastic isotropic media.

# Citation
If you use this software for your researches and prepare publications, please citing at least one of the following papers:

- Youshan Liu, Chenglong Wu, Tao Xu, Liang Zhao, Jiwen Teng. 2024. Efficient two-dimensional teleseismic wavefield hybrid simulation method 
   for receiver function analysis. Seismological Research Letters, undergoing review.
- Chenglong Wu, Tao Xu, Xiaobo Tian, Ross N. Mitchell, Jiyan Lin, Jianfeng Yang, Xin Wang, Zhanwu Lu. 2024. Underthrusting of Tarim lover crust beneath the Tibetan Plateau revealed by receiver function imaging. Geophysical Research Letters, 51, e2024GL108220.

# Prerequisites
Before you compile this software, you should first install the following packages.
- blas
- lapack
- fftw3

# Installation and Usuage
Before you run this software, you can type `make lib` to copy the library files in the `lib` folder to `/usr/lib64` with root permission
or add the directory to the environment variable by setting `export LD_LIBRARY_PATH=/your_path:$LD_LIBRARY_PATH`.
Now, you can install this software as follows.
- `cd ./src`
- `make install`
- `cd ../Tests/test1`
- `./run_this_example`

# Documentation
All input parameters are placed in a folder, such as `input`.
In this folder, there are six files.
- `inpar.data`: list of Filenames
- FD_model.dat       -> Parameters for FD modeling.
- Source.dat         -> Parameters for a earthquake. You can use the TauP software to compute the ray parameter of a event.
- Receiver.dat       -> Parameters for stations.
- FK_model_left.dat  -> 1D media on left side of the 2D model.
- FK_model_right.dat -> 1D media on right side of the 2D model.
- Vp_model.su        -> Vp velocity model in SU format.
- Vs_model.su        -> Vs velocity model in SU format.

The parameters are self-explaining. All physical quantities are international units.

# Utils
There are several auxiliary programs in `utils`.
- gauss_smooth.f90 -> A Fortran program to smooth the velocity model. Usually, the geological structure in real world is relatively smooth. You can use this program to smooth your velocity models. It provides two options to smooth, i.e., gaussian smooth and triangular smooth, the latter may has better structure-preserving characteristics. You can just smooth part of the model. You can see the usage by typing `./gauss_smooth`.
- readsu.m         -> A MATLAB script to read velocity model in SU format (i.e., the SEGY format after removing 3200 bytes textual header and 400 bytes binary header).
- writesu.m        -> A MATLAB script to write velocity model in SU format.
- makeRFitdecon.m  -> A MATLAB script to compute RFs using time-domain iterative deconvolution written by Bailey, et al. (2010).
- plot_seism.m     -> A MATLAB script to plot seismograms.
- plot_snap.m      -> A MATLAB script to plot snapshots.

# License
FDFK2D is a free software: you can redistribute it and/or modify it under the terms of the MIT License. A copy of this license is provided in LICENSE.
