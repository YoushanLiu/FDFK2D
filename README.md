# FDFK-2D
Efficient two-dimensional teleseismic wavefiled hybrid simulation method for receiver function analysis

Until now, simulation of high-frerquency seismic wavefield in global-scale is still challenging even using large-scale computer clusters. One feasible way to mitigate this issue is to adopt a so-called hybrid method. This method usually computes the wavefield from source to the boundaries of the target heterogeneous media using a fast 1D wavefield sovler (such as AxiSEM, DSM, general Ray method, Normal mode, FK, etc.), which is referred to background wavefield. Then, the complex wavefield in the localized target heterogeneous media is computed using 2D/3D wavefield solver (such as finite difference method, pseduospectral method, finite element method, spectral element method, etc.) with the background wavefield as incident wavefield, which is referred to complete wavefield. In the boundaries of the target heterogenous model, a interface between background wavefield and complete wavefiled is applied.

This software is a hybrid method of the FD and FK to simulate teleseismic wavefield.

FD: Seismic wavefield modeling using collocated grid finite difference method in two-dimensional elastic isotropic meida.

FK: Seismic wavefield modeling using the frequency-wavenumber or propagator matrix method in three-dimensional elastic isotropic media.

If this software is used in your research, please consider cite one of the following papers.

Youshan Liu, Chenglong Wu, Tao Xu, Liang Zhao, Jiwen Teng. 2024. Efficient two-dimensional teleseismic wavefield hybrid simulation method for receiver function analysis. Seismological Research Letters, undergoing review.
Chenglong Wu, Tao Xu, Xiaobo Tian, Ross N. Mitchell, Jiyan Lin, Jianfeng Yang, Xin Wang, Zhanwu Lu. 2024. Underthrusting of Tarim lover crust beneath the Tibetan Plateau revealed by receiver function imaging. Geophysical Research Letters, 51, e2024GL108220.
