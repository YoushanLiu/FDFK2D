module FDFK_par

! relative path of the input data
character(256) :: prefix = ''

! output path
character(256) :: outpath = ''

! prefix of filename for seismograms
character(256) :: fsis = 'seis'

! filenames of FK, Vp and Vs, models
character(256) :: fFKmodel_left = '', fFKmodel_right = '', fVpmodel = '', fVsmodel = ''

! filenames of FDmodel, Receiver and Source
character(256) :: fFDmodel = '', freceiver = '', fsource = ''

! filename of parfile
character(256) :: finpar

! version information
character(16) :: ver = 'v0.9'


! latitude, longitude, and azimuth of the reference point of the profile (x0, z0)
real(8) :: lat_org = 0.0, lon_org = 0.0, az_org = 0.0


! rotation matrix and inverse rotation matrix with respect to (lat_org, lon_org, az)
real(8), dimension(3,3) :: rotate_matrix, inv_rotate_matrix


end module FDFK_par

