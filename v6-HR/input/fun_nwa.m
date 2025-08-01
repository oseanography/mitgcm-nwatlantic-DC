
%% THIS FUNCTION RETURNS longitudes/latitudes/depths/T/S of NWARC

function [lon_nwa,lat_nwa,dep_nwa,tem_nwa,sal_nwa] = fun_nwa

% Notes
% -----
% 1) NWARC refers to the Northwest Atlantic regional climatology version 2
% 2) The domain of this climatology is between 80W-40W and 32N-65N

addpath ../data/
cidT = 'nwa_decav_t00_10.nc';
cidS = 'nwa_decav_s00_10.nc'; 

% Read in ncdf T/S files

lon_nwa = ncread(cidT,'lon'  );
lat_nwa = ncread(cidT,'lat'  );
dep_nwa = ncread(cidT,'depth');
tem_nwa = ncread(cidT,'t_an' );
sal_nwa = ncread(cidS,'s_an' );