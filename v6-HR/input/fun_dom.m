
% THIS FUNCTION RETURNS longitudes/latitudes/depths of the domain

function [lon1,lat1,dep1,lon2,lat2,dep2] = fun_dom


%% DOMAIN 1: [85W - 36W] - [25N - 60N]

% Read bathymetric data file

addpath ~omarchal/Data/Bathymetry/SS97/;

fid = fopen('ss97_data_85W_36W_25N_60N.cgi','r');
C = textscan(fid,'%f %f %f');
fclose(fid);


% Assign longitudes/latitudes/depths

lonC  = C{1}-360;
latC  = C{2};
depC  = C{3};

ntot = length(lonC);


% Longitudes [deg]

      nlon1   =      1;
 lon1(nlon1) = lonC(1);

for i=1+1:ntot
    if latC(i) == latC(i-1)
             nlon1  = nlon1 + 1;
        lon1(nlon1) = lonC(i);
    else
        break
    end
end


% Latitudes [deg]

     nlat1  =      1;
lat1(nlat1) = latC(1);
   
for j=1+1:ntot
    if latC(j) == latC(j-1)
        continue
    else
             nlat1  = nlat1 + 1;
        lat1(nlat1) = latC(j);
    end
end


% Depths [m]

k = 0;

for j=1:nlat1
for i=1:nlon1
    k         = k + 1;
    dep1(i,j) = depC(k);
end
end


%% DOMAIN 2: [82W - 51W] - [31N - 50N]

% Longitudes [deg]

lon_w = -82; % -81
lon_e = -51;
nlon2 = 0;

for i=1:nlon1
    if lon1(i) > lon_w & lon1(i) < lon_e    
            nlon2  = nlon2 + 1;
       lon2(nlon2) =  lon1(i) ;
    end
end


% Latitudes [deg]

lat_s =  31;
lat_n =  50; % 46
nlat2 = 0;

for j=1:nlat1
    if lat1(j) > lat_s & lat1(j) < lat_n
            nlat2  = nlat2 + 1;
       lat2(nlat2) =  lat1(j) ;
    end
end


% Depths [m]

nlon2 = 0;
nlat2 = 0;

for i=1:nlon1
    
    if lon1(i) > lon_w & lon1(i) < lon_e
       nlon2 = nlon2 + 1;
       nlat2 = 0;
       for j=1:nlat1
           if lat1(j) > lat_s & lat1(j) < lat_n
                         nlat2  = nlat2 + 1;
              dep2(nlon2,nlat2) = dep1(i,j);
           end
       end
    end
    
end
          
