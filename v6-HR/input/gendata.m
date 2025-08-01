
% THIS MATLAB SCRIPT GENERATES THE INPUT DATA FOR /sean_v3/

clear all
close all


%% CONSTANTS

% Earth radius [m]

aa = 6370e3;

% Conversion factor [rad/deg]

rad = pi/180;


%% GET ETOPO DATA

disp(' ')
disp('Get ETOPO data ...')


% Get longitudes, latitudes, depths from ETOPO

addpath ../data/;

prec = 'real*8';
ieee = 'b';

fid         = fopen ('ETOPO_2022_v2_sean_lon.d','r'); 
lon_etopo   = fread (fid,prec,ieee)                ; 
              fclose(fid)                          ; 
n_lon_etopo = size(lon_etopo,1);
       
fid         = fopen ('ETOPO_2022_v2_sean_lat.d','r'); 
lat_etopo   = fread (fid,prec,ieee)                ; 
              fclose(fid)                          ; 
n_lat_etopo = size(lat_etopo,1);
       
fid        = fopen ('ETOPO_2022_v2_sean_dep.d','r')          ; 
dep_etopo  = fread (fid,[n_lon_etopo n_lat_etopo],prec,ieee); 
             fclose(fid)                                    ;
       

% Plot

figure(1)
contour(lon_etopo,lat_etopo,dep_etopo',[0 0],'k'); hold on
xlabel('LONGITUDE [\circW]'); ylabel('LATITUDE [\circN]');
set(gca,'XTick',-80:5:-55,'XTickLabel',80:-5:55);
axis([-82 -51 31 49])


%% GET NWA DATA & HORIZONTAL MODEL GRID

disp(' ')
disp('Get NWA data & horizontal model grid ...')


% Read longitudes/latitudes/depths/T/S of NWA v2 (0.1 deg)

[lon_nwa,lat_nwa,dep_nwa,tem_nwa,sal_nwa] = fun_nwa;


% Making NWA v2 depths negative

dep_nwa = -abs(dep_nwa);


% Number of longitudes/latitudes/depths in NWA v2 

n_lon_nwa = length(lon_nwa);
n_lat_nwa = length(lat_nwa);
n_dep_nwa = length(dep_nwa);


% Zonal/meridional spacings between cell faces [deg]

dxSpacing = 0.05; 
dySpacing = 0.05;


% Bounding longitude/latitude of model domain [deg]
%
% Notes
% -----
% lon_nwa(1) = longitude of westernmost  model grid cell center
% lat_nwa(1) = latitude  of southernmost model grid cell center

xgOrigin = lon_nwa(1) - dxSpacing/2;
ygOrigin = lat_nwa(1) - dySpacing/2;


% Number of grid points in zonal/meridional directions

nx = 560;
ny = 315; 


% Longitudes/latitudes of model grid cell centers [deg]

for i=1:nx
for j=1:ny    
    XC(i,j) = xgOrigin + dxSpacing/2 + (i-1)*dxSpacing;
    YC(i,j) = ygOrigin + dySpacing/2 + (j-1)*dySpacing;
end
end


% Plot SST from NWA

figure
contour(lon_etopo,lat_etopo,dep_etopo',[0 0],'k'); hold on
xlabel('LONGITUDE [\circW]'); ylabel('LATITUDE [\circN]');
set(gca,'XTick',-80:5:-55,'XTickLabel',80:-5:55); 
[CC,hh]=contour(lon_nwa,lat_nwa,tem_nwa(:,:,1)',[-10:1:50]); clabel(CC,hh) 
title('ORIGINAL NWA SST [\circC]')
axis([-82 -51 31 49])


% Plot SSS - 35 from NWA

figure
contour(lon_etopo,lat_etopo,dep_etopo',[0 0],'k'); hold on
xlabel('LONGITUDE [\circW]'); ylabel('LATITUDE [\circN]');
set(gca,'XTick',-80:5:-55,'XTickLabel',80:-5:55); 
[CC,hh]=contour(lon_nwa,lat_nwa,sal_nwa(:,:,1)'-35,[-5:0.5:-0.5]); clabel(CC,hh) 
[CC,hh]=contour(lon_nwa,lat_nwa,sal_nwa(:,:,1)'-35,[0.5:0.5:2]); clabel(CC,hh) 
[CC,hh]=contour(lon_nwa,lat_nwa,sal_nwa(:,:,1)'-35,[2:0.2:5]); clabel(CC,hh) 
title('ORIGINAL NWA SSS - 35')
axis([-82 -51 31 49])
colormap(bluewhitered)


%% MODEL TOPOGRAPHY

disp(' ')
disp('Set bottom topography ...')


% Depths of model grid cell centers [m]

h(1:nx,1:ny) = 0;

for i1=1:nx
for j1=1:ny
    
    lon_cell_w = XC(i1,j1) - dxSpacing/2;
    lon_cell_e = XC(i1,j1) + dxSpacing/2;
    
    for i2=1:n_lon_etopo
    
        if lon_etopo(i2) >= lon_cell_w & lon_etopo(i2) <= lon_cell_e
    
           lat_cell_s = YC(i1,j1) - dySpacing/2;
           lat_cell_n = YC(i1,j1) + dySpacing/2;

            kk(i1,j1) = 0;
            h (i1,j1) = 0;

            for j2=1:n_lat_etopo
                if lat_etopo(j2) >= lat_cell_s & lat_etopo(j2) <= lat_cell_n
                   kk(i1,j1) = kk(i1,j1) + 1;
                   h (i1,j1) = h (i1,j1) + dep_etopo(i2,j2);
                end
            end
            
        end
        
    end

    if kk(i1,j1) > 0 
        h(i1,j1) = h(i1,j1) / kk(i1,j1);
    else
        disp(' ')
        disp('STOP: no h data in this grid cell')
        disp(' ')
        return
    end

end
end


% Set h(,) to 0 at grid cells where h(,) >= 0

h_wall = 0;

for i=1:nx
for j=1:ny
    if h(i,j) > h_wall
       h(i,j) = 0;
    end
end
end


% Dry grid cells in Gulf of St Lawrence

for i=1:nx
for j=1:ny
    
    if XC(i,j) <= -60.55 & YC(i,j) >= 46.00
        h(i,j) = 0;
    end

    if XC(i,j) >= -63.95 & XC(i,j) <= -61.50 & YC(i,j) >= 45.65
        h(i,j) = 0;
    end
   
end
end


% Dry grid cells in continental depression

for i=1:nx
for j=1:ny
    
    if XC(i,j) > -80 & XC(i,j) < -75 & ...
       YC(i,j) >  43 & YC(i,j) <  44
        h(i,j) = 0;
    end
    
end
end


% Dry all grid cells along western bry at i = 1

h(1,1:ny) = 0;


% Dry all grid cells along northern bry at j = ny

h(1:nx,ny) = 0; 


% Plot wet model grid points: All wet pts

xp(1:nx,1:ny) = NaN; 
yp(1:nx,1:ny) = NaN;
for i=1:1:nx
for j=1:1:ny
    if h(i,j) ~= 0
       xp(i,j) = XC(i,j);
       yp(i,j) = YC(i,j);
    end
end
end

figure(1)
plot(xp,yp,'.','Color',[0.85 0.85 0.85]); 


% Plot other isobaths

contour(lon_etopo,lat_etopo,dep_etopo',[0 -200 -1000 -2000 -3000 -4000 -5000],'k'); hold on


% Print names of geographic/topographic features

text(-78.70,33.00,'SAB','FontSize',8,'Rotation',45)
text(-78.5,35.50,'CH' ,'FontSize' ,8); 
text(-74.00,39.00,'MAB','FontSize',8,'Rotation',45)
text(-70.00,40.75,'NS' ,'FontSize',8)
text(-68.00,41.40,'GB' ,'FontSize',8)
text(-68.65,43.00,'GoM','FontSize',8)
text(-61.50,43.85,'SS' ,'FontSize',8)
text(-54.50,45.75,'NfS','FontSize',8)


% Write bottom topography

ieee = 'b';
prec = 'real*8';
fid  = fopen('topog.d','w',ieee); fwrite(fid,h,prec); fclose(fid);


%% VERTICAL MODEL GRID (interface centered approach)

disp(' ')
disp('Set vertical model grid ...')


% Select grid

igrid = 3;


% Grid 0

if igrid == 0
    
    % Number & depths of vertical levels (grid cell centers) [m]

      nz  =    1;
    z(nz) = -2.5;

    for k=1+1:n_dep_nwa
          nz  = nz + 1;
        z(nz) = dep_nwa(k);
    end

    % Vertical spacing between vertical levels [m] (=delRc in input/data)

    delRc(1) = 0 - z(1);

    for k=1+1:nz
        delRc(k) = z(k-1) - z(k);
    end
    
end



% Grid 1

if igrid == 1
    
    % Number & depths of levels (grid cell centers) for all z <= z0 [m]

    z0     =-1000;
      nz0  =    1;
    z(nz0) = -2.5;

    for k=1+1:n_dep_nwa
        if dep_nwa(k) >= z0
              nz0  = nz0 + 1;
            z(nz0) = dep_nwa(k);
        end
    end

    % Spacing between levels for all z <= z0 [m] (=delRc in input/data)

    delRc(1) = 0 - z(1);

    for k=1+1:nz0
        delRc(k) = z(k-1) - z(k);
    end


    % Number & depths of levels (grid cell centers) for all z > z0 [m]

    rr = 0.99;
    nz = nz0 ;
    while z(nz) >= dep_nwa(n_dep_nwa)
        nz        = nz + 1               ;
        delRc(nz) = rr      * delRc(nz-1);
        z    (nz) = z(nz-1) - delRc(nz  );
    end
    dum = delRc; clear delRc; delRc = dum(1:nz-1);
    dum = z    ; clear z    ; z     = dum(1:nz-1); nz = nz - 1;
    
end


% Grid 2

if igrid == 2
    
   nz       =   300   ;
   mean_z0  = -1000   ;
   var_z0   =  (1000/1)^2   ;
   delRc_fac   =  10000*4  ; % [m]
   z0   (1) = 0;
   delRc0   =  5*2;
   dz0      = dep_nwa(end)/(nz+1);

   for k=1:nz
       z0   (k) = k*dz0; 
       delRc(k) = delRc0 ...
                + delRc_fac/sqrt(2*pi*var_z0) * exp(-((z0(k)-mean_z0)^2)/(2*var_z0));
   end

%    figure
%    plot(delRc,z0,'+-'); xlabel('\Delta z_k [m]'); ylabel('z_0 [m]')
   
   z(1) = dep_nwa(1+1);
   for k=1+1:nz
       z(k) = z(k-1) - delRc(k);
   end
 
%    figure
%    plot(delRc(1:nz),z(1:nz),'+-'); xlabel('\Delta z_k [m]'); ylabel('z [m]')   
    
end


% Grid 3

if igrid ==3

    nz          = -dep_nwa(end)/10;
    delRc(1:nz) = -dep_nwa(end)/nz;

    z(1) = -delRc(1);
    for k=1:nz
        z(k) = -k*delRc(k);
    end

end

% Plot

figure
plot(z,'+-'); xlabel('k'); ylabel('z_k [m]');
figure
plot(delRc,z,'+-'); xlabel('\Delta z_k'); ylabel('z_k [m]');


% Increment delRc by 1 (error message if delRc(nz+1) not in input/data)

delRc(nz+1) = delRc(nz);


% Print out delRc

disp('   delRc to put in /input/data:')
fprintf('   delRc='); fprintf(' %4.3g,',delRc); fprintf('\n');
disp(' ')


% Write delRc 

ieee = 'b';
prec = 'real*8';
fid  = fopen('delRc.d','w',ieee); fwrite(fid,delRc,prec); fclose(fid);


% Vertical spacing between cell interfaces [m] 

z_w(1) = 0;
for k=1+1:nz
    z_w(k) = (z(k-1)+z(k))/2;
end
z_w(nz+1) = z(nz) - (z_w(nz)-z(nz));

for k=1:nz
    dz(k) = z_w(k) - z_w(k+1);
end


% Write vertical levels

ieee = 'b';
prec = 'real*8';
fid  = fopen('z.d','w',ieee); fwrite(fid,z,prec); fclose(fid);


%% VERTICAL INTERPOLATION OF NWA T,S ONTO VERTICAL MODEL GRID

disp('Vertical interpolation of NWA T,S onto vertical model grid ...')


% Interpolate

for k1=1:nz
    
    % Get k-index of NWA pt just above z(k1)

    dz_up = aa;
    for k2=1:n_dep_nwa
        if     abs(dep_nwa(k2)-z(k1)) < dz_up & ...
                   dep_nwa(k2)>z(k1)
           dz_up = dep_nwa(k2)-z(k1);
           k_up  = k2;
        end
    end

    % Get k-index of NWA pt just below z(k1)

    dz_do = aa;
    for k2=1:n_dep_nwa
        if     abs(z(k1)-dep_nwa(k2)) < dz_do & ...
                   z(k1)>dep_nwa(k2)
           dz_do = z(k1)-dep_nwa(k2);
           k_do  = k2;
        end
    end       
    
    % Interpolate
    
    tem_nwa_int(:,:,k1) = (dz_do*tem_nwa(:,:,k_up)  + ...
                           dz_up*tem_nwa(:,:,k_do)) / (dz_do+dz_up);
                       
    sal_nwa_int(:,:,k1) = (dz_do*sal_nwa(:,:,k_up)  + ...
                           dz_up*sal_nwa(:,:,k_do)) / (dz_do+dz_up);
   
end

    
% Clear original NWA T,S values
    
clear tem_nwa sal_nwa; 
             

% Plot vertically interpolated T from NWA

figure
contour(lon_etopo,lat_etopo,dep_etopo',[0 0],'k'); hold on
xlabel('LONGITUDE [\circW]'); ylabel('LATITUDE [\circN]');
set(gca,'XTick',-80:5:-55,'XTickLabel',80:-5:55); 
[CC,hh]=contour(lon_nwa,lat_nwa,tem_nwa_int(:,:,1)',[-10:1:50]); clabel(CC,hh) 
title('z-INTERPOLATED NWA T [\circC]')
axis([-82 -51 31 49])


% Plot vertically interpolated S-35 from NWA

figure
contour(lon_etopo,lat_etopo,dep_etopo',[0 0],'k'); hold on
xlabel('LONGITUDE [\circW]'); ylabel('LATITUDE [\circN]');
set(gca,'XTick',-80:5:-55,'XTickLabel',80:-5:55); 
[CC,hh]=contour(lon_nwa,lat_nwa,sal_nwa_int(:,:,1)'-35,[-5:0.5:-0.5]); clabel(CC,hh) 
[CC,hh]=contour(lon_nwa,lat_nwa,sal_nwa_int(:,:,1)'-35,[0.5:0.5:2]); clabel(CC,hh) 
[CC,hh]=contour(lon_nwa,lat_nwa,sal_nwa_int(:,:,1)'-35,[2:0.2:5]); clabel(CC,hh) 
title('z-INTERPOLATED NWA S - 35')
axis([-82 -51 31 49])
colormap(bluewhitered)


%% ATTRIBUTES OF GULF STREAM INFLOW (GSI)

disp(' ')
disp('Set attributes of GSI ...')


% Volume transport of GSI [m3/s]

V_gsi = 45e6; % CHG            


% Indices of westernmost/easternmost grid pts in GSI

for i=1:nx
    if XC(i,1) > -77.6-0.5 % CHG
       iw_gsi = i;
       break
    end
end

for i=1:nx
    if XC(i,1) < -76.6-0.5 % CHG
       ie_gsi = i;
    end
end


% Print grid pts in GSI

figure(1)
plot(XC(iw_gsi,1),YC(iw_gsi,1),'^','MarkerSize',5','Color','b');
plot(XC(ie_gsi,1),YC(ie_gsi,1),'^','MarkerSize',5','Color','b');


%% ATTRIBUTES OF GULF STREAM OUTFLOW (GSO)

disp(' ')
disp('Set attributes of GSO ...')


% Volume transport of GSO [m3/s]

U_gso = (95e6)*1; % CHG


% Indices of westernmost/easternmost grid pts in GSO (Hoog and Johns 1995)

for j=1:ny
    if YC(nx,j) > 38.5 % 39
       js_gso = j;
       break
    end
end

for j=1:ny
    if YC(nx,j) < 40.5 % 40
       jn_gso = j;
    end
end


% Print coordinates of GSO

% figure(1)
% plot(XC(nx,js_gso),YC(nx,js_gso),'>','MarkerSize',5,'Color','r');
% plot(XC(nx,jn_gso),YC(nx,jn_gso),'>','MarkerSize',5,'Color','r');


%% ATTRIBUTES OF SOUTHERN COMPENSATING INFLOW (SCI)

disp(' ')
disp('Set attributes of SCI ...')


% Volume transport of SCI [m3/s]

U_sci = (-22.5e6)*1; % CHG


% Indices of southern/northern grid pts in SCI

js_sci = 1;
jn_sci = js_gso - 1;


% Print coordinates of SCI

% figure(1)
% plot(XC(nx,js_sci),YC(nx,js_sci),'<','MarkerSize',5,'Color','g');
% plot(XC(nx,jn_sci),YC(nx,jn_sci),'<','MarkerSize',5,'Color','g'); 


%% ATTRIBUTES OF NORTHERN COMPENSATING INFLOW (NCI) OR DWBC

disp(' ')
disp('Set attributes of NCI ...')


% Volume transport of NCI [m3/s]

U_nci = -30e6; % CHG


% Indices of southern/northern grid pts in DWBC 

% for j=ny:-1:1
%     if h(nx,j) < -200
%         jn_nci = j;
%         break
%     end
% end
% js_nci = jn_nci - nearest(1/dySpacing); % 1 degree of latitude

% Set to a region between 1000 m and 4500 m depths

for j=ny:-1:1
    if h(nx,j) < -1000 
        jn_nci = j;
        break
    end
end

for j=ny:-1:1
    if h(nx,j) < -4500 
        js_nci = j-1;
        break
    end
end


% Print coordinates of DWBC

figure(1)
plot(XC(nx,js_nci),YC(nx,js_nci),'<','MarkerSize',5,'Color','r');
plot(XC(nx,jn_nci),YC(nx,jn_nci),'<','MarkerSize',5,'Color','r');


%% INITIAL CONDITIONS FOR T,S

disp(' ')
disp('Get initial T,S ...')


% Get i-index of NWA grid pts closest to model grid pts

for i=1:nx
    lon_mod(1:n_lon_nwa) = XC(i,1);
    [~,i_nwa(i)]         = min(abs(lon_nwa-lon_mod'));
end


% Get j-index of NWA grid pts closest to model grid pts
    
for j=1:ny
    lat_mod(1:n_lat_nwa) = YC(1,j);
    [~,j_nwa(j)]         = min(abs(lat_nwa-lat_mod'));
end
    

% Initial T,S at model grid pts with/without missing T,S in NWA

for i=1:nx
for j=1:ny
    t(i,j,:) = tem_nwa_int(i_nwa(i),j_nwa(j),:);
    s(i,j,:) = sal_nwa_int(i_nwa(i),j_nwa(j),:);
end
end


% Initial T at model grid pts with missing T in NWA

disp(' ')

for k=1:nz
    
    % Initialize number of missing T values at level k
    
    n_mis_t(k) = 0;
    
    % Loop over x,y
    
    for i1=1:nx
    for j1=1:ny
        
        % Test that grid pt is wet and has missing T
        
        if z(k) > h(i1,j1) & isnan(tem_nwa_int(i_nwa(i1),j_nwa(j1),k))==1

                        n_mis_t(k)  = n_mis_t(k) + 1;
            lon_mis_t(k,n_mis_t(k)) = XC(i1,j1);
            lat_mis_t(k,n_mis_t(k)) = YC(i1,j1);
            
            % Go south to find nearest available NWA T
            
            ts =       NaN;
            js = j_nwa(j1);
            while js > 1 & isnan(ts)==1
                js = js - 1;
                ts = tem_nwa_int(i_nwa(i1),js,k);
            end
            
            % Go north to find nearest available NWA T
            
            tn =       NaN;
            jn = j_nwa(j1);
            while jn < n_lat_nwa & isnan(tn)==1
                jn = jn + 1;
                tn = tem_nwa_int(i_nwa(i1),jn,k);
            end           
            
            % Go west to find nearest available NWA T
            
            tw =       NaN;
            iw = i_nwa(i1);
            while iw > 1 & isnan(tw)==1
                iw = iw - 1;
                tw = tem_nwa_int(iw,j_nwa(j1),k);
            end                   
            
            % Go east to find nearest available NWA T
            
            te =       NaN;
            ie = i_nwa(i1);
            while ie < n_lon_nwa & isnan(te)==1
                ie = ie + 1;
                te = tem_nwa_int(ie,j_nwa(j1),k);
            end                   
            
            % Fill in missing T with average of ts,tn,tw,te
            
            ifs = 1-isnan(ts); ts(isnan(ts)==1) = 0;
            ifn = 1-isnan(tn); tn(isnan(tn)==1) = 0;
            ifw = 1-isnan(tw); tw(isnan(tw)==1) = 0;
            ife = 1-isnan(te); te(isnan(te)==1) = 0;
            
            t(i1,j1,k) = (ifs*ts+ifn*tn+ifw*tw+ife*te) / ...
                         (ifs   +ifn   +ifw   +ife   );      
                     
            % Print on screen
            
%             disp(['missing T: k=',num2str(k), ...
%                   ' - n_mis='    ,num2str(n_mis_t(k)), ...
%                   ' - ts='       ,num2str(ts,'%.2f'), ...
%                   ' - tn='       ,num2str(tn,'%.2f'), ...
%                   ' - tw='       ,num2str(tw,'%.2f'), ...
%                   ' - te='       ,num2str(te,'%.2f')])
            
        end           
        
    end
    end
    
    % Print on screen
    
    disp(['   missing T: k=',num2str(k),' - n_mis=',num2str(n_mis_t(k))])
    
end


% Initial S at model grid pts with missing S in NWA

disp(' ')

for k=1:nz
    
    % Initialize number of missing S values at level k
    
    n_mis_s(k) = 0;
    
    % Loop over x,y
    
    for i1=1:nx
    for j1=1:ny
        
        % Test that grid pt is wet and has missing S
        
        if z(k) > h(i1,j1) & isnan(sal_nwa_int(i_nwa(i1),j_nwa(j1),k))==1

                        n_mis_s(k)  = n_mis_s(k) + 1;
            lon_mis_s(k,n_mis_s(k)) = XC(i1,j1);
            lat_mis_s(k,n_mis_s(k)) = YC(i1,j1);
            
            % Go south to find nearest available NWA S
            
            ss =       NaN;
            js = j_nwa(j1);
            while js > 1 & isnan(ss)==1
                js = js - 1;
                ss = sal_nwa_int(i_nwa(i1),js,k);
            end
            
            % Go north to find nearest available NWA S
            
            sn =       NaN;
            jn = j_nwa(j1);
            while jn < n_lat_nwa & isnan(sn)==1
                jn = jn + 1;
                sn = sal_nwa_int(i_nwa(i1),jn,k);
            end           
            
            % Go west to find nearest available NWA S
            
            sw =       NaN;
            iw = i_nwa(i1);
            while iw > 1 & isnan(sw)==1
                iw = iw - 1;
                sw = sal_nwa_int(iw,j_nwa(j1),k);
            end                   
            
            % Go east to find nearest available NWA S
            
            se =       NaN;
            ie = i_nwa(i1);
            while ie < n_lon_nwa & isnan(se)==1
                ie = ie + 1;
                se = sal_nwa_int(ie,j_nwa(j1),k);
            end                   
            
            % Fill in missing S with average of ss,sn,sw,se
            
            ifs = 1-isnan(ss); ss(isnan(ss)==1) = 0;
            ifn = 1-isnan(sn); sn(isnan(sn)==1) = 0;
            ifw = 1-isnan(sw); sw(isnan(sw)==1) = 0;
            ife = 1-isnan(se); se(isnan(se)==1) = 0;
            
            s(i1,j1,k) = (ifs*ss+ifn*sn+ifw*sw+ife*se) / ...
                         (ifs   +ifn   +ifw   +ife   );     
                     
            % Print on screen
            
%             disp(['missing S: k=',num2str(k), ...
%                   ' - n_mis='    ,num2str(n_mis_s(k)), ...
%                   ' - ss='       ,num2str(ss,'%.2f'), ...
%                   ' - sn='       ,num2str(sn,'%.2f'), ...
%                   ' - sw='       ,num2str(sw,'%.2f'), ...
%                   ' - se='       ,num2str(se,'%.2f')])
            
        end
                    
    end
    end
    
    disp(['   missing S: k=',num2str(k),' - n_mis=',num2str(n_mis_s(k))])
    
end


% Plot locations of missing T at selected level ksel

figure; ksel = 10;
contour(lon_etopo,lat_etopo,dep_etopo',[0 0],'k'); hold on
xlabel('LONGITUDE [\circW]'); ylabel('LATITUDE [\circN]');
set(gca,'XTick',-80:5:-55,'XTickLabel',80:-5:55); 
axis([-82 -51 31 49])

clear xp yp
for i=1:n_mis_t(ksel)
    xp(i)=lon_mis_t(ksel,i);
    yp(i)=lat_mis_t(ksel,i);
end
if n_mis_t(ksel) > 0
    plot(xp,yp,'*','MarkerSize',5,'Color',[0 0 1])
    title('MISSING T VALUES')
end


% Plot locations of missing S at selected level ksel

figure; ksel = 10;
contour(lon_etopo,lat_etopo,dep_etopo',[0 0],'k'); hold on
xlabel('LONGITUDE [\circW]'); ylabel('LATITUDE [\circN]');
set(gca,'XTick',-80:5:-55,'XTickLabel',80:-5:55); 
axis([-82 -51 31 49])

clear xp yp
for i=1:n_mis_s(ksel)
    xp(i)=lon_mis_s(ksel,i);
    yp(i)=lat_mis_s(ksel,i);
end
if n_mis_s(ksel) > 0
    plot(xp,yp,'*','MarkerSize',5,'Color',[0 0 1])
    title('MISSING S VALUES')
end


% Total number of wet grid pts

ntot = 0;

for i=1:nx
for j=1:ny
for k=1:nz
    if z(k) > h(i,j)
        ntot = ntot + 1;
    end
end
end
end


% Print on screen fraction of missing T,S values

disp(' ')
disp(['   Fraction of missing T values = ',num2str(sum(n_mis_t)/ntot)])
disp(['   Fraction of missing S values = ',num2str(sum(n_mis_s)/ntot)])


%% REFERENCE T,S PROFILES tRef,sRef
%
% Note
% ----
% Reference T,S profiles are domain-averaged profiles from NWA (0.1 deg)

disp(' ')
disp('Get tRef,sRef ...')


% tRef profile [deg]

for k=1:nz
    
    tRef(k) = 0;
    surf    = 0;
    
    for i=1:nx
    for j=1:ny

        if z(k) > h(i,j) & isnan(t(i,j,k))==0

           phi_s   = (YC(i,j)-dySpacing/2)*rad;
           phi_n   = (YC(i,j)+dySpacing/2)*rad;

           dsurf   = (sin(phi_n)-sin(phi_s)) * (dxSpacing*rad);

           tRef(k) = tRef(k) + t(i,j,k)*dsurf;
           surf    = surf    +          dsurf;
           
        end

    end
    end
    
    tRef(k) = tRef(k) / surf;

end


% sRef profile

for k=1:nz
    
    sRef(k) = 0;
    surf    = 0;
    
    for i=1:nx
    for j=1:ny

        if z(k) > h(i,j) & isnan(s(i,j,k))==0

           phi_s   = (YC(i,j)-dySpacing/2)*rad;
           phi_n   = (YC(i,j)+dySpacing/2)*rad;

           dsurf   = (sin(phi_n)-sin(phi_s)) * (dxSpacing*rad);

           sRef(k) = sRef(k) + s(i,j,k)*dsurf;
           surf    = surf    +          dsurf;
           
        end

    end
    end
    
    sRef(k) = sRef(k) / surf;

end


% Plot tRef,sRef profiles

figure
plot(tRef,z,'b+-'); xlabel('TEMPERATURE [\circC]'); ylabel('DEPTH [m]')
figure
plot(sRef,z,'r+-'); xlabel('SALINITY'); ylabel('DEPTH [m]')


% Print tRef,sRef profiles

disp('   tRef/sRef to put in /input/data:')
fprintf('   tRef='); fprintf(' %8.6g,',tRef); fprintf('\n');
fprintf('   sRef='); fprintf(' %8.6g,',sRef); fprintf('\n');


% Write tRef,sRef profiles

ieee = 'b';
prec = 'real*8';
fid  = fopen('tRef.d','w',ieee); fwrite(fid,tRef,prec); fclose(fid);
fid  = fopen('sRef.d','w',ieee); fwrite(fid,sRef,prec); fclose(fid);


%% REPLACE NaN VALUES IN INITIAL T,S FIELDS WITH tRef,sRef
% 
% Note
% ----
% This is done to prevent MITgcm to handle NaN values

disp(' ')
disp('Replace NaN values in initial T,S fields with tRef,sRef ...')


% Replace NaN T values

for i=1:nx
for j=1:ny
for k=1:nz
    
    if isnan(t(i,j,k))==1
             t(i,j,k) = tRef(k);
    end

end
end
end


% Replace NaN S values

for i=1:nx
for j=1:ny
for k=1:nz
    
    if isnan(s(i,j,k))==1
             s(i,j,k) = sRef(k);
    end

end
end
end


% Plot initial SST field

clear tp; tp(1:nx,1:ny) = NaN;
for i=1:nx
for j=1:ny
    if h(i,j) ~= 0
        tp(i,j) = t(i,j,1);
    end
end
end

figure
contour(lon_etopo,lat_etopo,dep_etopo',[0 0],'k'); hold on
xlabel('LONGITUDE [\circW]'); ylabel('LATITUDE [\circN]');
set(gca,'XTick',-80:5:-55,'XTickLabel',80:-5:55);
axis([-82 -51 31 49])
clear CC hh
[CC,hh] = contour(XC,YC,tp,[-10:1:40]); clabel(CC,hh,'FontSize',7);
title('INITIAL SST [\circC]')


% Plot initial SSS field

clear sp; sp(1:nx,1:ny) = NaN;
for i=1:nx
for j=1:ny
    if h(i,j) ~= 0
        sp(i,j) = s(i,j,1);
    end
end
end

figure
contour(lon_etopo,lat_etopo,dep_etopo',[0 0],'k'); hold on
xlabel('LONGITUDE [\circW]'); ylabel('LATITUDE [\circN]');
set(gca,'XTick',-80:5:-55,'XTickLabel',80:-5:55);
axis([-82 -51 31 49])
clear CC hh
[CC,hh] = contour(XC,YC,sp-35,[-5:0.5:-0.5]); clabel(CC,hh,'FontSize',7);
clear CC hh
[CC,hh] = contour(XC,YC,sp-35,[0.5:0.5:2]); clabel(CC,hh,'FontSize',7);
clear CC hh
[CC,hh] = contour(XC,YC,sp-35,[2:0.2:5]); clabel(CC,hh,'FontSize',7);
colormap(bluewhitered)
title('INITIAL SSS - 35')


% Write in files

fid = fopen('t.init','w',ieee); fwrite(fid,t,prec); fclose(fid);
fid = fopen('s.init','w',ieee); fwrite(fid,s,prec); fclose(fid);


%% OPEN BOUNDARY CONDITIONS: T,S

disp(' ')
disp('Set boundary conditions for T,S ...')


% T,S at western & eastern bries

for j=1:ny
for k=1:nz
    
    tWest(j,k) = t( 1,j,k);
    tEast(j,k) = t(nx,j,k);
    
    sWest(j,k) = s( 1,j,k);
    sEast(j,k) = s(nx,j,k);
        
end
end


% T,S at southern bry

for i=1:nx
for k=1:nz
        
    tSouth(i,k) = t(i,1,k);
    sSouth(i,k) = s(i,1,k);
      
end
end


% T,S at northern bry

for i=1:nx
for k=1:nz
    
    tNorth(i,k) = t(i,ny,k);
    sNorth(i,k) = s(i,ny,k);
    
end
end


%% OPEN BOUNDARY CONDITIONS: u,v (default)

disp(' ')
disp('Set boundary conditions for u,v: Default ...')


% u,v at all bries: Default values

uWest (1:ny,1:nz) = 0; uEast (1:ny,1:nz) = 0; uBaroEast(1:ny,1:nz) = 0;
vWest (1:ny,1:nz) = 0; vEast (1:ny,1:nz) = 0;

uSouth(1:nx,1:nz) = 0; uNorth(1:nx,1:nz) = 0;
vSouth(1:nx,1:nz) = 0; vNorth(1:nx,1:nz) = 0;


%% OPEN BOUNDARY CONDITIONS: v at GULF STREAM INFLOW

disp(' ')
disp('Set boundary conditions for u,v: v at GSI ...')


% Surface area of GSI [m2]

dz_gsi(1:nx,1:nz) = 0;
in_gsi(1:nx,1:nz) = 0;
surf              = 0;
acosphidl         = aa * cos(ygOrigin*rad) * (dxSpacing*rad);

for i=iw_gsi:ie_gsi
for k=1     :nz
              
    if     h(i,1) > z_w(k+1) & h(i,1) < z_w(k)
      dz_gsi(i,k) = z_w(k  ) - h(i,1);
      in_gsi(i,k) = 1;
    elseif h(i,1) < z_w(k+1)
      dz_gsi(i,k) = dz(k);
      in_gsi(i,k) = 1;
    end  

    surf = surf + acosphidl*dz_gsi(i,k);
    
end
end 


% Specify only barotropic velocity [m/s]

for i=iw_gsi:ie_gsi
for k=1     :nz
    vSouth(i,k) = (V_gsi/surf) * in_gsi(i,k); % CHG
end
end


%% OPEN BOUNDARY CONDITIONS: u at GULF STREAM OUTFLOW

disp(' ')
disp('Set boundary conditions for u,v: u at GSO ...')


% Surface area of GSO [m2]

surf = 0;

for j=js_gso:jn_gso
for k=1:nz
    if z(k) > h(nx,j) 
        surf = surf + aa*(dySpacing*rad)*dz(k);
    end
end
end


% Specify only barotropic velocity [m/s]

for j=js_gso:jn_gso
for k=1:nz
    if z(k) > h(nx,j) 
%         uEast(j,k) = U_gso/surf; % CHG
    end
end
end


%% OPEN BOUNDARY CONDITIONS: u at SOUTHERN COMPENSATING INFLOW

disp(' ')
disp('Set boundary conditions for u,v: u at SCI ...')


% Surface area of SCI [m2]

surf = 0;

for j=js_sci:jn_sci
for k=1:nz
    if z(k) > h(nx,j) 
       surf = surf + aa*(dySpacing*rad)*dz(k);
    end
end
end


% Specify only barotropic velocity [m/s]

for j=js_sci:jn_sci
for k=1:nz
    if z(k) > h(nx,j) 
%         uEast(j,k) = U_sci/surf; % CHG
    end
end
end


%% OPEN BOUNDARY CONDITIONS: u and v at NORTHERN COMPENSATING INFLOW OR DWBC

disp(' ')
disp('Set boundary conditions for u,v: u at NCI ...')


% Surface area of DWBC [m2]
% restrict DWBC to below 1000 m isobath 
nz_dwbc = 0;
for k=1:nz
   if z(k)<-1000
        nz_dwbc = k;
        break
    end
end

dz_nci(1:ny,1:nz) = 0;
in_nci(1:ny,1:nz) = 0;
surf              = 0;
adphi             = aa * (dySpacing*rad);

for j=js_nci:jn_nci
for k=nz_dwbc:nz
    
    if h(nx,j) > z_w(k+1) & h(nx,j) < z_w(k)
      dz_nci(j,k) = z_w(k   ) - h(nx,j);
      in_nci(j,k) = 1;
    elseif h(nx,j) < z_w(k+1)
      dz_nci(j,k) = dz(k);
      in_nci(j,k) = 1;
    end  

    surf = surf + adphi*dz_nci(j,k);
    
end
end 

% Specify only barotropic velocity [m/s]
% Calculate barotropic velocity and assign velocity vector components to
% create an along-isobath inflow for DWBC

m = -0.778; % isobath directional slope (y=mx+b)
angle = atan(m);

for j=js_nci:jn_nci
for k=nz_dwbc:nz
    uBaroEast(j,k) = U_nci/surf * in_nci(j,k); % CHG
    uEast(j,k) = uBaroEast(j,k) * cos(angle);
    vEast(j,k) = uBaroEast(j,k) * sin(angle);
end
end


%% CHECK: VOLUME TRANSPORTS ACROSS OPEN BOUNDARIES

disp(' ')
disp('Compute transports across open bries ...')


% Transport across GSI [m3/s]

V_gsi     = 0;
acosphidl = aa * cos(ygOrigin*rad) * (dxSpacing*rad);

for i=iw_gsi:ie_gsi
for k=1:nz
    V_gsi = V_gsi + vSouth(i,k)*acosphidl*dz_gsi(i,k);
end
end

disp(['   Transport of GSI = ',num2str(V_gsi*1e-6),' Sv']);


% Transport across GSO [m3/s]

% U_gso = 0;
% adphi = aa * (dySpacing*rad);
% 
% for j=js_gso:jn_gso
% for k=1:nz
%     U_gso = U_gso + uEast(j,k)*adphi*dz(k);
% end
% end
% 
% disp(['   Transport of GSO = ',num2str(U_gso*1e-6),' Sv']);


% Transport across SCI [m3/s]

% U_sci = 0;
% adphi = aa * (dySpacing*rad);
% 
% for j=js_sci:jn_sci
% for k=1:nz
%     U_sci = U_sci + uEast(j,k)*adphi*dz(k);
% end
% end
% 
% disp(['   Transport of SCI = ',num2str(U_sci*1e-6),' Sv']);


% Transport across DWBC [m3/s]

U_nci = 0;
adphi = aa * (dySpacing*rad);

for j=js_nci:jn_nci
for k=nz_dwbc:nz
    U_nci = U_nci + uBaroEast(j,k)*adphi*dz_nci(j,k);
end
end

disp(['   Transport of DWBC = ',num2str(U_nci*1e-6),' Sv']);


%% PLOT SECTIONS ALONG SOUTHERN BRY

disp(' ')
disp('Plot sections along southern bry ...')

clear xp zp

% T section

xp(1:nx,1:nz) = NaN;
zp(1:nx,1:nz) = NaN;
cp(1:nx,1:nz) = NaN;

for i=1:nx
for k=1:nz
    if (z(k) >= h(i,1)) 
        xp(i,k) = XC(i,1);
        zp(i,k) = z (k  );
        cp(i,k) = tSouth(i,k);
    end
end
end

figure
[CC,hh] = contour(xp,zp,cp); clabel(CC,hh); hold on;
plot(XC(:,1),h(:,1),'k+-')
title('T SECTION ALONG SOUTHERN BRY')
xlabel('LONGITUDE [\circW]'); ylabel('DEPTH [m]');


% S section

xp(1:nx,1:nz) = NaN;
zp(1:nx,1:nz) = NaN;
cp(1:nx,1:nz) = NaN;

for i=1:nx
for k=1:nz
    if (z(k) >= h(i,1)) 
        xp(i,k) = XC(i,1);
        zp(i,k) = z (k  );
        cp(i,k) = sSouth(i,k);
    end
end
end

figure
[CC,hh] = contour(xp,zp,cp); clabel(CC,hh); hold on;
plot(XC(:,1),h(:,1),'k+-')
title('S SECTION ALONG SOUTHERN BRY')
xlabel('LONGITUDE [\circW]'); ylabel('DEPTH [m]');


% v section

xp(1:nx,1:nz) = NaN;
zp(1:nx,1:nz) = NaN;
cp(1:nx,1:nz) = NaN;

for i=1:nx
for k=1:nz
    if (z(k) >= h(i,1)) 
        xp(i,k) = XC(i,1);
        zp(i,k) = z (k  );
        cp(i,k) = vSouth(i,k);
    end
end
end

figure
[CC,hh] = contour(xp,zp,cp); clabel(CC,hh); hold on;
plot(XC(:,1),h(:,1),'k+-')
title('v SECTION ALONG SOUTHERN BRY')
xlabel('LONGITUDE [\circW]'); ylabel('DEPTH [m]');


%% PLOT SECTIONS OF T,S,u,v ALONG EASTERN BRY

disp(' ')
disp('Plot sections along eastern bry ...')

clear yp zp cp

% T section along easthern bry

yp(1:ny,1:nz) = NaN;
zp(1:ny,1:nz) = NaN;
cp(1:ny,1:nz) = NaN;

for j=1:ny
for k=1:nz
    if (z(k) >= h(nx,j)) 
        yp(j,k) = YC   (nx,j);
        zp(j,k) = z    (   k);
        cp(j,k) = tEast(j ,k);
    end
end
end

figure
[CC,hh] = contour(yp,zp,cp); clabel(CC,hh); hold on;
plot(YC(nx,:),h(nx,:),'k+-')
title('T SECTION ALONG EASTERN BRY')
xlabel('LATITUDE [\circN]'); ylabel('DEPTH [m]');


% S section along eastern bry

yp(1:ny,1:nz) = NaN;
zp(1:ny,1:nz) = NaN;
cp(1:ny,1:nz) = NaN;

for j=1:ny
for k=1:nz
    if (z(k) >= h(nx,j)) 
        yp(j,k) = YC   (nx,j);
        zp(j,k) = z    (   k);
        cp(j,k) = sEast(j ,k);
    end
end
end

figure
[CC,hh] = contour(yp,zp,cp); clabel(CC,hh); hold on;
plot(YC(nx,:),h(nx,:),'k+-')
title('S SECTION ALONG EASTERN BRY')
xlabel('LATITUDE [\circN]'); ylabel('DEPTH [m]');


% u section along eastern bry

yp(1:ny,1:nz) = NaN;
zp(1:ny,1:nz) = NaN;
cp(1:ny,1:nz) = NaN;

for j=1:ny
for k=1:nz
    if (z(k) >= h(nx,j)) 
        yp(j,k) = YC   (nx,j);
        zp(j,k) = z    (   k);
        cp(j,k) = uEast(j ,k);
    end
end
end

figure
[CC,hh] = contour(yp,zp,cp); clabel(CC,hh); hold on;
plot(YC(nx,:),h(nx,:),'k+-')
title('u SECTION ALONG EASTERN BRY')
xlabel('LATITUDE [^o]'); ylabel('DEPTH [m]');

% v section along eastern bry

yp(1:ny,1:nz) = NaN;
zp(1:ny,1:nz) = NaN;
cp(1:ny,1:nz) = NaN;

for j=1:ny
for k=1:nz
    if (z(k) >= h(nx,j)) 
        yp(j,k) = YC   (nx,j);
        zp(j,k) = z    (   k);
        cp(j,k) = vEast(j ,k);
    end
end
end

figure
[CC,hh] = contour(yp,zp,cp); clabel(CC,hh); hold on;
plot(YC(nx,:),h(nx,:),'k+-')
title('v SECTION ALONG EASTERN BRY')
xlabel('LATITUDE [^o]'); ylabel('DEPTH [m]');


%% OPEN BOUNDARY CONDITIONS: Write in binary files

disp(' ')
disp('Write boundary conditions in binary files ...')

% Format parameters

ieee = 'b';
prec = 'real*8';

% Velocity component u

fid=fopen('uWest.bin' ,'w',ieee); fwrite(fid,uWest ,prec); fclose(fid);
fid=fopen('uEast.bin' ,'w',ieee); fwrite(fid,uEast ,prec); fclose(fid);
fid=fopen('uSouth.bin','w',ieee); fwrite(fid,uSouth,prec); fclose(fid);
fid=fopen('uNorth.bin','w',ieee); fwrite(fid,uNorth,prec); fclose(fid);

% Velocity component v

fid=fopen('vWest.bin' ,'w',ieee); fwrite(fid,vWest ,prec); fclose(fid);
fid=fopen('vEast.bin' ,'w',ieee); fwrite(fid,vEast ,prec); fclose(fid);
fid=fopen('vSouth.bin','w',ieee); fwrite(fid,vSouth,prec); fclose(fid);
fid=fopen('vNorth.bin','w',ieee); fwrite(fid,vNorth,prec); fclose(fid);

% Temperature

fid=fopen('tWest.bin' ,'w',ieee); fwrite(fid,tWest ,prec); fclose(fid);
fid=fopen('tEast.bin' ,'w',ieee); fwrite(fid,tEast ,prec); fclose(fid);
fid=fopen('tSouth.bin','w',ieee); fwrite(fid,tSouth,prec); fclose(fid);
fid=fopen('tNorth.bin','w',ieee); fwrite(fid,tNorth,prec); fclose(fid);

% Salinity

fid=fopen('sWest.bin' ,'w',ieee); fwrite(fid,sWest ,prec); fclose(fid);
fid=fopen('sEast.bin' ,'w',ieee); fwrite(fid,sEast ,prec); fclose(fid);
fid=fopen('sSouth.bin','w',ieee); fwrite(fid,sSouth,prec); fclose(fid);
fid=fopen('sNorth.bin','w',ieee); fwrite(fid,sNorth,prec); fclose(fid);

% Free surface elevation

% fid=fopen('etaEast.bin' ,'w',ieee); fwrite(fid,etaEast ,prec); fclose(fid);
% fid=fopen('etaSouth.bin','w',ieee); fwrite(fid,etaSouth,prec); fclose(fid);


%% ENDING STATEMENT

disp(' ')
disp('End of run')
disp(' ')
