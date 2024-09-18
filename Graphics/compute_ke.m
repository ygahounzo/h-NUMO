% This program compute the total kinetic energy
% Written by : Yao Gahounzo

clear all
warning('off', 'all')

file_num_start = 1;     % start animating with this file number
file_num_end = 720;       % end animating with this file number
delta_file_num = 1;      % what is the skip between the file numbers

n = ceil(file_num_end / delta_file_num);

ke = zeros(n,1);
ke1 = zeros(n,1);
ke2 = zeros(n,1);

% get grid
kdm = 2;
idm = 101; jdm = 101;

bc = 'freeslip'; % or noslip
nop = 4; % or 2
nu = 50; % or 50 , viscosity
Ne = 25;

% area
L = 2e6; % Length of the simulation domain (here [0,2e6]x[0,2e6])
dx = L/(Ne*nop+1);
dy = dx;
area_numo = dx * dy;

save_folder = './KE_nop/NUMO_KE_N4';

if ~exist(save_folder, 'dir')
    mkdir(save_folder);
end

name_root = './Double_gyre/'; % make to put the correct path to your output files

ii = 0;

%READ DATA FROM MATLAB OUTPUT FILE

name_file = [name_root, sprintf('mlswe%04d', 1)];
[~,pb,ubp,vbp,dp_df,u_df,v_df,coord,dt,nk, dt_btp] = load_data_numo(name_file);

xmin=min(coord(1,:)); xmax=max(coord(1,:));
ymin=min(coord(2,:)); ymax=max(coord(2,:));
xe = coord(1,:);
ye = coord(2,:);
nxx = 200; nyy = 200;
dx = (xmax-xmin)/nxx;
dy = (ymax-ymin)/nyy;
[xi,yi] = meshgrid(xmin:dx:xmax,ymin:dy:ymax);

u = zeros(nxx+1,nxx+1,nk);
v = zeros(nxx+1,nxx+1,nk);
dp = zeros(nxx+1,nxx+1,nk);

s = zeros(nxx+1,nxx+1,kdm);
vol = zeros(nxx+1,nxx+1,kdm);

for ifile = file_num_start:delta_file_num:file_num_end

    ii = ii + 1

    
    %READ DATA FROM MATLAB OUTPUT FILE
    name_file = [name_root, sprintf('mlswe%04d', ifile)];
    [~,pb,ubp,vbp,dp_df,u_df,v_df,coord,dt,nk, dt_btp] = load_data_numo(name_file);


    u(:,:,1) = griddata(xe,ye,u_df(:,1),xi,yi,'cubic');
    u(:,:,2) = griddata(xe,ye,u_df(:,2),xi,yi,'cubic');

    v(:,:,1) = griddata(xe,ye,v_df(:,1),xi,yi,'cubic');
    v(:,:,2) = griddata(xe,ye,v_df(:,2),xi,yi,'cubic');

    dp(:,:,1) = griddata(xe,ye,dp_df(:,1),xi,yi,'cubic');
    dp(:,:,2) = griddata(xe,ye,dp_df(:,2),xi,yi,'cubic');

    vol(:,:,1) = area_numo.*dp(:,:,1);
    vol(:,:,2) = area_numo.*dp(:,:,2);

    s(:,:,1) = 0.5*(u(:,:,1).^2 + v(:,:,1).^2).*vol(:,:,1);
    s(:,:,2) = 0.5*(u(:,:,2).^2 + v(:,:,2).^2).*vol(:,:,2);
    
    vol1 = vol(:,:,1);
    tmp_vol1 = sum(vol1(:),"omitnan");
    s1 = s(:,:,1);
    ke1(ii) = 10e3*sum(s1(:),"omitnan") / tmp_vol1;

    vol2 = vol(:,:,2);
    tmp_vol2 = sum(vol2(:),"omitnan");
    s2 = s(:,:,2);
    ke2(ii) = 10e3*sum(s2(:),"omitnan") / tmp_vol2;

    ke(ii) = ke1(ii) + ke2(ii);

end 

% Save data to file for plotting

Time = linspace(0,file_num_end/36,n);

temp = [n;Time'; ke1; ke2; ke];

file = [save_folder,sprintf('/ke_numo_20kmv%d_%s', nu,bc)];

eval( ['save ',  file, ' temp -ascii -double'] )

