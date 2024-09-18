% This program compute the SSH
% Written by : Yao Gahounzo

clear all
warning('off', 'all')

file_num_start = 1;     % start animating with this file number
file_num_end = 720;       % end animating with this file number
delta_file_num = 1;      % what is the skip between the file numbers


bc = 'freeslip'; % or noslip
nop = 4; % polynomial order
nu = 50; % viscosity
Ne = 25; % number of elements 

name_root = './Double_gyre/';

ii = 0;

% Bottom depth
layer_dz_eq = [1489.5,8438.5];
depth = sum(layer_dz_eq);
save_folder = sprintf('./NUMO_SSH/ssh_numo_v%d_%s/',nu,bc);

if ~exist(save_folder, 'dir')
    mkdir(save_folder);
end

for ifile = file_num_start:delta_file_num:file_num_end

    ii = ii + 1
    
    iloop = 3 - floor(log10(ii));
    num = repmat('0', 1, iloop);
    num = [num, num2str(ii)];

    file2 = strcat(save_folder,'numo', num);
    fileID1 = fopen(file2,'W');

    % Compute KE for NUMO

    name_file = [name_root, sprintf('mlswe%04d', ifile)];

    [npoin,pb,ub,vb,h,u,v,coord,dt,nk, dt_btp] = load_data_numo(name_file);

    z = zeros(npoin,nk+1);

    z(:,nk+1) = -depth;

    for k = nk:-1:1
        z(:,k) = z(:,k+1) + h(:,k);
    end

    ssh = z(:,1);

    % Interpolate and save to a file

    xmin=min(coord(1,:)); xmax=max(coord(1,:));
    ymin=min(coord(2,:)); ymax=max(coord(2,:));
    xe = coord(1,:);
    ye = coord(2,:);
    dx = 20e3; % put the appropriate grid resolution
    dy = 20e3; % put the appropriate grid resolution
    [xi,yi] = meshgrid(xmin:dx:xmax,ymin:dy:ymax);

    sshi = griddata(xe,ye,ssh,xi,yi,'cubic');

    fprintf(fileID1,'%23.16f\n', xi./1e3);
    fprintf(fileID1,'%23.16f\n', yi./1e3);
    fprintf(fileID1,'%23.16f\n', sshi); 

    fclose(fileID1);

end 



