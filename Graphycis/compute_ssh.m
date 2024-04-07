% This program compute the SSH
% Written by : Yao Gahounzo

clear all
warning('off', 'all')

file_num_start = 1;     % start animating with this file number
file_num_end = 720;       % end animating with this file number
delta_file_num = 1;      % what is the skip between the file numbers


bc = 'noslip'; % or noslip
nop = 2; % 2, 4
nu = 50; % 50, 500 viscosity
Ne = 50;

name_root = sprintf('./FALCON/N%d_Ne%d_v%d_%s/outbora',nop,Ne,nu,bc);
% name_root = sprintf('/Volumes/GLACIE/NUMO_Data_viz/N%d_Ne%d_v%d_%s/outbora',nop,Ne,nu,bc);

ii = 0;

% Bottom depth
layer_dz_eq = [1489.5,8438.5];
depth = sum(layer_dz_eq);

% save_folder = sprintf('./NUMO_SSH/ssh_numo_ne%dv%d_N%d_%s/', Ne,nu, nop,bc);
save_folder = sprintf('./NUMO_SSH_v1/ssh_numo_v%d_%s_unstruct/',nu,bc);

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

    name_file = [name_root, sprintf('%04d', ifile)];

    [npoin,pb,ub,vb,h,u,v,coord,dt,nk, dt_btp,~] = load_data_numo1(name_file);

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
    %idm = 101; jdm = 101;
    dx = 20e3;
    dy = 20e3;
    [xi,yi] = meshgrid(xmin:dx:xmax,ymin:dy:ymax);

    sshi = griddata(xe,ye,ssh,xi,yi,'cubic');

    fprintf(fileID1,'%23.16f\n', xi./1e3);
    fprintf(fileID1,'%23.16f\n', yi./1e3);
    fprintf(fileID1,'%23.16f\n', sshi); 

    fclose(fileID1);

end 



