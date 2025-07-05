% This program plot the SSH at a given time
% Written by : Yao Gahounzo

clear all
warning('off', 'all')

nop = 4; % polynomial order
nu = 50; % viscosity
Ne = 25; % number of elements

idm = 101;
jdm = 101;
nlayers = 2; % total number of layers

res = 20; % grid resolution
year = 20; % snapshot at that year

save_folder = './SSH';

if ~exist(save_folder, 'dir')
    mkdir(save_folder);
end

ifile = 36*year;
   
name_root = sprintf('./NUMO_SSH/ssh_numo/numo');  % numo       

num_level = 20;
levels = -0.26:0.02:0.26;
levels = 1e2*levels;
v = -0.26:0.04:0.26;
v = 1e2*v; 

name_numo = [name_root, sprintf('%04d', ifile)];
[X_numo,Y_numo,ssh_numo] = load_data_ssh(name_numo,idm,jdm);

ssh_numo = 1e2*ssh_numo;

fig = figure('Position',[1 1 1200 1200]);

[C,h] = contourf(X_numo,Y_numo,ssh_numo,levels);
clabel(C,h,v);
clim([-25 25])
xticks(linspace(0,2000,11))
yticks(linspace(0,2000,11))
axis equal
% xlabel(sprintf("X [km]"))
% ylabel(sprintf("Y [km]"))
% set(gca,'XTick','')
set(gca,'Yticklabel',[]) 
set(gca,'Xticklabel',[])

set(findall(fig,'-property','FontSize'),'FontSize',18,'fontweight', 'bold','FontName','Times')

clear gcf


function [plon,plat,ssh] = load_data_ssh(name_file,idm,jdm)

    %   Load the data written by the Fortran DG code.

    temp = load(name_file, '-ascii');

    count = 1;

    dim = [idm,jdm];
    plon = temp(count: (count+prod(dim)-1));  count=count+prod(dim);
    plon = reshape(plon, dim);

    dim = [idm,jdm];
    plat = temp(count: (count+prod(dim)-1));  count=count+prod(dim);
    plat = reshape(plat, dim);
    
    dim = [idm,jdm];
    ssh = temp(count: (count+prod(dim)-1));  count=count+prod(dim);
    ssh = reshape(ssh, dim);


end

