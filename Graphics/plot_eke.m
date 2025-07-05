%This program plot eddy kinetic energy EKE
% Written by : Yao Gahounzo

clear all

nu = 50;  % viscosity
res = 20; % grid resolution
lim_eke = 400; % limit for EKE

ilayer = 1;  % layer to visualize

name_root = sprintf('./EKE/eke_%dkm_l%d',res,ilayer);

[idm,jdm,mke_numo,eke_numo] = load_data_eke(name_root);

x = linspace(0, 2000, idm);
[plon,plat] = meshgrid(x, x);

% fig = figure('Position',[1 1 1200 800]);
fig = figure('units','inches','Position',[10 10 8 8]);
tiledlayout(1,2,'TileSpacing','Compact','Padding','Compact');

vals = get_cmapeke();
levels = 0:40:lim_eke;
% Eddy Kinetic Energy

contourf(plon,plat,eke_numo,levels);
% colorbar('eastOutside')
colormap(vals)
caxis([0 lim_eke])
xticks(linspace(0,2000,5))
yticks(linspace(0,2000,5))
xlabel('X [km]')
ylabel('Y [km]')
title('h-NUMO: EKE (cm^2/s^2), 15-20 years');
hC = colorbar('eastOutside');
ylabel(hC,'EKE (cm^2/s^2)')
axis equal

set(findall(fig,'-property','FontSize'),'FontSize',18,'FontName','Times')


function [idm,jdm,kem,eke] = load_data_eke(name_file)

    %   Load the data written by the Fortran DG code.

    temp = load(name_file, '-ascii');

    count = 1;
    idm = temp(count);  count=count+1;
    jdm = temp(count);  count=count+1;

    dim = [idm,jdm];
    kem = temp(count: (count+prod(dim)-1));  count=count+prod(dim);
    kem = reshape(kem, dim);

    dim = [idm,jdm];
    eke = temp(count: (count+prod(dim)-1));  count=count+prod(dim);
    eke = reshape(eke, dim);

end
