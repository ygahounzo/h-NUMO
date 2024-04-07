%This program plot eddy kinetic energy EKE
% Written by : Yao Gahounzo

clear all

kdm = 2;  % number of layers

nu = 50;  % viscosity
bc = 'noslip';  % boundary condition
bc1 = 'free-slip';

res = 10;

lim_kem = 400;  % limit for mean KE
lim_eke = 400; % limit for EKE

ilayer = 1;  % layer to visualize

name_root = sprintf('./EKE_numo/EKE_N2/eke_%dkmv%d_%s_l%d',res,nu,bc,ilayer);
name_rooth = sprintf('./EKE_HYCOM/eke_%dkmv%d_%s_l%d',res,nu,bc,ilayer);

[idm,jdm,mke_numo,eke_numo] = load_data_eke(name_root);
[idm1,jdm1,mke_hycom,eke_hycom] = load_data_eke(name_rooth);


x = linspace(0, 2000, idm);
[plon,plat] = meshgrid(x, x);

x = linspace(0, 2000, idm1);
[plonh,plath] = meshgrid(x, x);


% fig = figure('Position',[1 1 1200 800]);
fig = figure('units','inches','Position',[10 10 8 8]);
tiledlayout(1,2,'TileSpacing','Compact','Padding','Compact');

vals = get_cmapeke();

% Eddy Kinetic Energy

nexttile
contourf(plonh,plath,eke_hycom);
% hC = colorbar('eastOutside');
colormap(vals)
caxis([0 lim_eke])
xticks(linspace(0,2000,5))
yticks(linspace(0,2000,5))
xlabel('X [km]')
ylabel('Y [km]')
set(gca,'Yticklabel',[]) 
set(gca,'Xticklabel',[])
title('HYCOM: EKE (cm^2/s^2), 15-20 years');
set(gca, 'FontSize', 12);
ylabel(hC,'EKE (cm^2/s^2)')
axis equal

nexttile
contourf(plon,plat,eke_numo6);
% colorbar('eastOutside')
colormap(vals)
caxis([0 lim_kem])
xticks(linspace(0,2000,5))
yticks(linspace(0,2000,5))
xlabel('X [km]')
ylabel('Y [km]')
title('h-NUMO: EKE (cm^2/s^2), 15-20 years');
hC = colorbar('eastOutside');
ylabel(hC,'EKE (cm^2/s^2)')
axis equal

set(findall(fig,'-property','FontSize'),'FontSize',18,'FontName','Times')
% print('-vector',sprintf('./KET/eke_N4_%s.eps',bc),'-depsc')


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


