% This program plot the kinetic energy over time
% Written by : Yao Gahounzo

clear all
warning('off', 'all')

kemax = 45; % limit for KE

nu = 50; % viscosity
nop = 4; % polynomial order
        
% 20 km
res = 20; % grid resolution
name_root = sprintf('./KE/ke_numo_20km');  % numo ke

[Time,ke_numo_l1,ke_numo_l2,ket_numo] = load_data_ke(name_root);

fig = figure('Position',[1 1 1000 500]);

plot(Time,ket_numo,'LineWidth',2, 'Color','#0072BD', 'LineStyle','--'); hold on

ylim([0,kemax])
grid("on");
xlabel(sprintf("Time [years]"))
ylabel(sprintf("KE [cm^2/s^2]"))

set(findall(fig,'-property','FontSize'),'FontSize',18,'FontName','Times')


function [time,ke1,ke2,ket] = load_data_ke(name_file)

    %   Load the data written by the Fortran DG code.

    temp = load(name_file, '-ascii');

    count = 1;
    n = temp(count);  count=count+1;
    dim = n;
    time = temp(count: (count+prod(dim)-1));  count=count+prod(dim);

    dim = n;
    ke1 = temp(count: (count+prod(dim)-1));  count=count+prod(dim);

    dim = n;
    ke2 = temp(count: (count+prod(dim)-1));  count=count+prod(dim);

    dim = n;
    ket = temp(count: (count+prod(dim)-1));  count=count+prod(dim);

end
