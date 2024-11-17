% This program plot the kinetic energy over time
% Written by : Yao Gahounzo

clear all
warning('off', 'all')

bc = 'freeslip';


kemax = 90; % limit for KE

% viscosity nu = 50 or 500 m^2/s
nu = 50;
nop = 4;
        
% 20 km
res = 20;
name_root = sprintf('./KE_nop/NUMO_KE_N%d/ke_numo_%dkmv%d_%s_lfr',nop,res,nu,bc);  % numo ke

[Time,ke_numo_l1,ke_numo_l2,ket_numo] = load_data_ke(name_root);

fig = figure('Position',[1 1 1000 500]);

plot(Time,ket_numo,'LineWidth',2, 'Color','#0072BD', 'LineStyle','--'); hold on

ylim([0,kemax])
grid("on");
xlabel(sprintf("Time [years]"))
ylabel(sprintf("KE [cm^2/s^2]"))


lgd = legend('h-numo: 20km', 'Location','northoutside', 'Orientation','horizontal');

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
