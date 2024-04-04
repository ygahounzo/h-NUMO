% This program plot the kinetic energy over time
% Written by : Yao Gahounzo

clear all
warning('off', 'all')

bc = 'noslip';

kemax = 25; % limit for KE

% viscosity nu = 50 or 500 m^2/s
nu = 50;

Nu = [50,500];
Nop = [2,4,6];


for inu = 1:length(Nu)
    nu = Nu(inu)

    if(strcmp(bc,'noslip') && nu == 50)
        kemax = 25;
    elseif(strcmp(bc,'noslip') && nu == 500) 
        kemax = 15;
    elseif(strcmp(bc,'freeslip') && nu == 50)
        kemax = 45;
    elseif(strcmp(bc,'freeslip') && nu == 500)
        kemax = 25;
    end

    for inop = 1:length(Nop)
        nop = Nop(inop)

        res = 10;
    
        name_root = sprintf('./KE_nop/NUMO_KE_N%d/ke_bb86_%dkmv%d_%s',nop,res,nu,bc);  % numo ke
        name_root2 = sprintf('./HYCOM_KE_v1/ke_bb86_%dkmv%d_%s',res,nu,bc); % hycom ke
        
        [Time,ke_numo_l1,ke_numo_l2,ket_numo] = load_data_ke(name_root);
        [Timeh,ke_hycom_l1,ke_hycom_l2,ket_hycom] = load_data_ke(name_root2);

        fig = figure('units','points','Position',[0 0 600 300]);
        
        plot(Time,ket_numo,'LineWidth',2, 'Color','#0072BD', 'LineStyle','-'); hold on
        plot(Timeh,ket_hycom,'LineWidth',2, 'Color',"#D95319", 'LineStyle','-');
        
        % 20 km
        res = 20;
        name_root = sprintf('./KE_nop/NUMO_KE_N%d/ke_bb86_%dkmv%d_%s',nop,res,nu,bc);  % numo ke
        name_root2 = sprintf('./HYCOM_KE_v1/ke_bb86_%dkmv%d_%s',res,nu,bc); % hycom ke
        
        [Time,ke_numo_l1,ke_numo_l2,ket_numo] = load_data_ke(name_root);
        [Timeh,ke_hycom_l1,ke_hycom_l2,ket_hycom] = load_data_ke(name_root2);
        
        plot(Time,ket_numo,'LineWidth',2, 'Color','#0072BD', 'LineStyle','--'); hold on
        plot(Timeh,ket_hycom,'LineWidth',2, 'Color',"#D95319", 'LineStyle','--'); hold off
        
        ylim([0,kemax])
        grid("on");
        xlabel(sprintf("Time [years]"))
        ylabel(sprintf("KE [cm^2/s^2]"))
        % title(sprintf("\\nu = %d m^2/s, h-NUMO: 2nd order polynomial",nu),'Fontweight','bold')
        % set(gca,'Xticklabel',[]) ;

        if(nop == 2 && nu == 50)

            lgd = legend('h-numo: 10km','hycom: 10km', 'h-numo: 20km','hycom: 20km', ...
                'Location','northoutside', 'Orientation','horizontal');
            
            hold off
            
            lgd.NumColumns = 2;
            
            lgd.Position(1) = 0.22;
            lgd.Position(2) = 0.79;

        end

        set(findall(fig,'-property','FontSize'),'FontSize',18,'FontName','Times')
        % print('-vector',sprintf('./KET/ke_10_20km_v%d_N%d_%s.eps',nu,nop,bc),'-depsc')

    end

end

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
