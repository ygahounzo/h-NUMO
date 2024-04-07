% This program plot time-averaged total kinetic energy from 15-20 years
% Written by : Yao Gahounzo

clear all
warning('off', 'all')
% fig = figure('Position',[1 1 1200 1000]);  
fig = figure('units','inches','Position',[10 10 9 7]);

bc = 'noslip';
bc1 = 'no-slip';

% viscosity nu = 500 m^2/s and noslip
nu = 500;

Ne_2 = [100,67,50,40,33];
Ne_4 = [50,34,25,20,17];
Ne_6 = [32,22,16,14,11];

nop = [2,4,6];

nyears_avg = 5; % last 5 years

Nu = [50,500];

for iNu = 1:length(Nu)

    nu = Nu(iNu);

    if(nu == 50) 
        Ne_2 = [200,100,67,50,40];
        Ne_4 = [100,50,34,25,20];
        Ne_6 = [66,32,22,16,14];

        res = [5,10,15,20,25];

        cc = ['b','k','r'];
        mm = '--s';
        ms = 10;
    else
        Ne_2 = [200,100,67,50,40];
        Ne_4 = [100,50,34,25,20];
        Ne_6 = [66,32,22,16,14];

        res = [5,10,15,20,25];

        cc = ['b','k','r'];
        mm = '--.';
        ms = 40;
    end

    NE = [Ne_2;Ne_4;Ne_6];

    for inop = 1:length(nop)
    
        N = nop(inop);
        
        nne = length(res);
        KE_avg = zeros(nne,1);
    
        for ie = 1:nne
        
            ires = res(ie);
    
            if(N == 6 && nu == 500 && ires == 5)
                % continue
                total_years = 17; % total rum model years
            else
                total_years = 20; % total rum model years
            end

            % if(N==2)
            %     name_root_N = sprintf('./KE_nop2/NUMO_KE_N%d/ke_bb86_%dkmv%d_%s',N,ires,nu,bc);  % numo ke
            % else
                name_root_N = sprintf('./KE_nop/NUMO_KE_N%d/ke_bb86_%dkmv%d_%s',N,ires,nu,bc);  % numo ke
                
            % end

            
            [Time,ke_l1,ke_l2,ket] = load_data_ke(name_root_N);
    
            
            n1 = (total_years-nyears_avg)*36; % get the index of the start from years
            n2 = total_years*36-n1 + 1;  % index of the last years

            NNN = total_years*36;
            
            KE_avg(ie) = sum(ket(n1:NNN))/n2;
    
        end
        
        KE_avg
    
        Dx = 2000./(N*NE(inop,:));
        if(N == 2)
            nn = nne;
            if(nne > 5); nn = nne-1; end 
            plot(Dx(1:nn),KE_avg(1:nn),mm,'Color',cc(inop),'MarkerSize',ms, 'LineWidth',2); hold on
        % elseif(N == 6 && nu ==500)
        %     plot(Dx(2:end),KE_avg(2:end),mm,'Color',cc(inop),'MarkerSize',ms, 'LineWidth',2); hold on
        else
            plot(Dx(1:nne),KE_avg,mm,'Color',cc(inop),'MarkerSize',ms, 'LineWidth',2); hold on
        end
    
    end

    % HYCOM
    KE_avg_hycom = zeros(3,1);
    
    total_years = 20;

    n1 = (total_years-nyears_avg)*36;
    n2 = total_years*36-n1 + 1; 

    name_root_hycom = sprintf('./HYCOM_KE_v1/ke_bb86_4kmv%d_%s',nu,bc); % hycom ke     
    [Time,ke_l1,ke_l2,ket_4km] = load_data_ke(name_root_hycom);
    KE_avg_hycom(1) = sum(ket_4km(n1:end))/n2;

    name_root_hycom = sprintf('./HYCOM_KE_v1/ke_bb86_10kmv%d_%s',nu,bc); % hycom ke     
    [Time,ke_l1,ke_l2,ket_10km] = load_data_ke(name_root_hycom);
    
    KE_avg_hycom(2) = sum(ket_10km(n1:end))/n2;
    
    name_root_hycom = sprintf('./HYCOM_KE_v1/ke_bb86_20kmv%d_%s',nu,bc); % hycom ke     
    [Time,ke_l1,ke_l2,ket_20km] = load_data_ke(name_root_hycom);
    KE_avg_hycom(3) = sum(ket_20km(n1:end))/n2;
    
    Dx_hycom = [4,10,20];
    plot(Dx_hycom,KE_avg_hycom,mm,'Color',"#77AC30",'MarkerSize',ms, 'LineWidth',2); hold on

end

grid on
xlim([0,30])
xlabel('Effective resolution \Delta x')
ylabel(sprintf('KE [cm^2/s^2]: last %d years',nyears_avg))

gd = {'\nu = 50, N = 2', '\nu = 50, N = 4', '\nu = 50, N = 6',...
    '\nu = 50, hycom', ...
    '\nu = 500, N = 2', '\nu = 500, N = 4', '\nu = 500, N = 6',...
    '\nu = 500 hycom'};
lgd = legend(gd,'Location','best');
% lgd.NumColumns = 4;
% 
% lgd.Position(1) = 0.22;
% lgd.Position(2) = 0.86;

% sgtitle(sprintf("Mean Kinetic Energy: %s",bc1))
set(findall(fig,'-property','FontSize'),'FontSize',18,'FontName','Times')
% print('-vector',sprintf('./KET/ke_convg_%s.eps',bc),'-depsc')




