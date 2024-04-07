% This program plot the SSH at a given time
% Written by : Yao Gahounzo

clear all
warning('off', 'all')
file_num_start = 360;     % start animating with this file number
file_num_end = 360;       % end animating with this file number
delta_file_num = 1;      % what is the skip between the file numbers


bc = 'noslip';
% bc1 = 'no-slip';

nop = 4; % 2, 4
nu = 50; % 50, 500 viscosity

Res = [10,20];
Years = [10,20];

code = 'hycom';

% save_folder = './SSH/NUMO';
save_folder = './SSH/HYCOM';

if ~exist(save_folder, 'dir')
    mkdir(save_folder);
end

for iy = 1:length(Years)
    nyrs = Years(iy)
    ifile = 36*nyrs;

    for ires = 1:length(Res)
        res = Res(ires)
        
        if(res == 20)
            idm = 101;
            jdm = 101;
            Ne = 25;
        else
            idm = 201;
            jdm = 201;
            Ne = 50;
        end

        kdm = 2;
        
        if strcmp(code,'numo')
            name_root = sprintf('./NUMO_SSH/ssh_numo_ne%dv%d_N%d_%s/numo',Ne,nu,nop,bc);  % numo 
        elseif strcmp(code,'hycom')

            name_root = sprintf('./HYCOM_SSH/hycom_ssh_%dkmv%d_%s/hycom',res,nu,bc); % hycom
        end
        
        
        num_level = 20;
        levels = -0.26:0.02:0.26;
        levels = 1e2*levels;
        v = -0.26:0.04:0.26;
        v = 1e2*v; 
        
        name_numo = [name_root, sprintf('%04d', ifile)];
        [X_numo,Y_numo,ssh_numo] = load_data_ssh(name_numo,idm,jdm);
        
        ssh_numo = 1e2*ssh_numo;
        
        fig = figure('units','inches','Position',[10 10 4 4]);
        % fig = figure('Position',[1 1 1200 1200]);

        [C,h] = contourf(X_numo,Y_numo,ssh_numo,levels);
        % h.LineWidth = 0.1;
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
        % title(sprintf('%s: \\nu = %d m^2/s, res = %d, y = %d',code,nu,res,nyrs))
        
        % hold off
      
        set(findall(fig,'-property','FontSize'),'FontSize',18,'fontweight', 'bold','FontName','Times')
    
        if strcmp(code,'numo')
            print('-vector',sprintf('./%s/ssh_%dkmv%d_N%d_%dy_%s.eps',save_folder,res,nu,nop,nyrs,bc),'-depsc')
        elseif strcmp(code,'hycom')

            % saveas(gcf,sprintf('./%s/ssh_%dkmv%d_%dy_%s.eps',save_folder,res,nu,nyrs,bc),'epsc')
            print('-vector',sprintf('./%s/ssh_%dkmv%d_%dy_%s.eps',save_folder,res,nu,nyrs,bc),'-depsc')
        end
        clear gcf
    end
end


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

