clear all

% Compute the mean eddy kinetic energy over the last 5 years (15-20 years)
% Written by : Yao Gahounzo

file_num_start = 540;     % start animating with this file number
file_num_end = 720;       % end animating with this file number
delta_file_num = 1;      % what is the skip between the file numbers

n = ceil(file_num_end / delta_file_num);

g = 9.806; 
alpha = [9.7370e-04; 9.7350e-04];

idm = 201; 
jdm = 201;
kdm = 2;

bc = 'noslip'; % or noslip
nop = 4; % or 2
nu = 50; % or 50 , viscosity
Ne = 25;

res = 10;

code = 'numo';

% constants
rho = 1000.;    % reference density

ii = 0;
nlayers = 2;

Nop = [2,4,6];
% NEE = [50,25,16];
NEE = [100,50,32];

BC = ['noslip','freeslip'];

for inop = 1:length(Nop)

    nop = Nop(inop);
    Ne = NEE(inop);

    for ibc = 1:2

        if(ibc == 1); bc = 'noslip'; end
        if(ibc == 2); bc = 'freeslip'; end

        name_root = sprintf('/Volumes/GLACIE/NUMO/N%d_Ne%d_v%d_%s/outbora',nop,Ne,nu,bc);
        % name_root = sprintf('../BB86_PACKAGE/MATLAB/hycom_%dkmv%d_20y_%s/hycom',res,nu,bc);
    
        save_folder = sprintf('./EKE_numo/EKE_N%d',nop);
        % save_folder = './EKE_HYCOM';
        
        if ~exist(save_folder, 'dir')
            mkdir(save_folder);
        end
    
        for ilayer = 1:nlayers
        
            um = zeros(idm,jdm);
            vm = zeros(idm,jdm);
            dpm = zeros(idm,jdm);
            ke = zeros(idm,jdm);
        
            for ifile = file_num_start:delta_file_num:file_num_end
            
                ii = ii + 1
            
                % Compute KE for NUMO
            
                if strcmp(code,'numo')
            
                    name_file = [name_root, sprintf('%04d', ifile)];
                
                    [~,pb,ubp,vbp,dp_df,udp_df,vdp_df,coord,dt,nk, dt_btp] = load_data_numo(name_file);
                
                    %Interpolate solution
                    xmin = min(coord(1,:)); xmax = max(coord(1,:));
                    ymin = min(coord(2,:)); ymax = max(coord(2,:));
                    dx = (xmax-xmin)/(idm-1);
                    dy = (ymax-ymin)/(jdm-1);
                    xe = coord(1,:);
                    ye = coord(2,:);
                    [xi,yi] = meshgrid(xmin:dx:xmax,ymin:dy:ymax);
                
                    dpi = griddata(xe,ye,dp_df(:,ilayer),xi,yi,'linear');
                    ui = griddata(xe,ye,udp_df(:,ilayer),xi,yi,'linear');
                    vi = griddata(xe,ye,vdp_df(:,ilayer),xi,yi,'linear');
            
                    um = um + ui.*dpi;
                    vm = vm + vi.*dpi;
                    dpm = dpm + dpi;
            
                    ke = ke + 0.5*(ui.^2 + vi.^2).*dpi;
            
                else % HYCOM
            
                    name_data_file2 = [name_root, sprintf('%04d', ifile)];
                    [dph,uh,vh, plon,plat] = load_data_hycom(name_data_file2, idm, jdm, kdm);
            
                    um = um + uh(:,:,ilayer).*dph(:,:,ilayer);
                    vm = vm + vh(:,:,ilayer).*dph(:,:,ilayer);
                    dpm = dpm + dph(:,:,ilayer);
            
                    ke = ke + 0.5*(uh(:,:,ilayer).^2 + vh(:,:,ilayer).^2).*dph(:,:,ilayer);
            
                end
            
            end 
        
            % ke of mean flow
            um = um ./dpm;
            vm = vm ./dpm;
            kem = 0.5*(um.^2 + vm.^2);
            
            % mean of ke
            mke = ke ./ dpm;
            
            % eke = mean of ke - ke of mean flow
            eke = mke - kem;
        
            mke = 1e4*mke; % convert to cm^2/s^2
            eke = 1e4*eke; %convert to cm^2/s^2

            % save data to file
            
            N = idm*jdm;
            temp = [idm;jdm;reshape(mke,N,1);reshape(eke,N,1)];
            
            file = [save_folder,sprintf('/eke_%dkmv%d_%s_l%d_unstruct', res, nu,bc,ilayer)];
            
            eval( ['save ',  file, ' temp -ascii -double'] )
        end 
    end
end


