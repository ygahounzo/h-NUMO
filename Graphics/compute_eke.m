clear all

% Compute the mean eddy kinetic energy over the last 5 years (15-20 years)
% Written by : Yao Gahounzo

file_num_start = 540;     % start animating with this file number
file_num_end = 720;       % end animating with this file number
delta_file_num = 1;      % what is the skip between the file numbers

n = ceil(file_num_end / delta_file_num);

g = 9.806; 
alpha = [9.7370e-04; 9.7350e-04];

idm = 101; 
jdm = 101;
nlayers = 2; % number of layers

nop = 4; % ploynomial order
nu = 50; % viscosity
Ne = 25; % number of elements

res = 20; % grid resolution

ii = 0;

name_root = './Double_gyre/'; % make sure to put the correct path to your output files

save_folder = sprintf('./EKE');

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

        name_file = [name_root, sprintf('mlswe%04d', ifile)];
    
        [~,pb,ub,vb,h,u,v,coord,dt,nk, dt_btp] = load_data_numo(name_file);
    
        %Interpolate solution
        xmin = min(coord(1,:)); xmax = max(coord(1,:));
        ymin = min(coord(2,:)); ymax = max(coord(2,:));
        dx = (xmax-xmin)/(idm-1);
        dy = (ymax-ymin)/(jdm-1);
        xe = coord(1,:);
        ye = coord(2,:);
        [xi,yi] = meshgrid(xmin:dx:xmax,ymin:dy:ymax);
    
        hi = griddata(xe,ye,h(:,ilayer),xi,yi,'linear');
        ui = griddata(xe,ye,u(:,ilayer),xi,yi,'linear');
        vi = griddata(xe,ye,v(:,ilayer),xi,yi,'linear');

        um = um + ui.*hi;
        vm = vm + vi.*hi;
        dpm = dpm + hi;

        ke = ke + 0.5*(ui.^2 + vi.^2).*hi;
    
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
    
    file = [save_folder,sprintf('/eke_%dkm_l%d', res, ilayer)];
    
    eval( ['save ',  file, ' temp -ascii -double'] )
end


