% This program plot animation of SSH over a given period
% Written by : Yao Gahounzo

clear all
warning('off', 'all')
file_num_start = 1;     % start animating with this file number
file_num_end = 720;       % end animating with this file number
delta_file_num = 1;      % what is the skip between the file numbers

fig = figure('Position',[10 10 800 400]); 

bc = 'freeslip';
bc1 = 'free-slip';

nop = 4; % polynomialorder
nu = 50; % viscosity
Ne = 50; % number of element in the DG 

idm = 101; % number of grid points in x
jdm = 101; % number of grid points in y
kdm = 2;  % number of layers

name_root = sprintf('./NUMO_SSH/ssh_numo_ne%dv50_N%d_%s/numo',Ne,nop,bc);  % numo  

save_folder = './SSH_movie';

if ~exist(save_folder, 'dir')
    mkdir(save_folder);
end

f = sprintf('%s/ssh_20km_numo_%s_N%d',save_folder,bc,nop);

myVideo = VideoWriter(f,'MPEG-4'); %open video file
myVideo.FrameRate = 10; 
open(myVideo)

ii = 0;

for ifile = file_num_start:delta_file_num:file_num_end
 
    ii = ii + 1

    nstep = ifile;

    num_level = 20;
    levels = -0.26:0.02:0.26;
    levels = 1e2*levels;
    v = -0.26:0.04:0.26;
    v = 1e2*v;

    % nu = 50

    name_numo = [name_root, sprintf('%04d', ifile)];
    [X_numo,Y_numo,ssh_numo] = load_data_ssh(name_numo,idm,jdm);

    name_hycom = [name_root_hycom, sprintf('%04d', ifile)];
    [X_hycom,Y_hycom,ssh_hycom] = load_data_ssh(name_hycom,idm,jdm);
    
    ssh_numo = 1e2*ssh_numo;
    ssh_hycom = 1e2*ssh_hycom;

    tcl = tiledlayout(1,1,"TileSpacing","compact");

    [C,h] = contourf(X_numo,Y_numo,ssh_numo,levels);
    clabel(C,h,v);
    % colorbar()
    clim([-25 25])

    % xlabel(sprintf("X [km]"))
    % ylabel(sprintf("Y [km]"))
    % set(gca,'XTick',[])
    set(gca,'Xticklabel',[]) ;
    set(gca,'Yticklabel',[]) ;
    title("h-NUMO")
    axis equal
    hold off

    sgtitle(sprintf('SSH: \\nu = 50 m^2/s, day = %d',ifile*10),'FontSize',15, 'fontweight', 'bold')
    set(findall(fig,'-property','FontSize'),'FontSize',15, 'fontweight', 'bold', 'LineWidth', 1)
    
    drawnow limitrate

    pause(0.001) %Pause and grab frame
    frame = getframe(gcf); %get frame
    writeVideo(myVideo, frame);

end 

close(myVideo)


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