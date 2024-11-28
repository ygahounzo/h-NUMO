
warning('off', 'all')
clear all

fig = figure('Position',[1 1 1600 1200]);

file_num_start = 0;     % start animating with this file number
file_num_end = 108;       % end animating with this file number
delta_file_num = 1;      % what is the skip between the file numbers
                          % (must match the actual file numbers)

name_root = './';  % path to the output files

ii = 0;

ilayer = 1;

save_folder = './bump'; % directory to save the animation

if ~exist(save_folder, 'dir')
    mkdir(save_folder);
end

f = sprintf('%s/bump',save_folder);

myVideo = VideoWriter(f,'MPEG-4'); %open video file
myVideo.FrameRate = 10; 
open(myVideo)


for ifile = file_num_start:delta_file_num:file_num_end

    ii = ii + 1;

    nstep = ifile;

    %READ DATA FROM MATLAB OUTPUT FILE
    name_data_file = [name_root, sprintf('mlswe%04d', ifile)];

    [npoin,pb,ub,vb,dp,u,v,coord,z,nk,dt] = load_data(name_data_file);

    xmin=min(coord(1,:)); xmax=max(coord(1,:));
    ymin=min(coord(2,:)); ymax=max(coord(2,:));
    xe = coord(1,:); 
    ye = coord(2,:);
    nxx = 201; nyy = 201;
    dx = (xmax-xmin)/nxx;
    dy = (ymax-ymin)/nyy;

    x = xmin:dx:xmax;
    [xi,yi] = meshgrid(xmin:dx:xmax,ymin:dy:ymax);

    z1 = z(:,1);
    z2 = z(:,2) + 20; z2 = z2*20-20;
    z3 = z(:,3);


    qi = griddata(xe,ye,z1,xi,yi,'cubic');

    C0 = zeros([size(qi),3]);
    C0(:,:,1) = zeros(size(qi));
    C0(:,:,2) = 0.5*ones(size(qi)) + 1e5*qi; 
    C0(:,:,3) = ones(size(qi));

    surf(xi,yi,qi,C0,'FaceColor', 'interp','EdgeColor','none',FaceAlpha=0.7); hold on
    
    qi = griddata(xe,ye,z2,xi,yi,'cubic');
    mesh(xi,yi,qi,'FaceColor', 'interp','EdgeColor','none'); hold on
    clim([-20 -15])

    qi = griddata(xe,ye,z3,xi,yi,'linear');
    mesh(xi,yi,qi,'FaceColor', 'interp','EdgeColor','none');
    zlim([-40,0])
    xlabel('X [m]');
    ylabel('Y [m]');

    hold off
    
    sgtitle(sprintf('Layer interfaces at %.f min',ifile*dt/60))
    set(findall(fig,'-property','FontSize'),'FontSize',16,'FontName','Times')
    view(-15,15)
    drawnow

    pause(0.001)
    frame = getframe(gcf);
    writeVideo(myVideo,frame);

end

close(myVideo)

function [npoin,pb,ub,vb,dp,u,v,coord,z,nk,dt] = load_data(name_data_file)

    %   Load the data written by the Fortran DG code.

    temp = load(name_data_file, '-ascii');

    count = 1;
    nk = temp(count);  count=count+1;
    npoin = temp(count);  count=count+1;
    dt = temp(count);  count=count+1;
    dt_btp = temp(count);  count= count+1;
    dim = [2,npoin];
    coord = temp(count: (count+prod(dim)-1));  count=count+prod(dim);
    coord = reshape(coord, dim);

    dim = npoin;
    pb = temp(count: (count+prod(dim)-1));  count=count+prod(dim);

    dim = npoin;
    upb = temp(count: (count+prod(dim)-1));  count=count+prod(dim);
    ub = upb./pb;

    dim = npoin;
    vpb = temp(count: (count+prod(dim)-1));  count=count+prod(dim);
    vb = vpb./pb;

    dim = [npoin,nk];
    dp = temp(count: (count+prod(dim)-1));  count=count+prod(dim);
    dp = reshape(dp, dim);

    dim = [npoin,nk];
    u = temp(count: (count+prod(dim)-1));  count=count+prod(dim);
    u = reshape(u, dim);

    dim = [npoin,nk];
    v = temp(count: (count+prod(dim)-1));  count=count+prod(dim);
    v = reshape(v, dim);

    dim = [npoin,nk+1];
    z = temp(count: (count+prod(dim)-1));  count=count+prod(dim);
    z = reshape(z, dim);

end




