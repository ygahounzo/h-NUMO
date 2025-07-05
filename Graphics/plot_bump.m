
warning('off', 'all')
file_num_start = 0;     % start animating with this file number
file_num_end = 100;       % end animating with this file number
delta_file_num = 1;      % what is the skip between the file numbers
                          % (must match the actual file numbers)

name_root = './bump/mlswe';

save_folder = './Ani_bump2';

if ~exist(save_folder, 'dir')
    mkdir(save_folder);
end

layer_dz_eq = [20, 20];

ii = 0;

ilayer = 1;
nd_view = 2;

name_fortran_data_file = [name_root, sprintf('%04d', 0)];

[npoin,pb,ubp,vbp,dp_df,udp_df,vdp_df,coord,z] = load_data(name_fortran_data_file);

fig = figure('Position',[1 1 1600 1200]);
% fig = figure('units','inches','Position',[10 10 8 4]);

f = sprintf('%s/bump_layer',save_folder);

myVideo = VideoWriter(f,'MPEG-4'); %open video file
myVideo.FrameRate = 10; 
open(myVideo)

for ifile = file_num_start:delta_file_num:file_num_end

    ii = ii + 1;

    nstep = ifile;

    %READ DATA FROM MATLAB OUTPUT FILE
    name_fortran_data_file = [name_root, sprintf('%04d', ifile)];
    [npoin,pb,ubp,vbp,dp_df,udp_df,vdp_df,coord,z,dt,nk] = load_data(name_fortran_data_file);

    %Plot Exact Solution
    xmin=min(coord(1,:)); xmax=max(coord(1,:));
    ymin=min(coord(2,:)); ymax=max(coord(2,:));
    xe = coord(1,:);
    ye = coord(2,:);
    nxx = 200; nyy = 200;
    dx = (xmax-xmin)/nxx;
    dy = (ymax-ymin)/nyy;
    [xi,yi] = meshgrid(xmin:dx:xmax,ymin:dy:ymax);

    z1 = z(:,1);
    z2 = z(:,2) + 20; z2 = z2*20-20;
    z3 = z(:,3);

    qi = griddata(xe,ye,z1,xi,yi,'cubic');

    C0 = zeros([size(qi),3]);
    C0(:,:,1) = zeros(size(qi));
    C0(:,:,2) = 0.5*ones(size(qi))+1e4*qi;
    C0(:,:,3) = ones(size(qi));

    xii = xi./1e3; yii = yi./1e3;

    surf(xii,yii,qi,C0,'FaceColor', 'interp','EdgeColor','none',FaceAlpha=0.7);
    hold on

    qi = griddata(xe,ye,z2,xi,yi,'cubic');
    mesh(xii,yii,qi,'FaceColor', 'interp','EdgeColor','none'); hold on
    clim([-20 -15])

    % hold on
    qi = griddata(xe,ye,z3,xi,yi,'cubic');
    mesh(xii,yii,qi,'FaceColor', 'interp','EdgeColor','none');

    xlabel('X [m]', 'FontSize',20);
    ylabel('Y [m]', 'FontSize',20);
    hold off
    zlim([-40,0])

    ax = gca;
    ax.XAxis.FontSize = 20;
    ax.YAxis.FontSize = 20;
    ax.ZAxis.FontSize = 20;
    ax.FontWeight = 'bold';

    sgtitle(sprintf('Layer interfaces at %.f min',ifile*dt/60),'FontSize',20)
    % set(findall(fig,'-property','FontSize'),'FontSize',14,'FontName','Times')
    set(findall(fig,'-property','FontSize'),'FontSize',14, 'fontweight', 'bold', 'LineWidth', 1)

    view(-15,15)
    drawnow % limitrate

    pause(0.001) %Pause and grab frame
    frame = getframe(gcf); %get frame
    writeVideo(myVideo, frame);

end

close(myVideo)

function [npoin,pb,ubp,vbp,dp_df,udp_df,vdp_df,coord,z,dt,nk] = load_data(name_fortran_data_file)

    %   Load the data written by the Fortran DG code.

    temp = load(name_fortran_data_file, '-ascii');

    count = 1;
    nk = temp(count);  count=count+1;
    npoin = temp(count);  count=count+1;
    dt = temp(count);  count=count+1;
    dt_btp = temp(count);  count=count+1;
    dim = [2,npoin];
    coord = temp(count: (count+prod(dim)-1));  count=count+prod(dim);
    coord = reshape(coord, dim);

    dim = npoin;
    pb = temp(count: (count+prod(dim)-1));  count=count+prod(dim);
    %pb = reshape(pb, dim);

    dim = npoin;
    ubp = temp(count: (count+prod(dim)-1));  count=count+prod(dim);
    %ubp = reshape(ubp, dim);

    dim = npoin;
    vbp = temp(count: (count+prod(dim)-1));  count=count+prod(dim);
    %vbp = reshape(vbp, dim);

    dim = [npoin,nk];
    dp_df = temp(count: (count+prod(dim)-1));  count=count+prod(dim);
    dp_df = reshape(dp_df, dim);

    dim = [npoin,nk];
    udp_df = temp(count: (count+prod(dim)-1));  count=count+prod(dim);
    udp_df = reshape(udp_df, dim);

    dim = [npoin,nk];
    vdp_df = temp(count: (count+prod(dim)-1));  count=count+prod(dim);
    vdp_df = reshape(vdp_df, dim);

    dim = [npoin,nk+1];
    z = temp(count: (count+prod(dim)-1));  count=count+prod(dim);
    z = reshape(z, dim);

end
