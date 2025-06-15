
warning('off', 'all')
file_num_start = 720;     % start animating with this file number
file_num_end = 720;       % end animating with this file number
delta_file_num = 10;      % what is the skip between the file numbers
                          % (must match the actual file numbers)

name_root = './';     % path to the output files

save_folder = './BB86'; % path to where to save the plots

if ~exist(save_folder, 'dir')
    mkdir(save_folder);
end

layer_dz_eq = [1489.5,8438.4];

ii = 0;

ilayer = 1;
nd_view = 2;

figure('Position',[1 1 1000 500])

f = sprintf('%s/bb86',save_folder);

myVideo = VideoWriter(f,'MPEG-4'); %open video file
myVideo.FrameRate = 10; 
open(myVideo)


for ifile = file_num_start:delta_file_num:file_num_end
 
    ii = ii + 1;

    nstep = ifile;

    levels = -1200:100:600;
    vl = -1300:400:600;

    %READ DATA FROM MATLAB OUTPUT FILE
    name_fortran_data_file = [name_root, sprintf('mlswe%04d', ifile)];

    [npoin,pb,ub,vb,dp,u,v,coord,dt,nk] = load_data_numo(name_fortran_data_file);
  
    % Interpolate the solution
    xmin=min(coord(1,:)); xmax=max(coord(1,:));
    ymin=min(coord(2,:)); ymax=max(coord(2,:));
    xe = coord(1,:);
    ye = coord(2,:);
    nxx = 401; nyy = 401;
    dx = (xmax-xmin)/nxx;
    dy = (ymax-ymin)/nyy;
    [xi,yi] = meshgrid(xmin:dx:xmax,ymin:dy:ymax);

    num_level = 10;

    subplot(1,2,1);
    
    uu = u(:,ilayer);
    vv = v(:,ilayer);
    uui = griddata(xe,ye,uu,xi,yi,'cubic');
    vvi = griddata(xe,ye,vv,xi,yi,'cubic');
    [n,m] = size(uui);
    xii = xi ./1e5; yii = yi ./1e5;
    
    contourf(xii,yii,uui,num_level);hold on
    colorbar()
    normeref_hyc = 0.02;
    normeref = normeref_hyc ; %% in m
    norme =sqrt(uui.^2+vvi.^2) ;
    normemax = max(max(norme))/normeref;
    sc= 4;
    idx = 1:sc:n;
    idy = 1:sc:m;
    axis equal
    hold off

    subplot(1,2,2);

    dp = dp(:,ilayer) - 1489.4;
    qi = griddata(xe,ye,dp,xi,yi,'cubic');
    xii = xi ./1e5; yii = yi ./1e5;
    [C,h] = contourf(xii,yii,qi,levels);
    clabel(C,h,vl);
    clim([-1000,1000])
    colorbar()
    axis equal
    hold off

    sgtitle(sprintf('MLSWE dt = %.f s at %d days',dt,ifile*10),'FontSize',28)
    
    
end


