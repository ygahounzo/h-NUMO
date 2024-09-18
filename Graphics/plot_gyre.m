
warning('off', 'all')
file_num_start = 720;     % start animating with this file number
file_num_end = 720;       % end animating with this file number
delta_file_num = 10;      % what is the skip between the file numbers
                          % (must match the actual file numbers)

name_root = './Double_gyre/';

save_folder = './BB86';

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
    v = -1300:400:600;

    %READ DATA FROM MATLAB OUTPUT FILE
    name_fortran_data_file = [name_root, sprintf('mlswe%04d', ifile)];

    [npoin,pb,ubp,vbp,dp_df,udp_df,vdp_df,coord,dt,nk] = load_data_Higdon(name_fortran_data_file);
  
    % Interpolate the solution
    xmin=min(coord(1,:)); xmax=max(coord(1,:));
    ymin=min(coord(2,:)); ymax=max(coord(2,:));
    xe = coord(1,:);
    ye = coord(2,:);
    nxx = 400; nyy = 400;
    dx = (xmax-xmin)/nxx;
    dy = (ymax-ymin)/nyy;
    [xi,yi] = meshgrid(xmin:dx:xmax,ymin:dy:ymax);

    num_level = 10;

    subplot(1,2,1);
    
    uu = udp_df(:,ilayer);
    vv = vdp_df(:,ilayer);
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

    dp = dp_df(:,ilayer) - 1489.4;
    qi = griddata(xe,ye,dp,xi,yi,'cubic');
    xii = xi ./1e5; yii = yi ./1e5;
    [C,h] = contourf(xii,yii,qi,levels);
    clabel(C,h,v);
    clim([-1000,1000])
    colorbar()
    axis equal
    hold off

    sgtitle(sprintf('MLSWE dt = %.f s at %d days',dt,ifile*10),'FontSize',28)
    
    
end


function [npoin,pb,ubp,vbp,dp_df,udp_df,vdp_df,coord,dt,nk] = load_data_Higdon(name_fortran_data_file)

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

%     dim = [npoin,nk+1];
%     z = temp(count: (count+prod(dim)-1));  count=count+prod(dim);
%     z = reshape(z, dim);

end




