
warning('off', 'all')
clear all

fig = figure('units','inches','Position',[10 10 8 5]);

file_num_start = 108;     % start animating with this file number
file_num_end = 108;       % end animating with this file number
delta_file_num = 10;      % what is the skip between the file numbers
                          % (must match the actual file numbers)

name_root = './lake/mlswe';

ii = 0;

ilayer = 1;


for ifile = file_num_start:delta_file_num:file_num_end

    ii = ii + 1;

    nstep = ifile;

    %READ DATA FROM MATLAB OUTPUT FILE
    name_fortran_data_file = [name_root, sprintf('%04d', ifile)];

    [npoin,pb,ub,vb,dp,u,v,coord,z,nk] = load_data(name_fortran_data_file);

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
    z2 = z(:,2);
    z3 = z(:,3);

    for k = 1:nk+1

        zk = z(:,k);

        qi = griddata(xe,ye,zk,xi,yi,'cubic');

        C0 = zeros([size(qi),3]);
        C0(:,:,1) = zeros(size(qi));
        C0(:,:,2) = 0.5*ones(size(qi));
        C0(:,:,3) = ones(size(qi));

        if(k == 1) 
            surf(xi,yi,qi,C0,'FaceColor', 'interp','EdgeColor','none',FaceAlpha=0.7);

            CC = C0;

        elseif(k == nk+1)

            z3 = z(:,nk+1);
            qi = griddata(xe,ye,zk,xi,yi,'linear');

            mesh(xi,yi,qi,'FaceColor', 'interp','EdgeColor','none');
            cmocean('topo')
        else

            mesh(xi,yi,qi,CC,'FaceColor', 'interp','EdgeColor','none');

        end

        hold on
    end
    zlim([-40,0])

    xlabel('X [m]');
    ylabel('Y [m]');
    hold off
    set(findall(fig,'-property','FontSize'),'FontSize',16,'FontName','Times')
    view(-15,15)

end


function [npoin,pb,ub,vb,dp,u,v,coord,z,nk] = load_data(name_fortran_data_file)

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




