
clear all
warning('off', 'all')
file_num_start = 0;     % start animating with this file number
file_num_end = 108;       % end animating with this file number
delta_file_num = 1;      % what is the skip between the file numbers
                          % (must match the actual file numbers)

ilayer = 1;
n = ceil(file_num_end/delta_file_num) + 1;

fig = figure('units','inches','Position',[10 10 8 5]);

name_root = './lake/mlswe';
ifile = 108;

name_fortran_data_file = [name_root, sprintf('%04d', ifile)];
[npoin,pb,ub,vb,dp,u,v,coord,z,dt,nk] = load_data(name_fortran_data_file);

xmin=min(coord(1,:)); xmax=max(coord(1,:));
ymin=min(coord(2,:)); ymax=max(coord(2,:));
nxx = 101; nyy = 101;
dx = (xmax-xmin)/nxx;
dy = (ymax-ymin)/nyy;
[xi,yi] = meshgrid(xmin:dx:xmax,ymin:dy:ymax);

x = xmin:dx:xmax;
xe = coord(1,:); 
ye = coord(2,:);

J = 1:2:length(x);

qi = griddata(xe,ye,u(:,1),xi,yi,'natural');
z1 = qi(:,ceil(nyy/2));

plot(x(J),z1(J),'-.', "LineWidth",3, 'MarkerSize',8,'Marker', 's'); hold on

hold on 

legend('2 layers', 'Location','northwest')
grid on
set(findall(fig,'-property','FontSize'),'FontSize',16,'FontName','Times')%,'fontweight', 'bold')

ax = gca;
ax.XAxis.FontSize = 20;
ax.YAxis.FontSize = 20;

function [npoin,pb,ub,vb,dp,u,v,coord,z,dt,nk] = load_data(name_fortran_data_file)

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
