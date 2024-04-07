
function [ssh,plon,plat] = load_data_hycom_ssh(name_data_file,ni,nj)

    %   Load the data written by the Fortran DG code.

    temp = load(name_data_file, '-ascii');

    count = 1;

    dim = [ni,nj];
    plon = temp(count: (count+prod(dim)-1));  count=count+prod(dim);
    plon = reshape(plon, dim);

    dim = [ni,nj];
    plat = temp(count: (count+prod(dim)-1));  count=count+prod(dim);
    plat = reshape(plat, dim);
    
    dim = [ni,nj];
    ssh = temp(count: (count+prod(dim)-1));  count=count+prod(dim);
    ssh = reshape(ssh, dim);


end