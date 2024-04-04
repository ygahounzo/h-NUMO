
function [dp,u,v,plon,plat] = load_data_hycom(name_data_file,ni,nj,kdm)

    %   Load the data written by the Fortran DG code.

    temp = load(name_data_file, '-ascii');

    count = 1;

    dim = [ni,nj];
    plon = temp(count: (count+prod(dim)-1));  count=count+prod(dim);
    plon = reshape(plon, dim);

    dim = [ni,nj];
    plat = temp(count: (count+prod(dim)-1));  count=count+prod(dim);
    plat = reshape(plat, dim);
    
    dim = [ni,nj,kdm];
    dp = temp(count: (count+prod(dim)-1));  count=count+prod(dim);
    dp = reshape(dp, dim);

    dim = [ni,nj,kdm];
    u = temp(count: (count+prod(dim)-1));  count=count+prod(dim);
    u = reshape(u, dim);

    dim = [ni,nj,kdm];
    v = temp(count: (count+prod(dim)-1));  count=count+prod(dim);
    v = reshape(v, dim);

end