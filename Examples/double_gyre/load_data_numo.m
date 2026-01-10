function [npoin,pb,ub,vb,dp,u,v,coord,dt,nk, dt_btp,z] = load_data_numo(name_fortran_data_file)

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