function [npoin,pb,ubp,vbp,dp_df,udp_df,vdp_df,coord,dt,nk, dt_btp] = load_data_numo(name_fortran_data_file)

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