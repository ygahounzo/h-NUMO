function [time,ke1,ke2,ket] = load_data_ke(name_file)

    %   Load the data written by the Fortran DG code.

    temp = load(name_file, '-ascii');

    count = 1;
    n = temp(count);  count=count+1;
    dim = n;
    time = temp(count: (count+prod(dim)-1));  count=count+dim;

    dim = n;
    ke1 = temp(count: (count+prod(dim)-1));  count=count+dim;

    dim = n;
    ke2 = temp(count: (count+prod(dim)-1));  count=count+dim;

    dim = n;
    ket = temp(count: (count+prod(dim)-1));  count=count+dim;


end