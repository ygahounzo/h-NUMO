
function vals = get_cmapeke()
    N = 30;
    vals = ones(N, 3);
    
    % Read cmap_eke
    fid = fopen('./cmap_eke.txt', 'r');
    
    % Loop over lines and extract variables of interest
    t = 1;
    while ~feof(fid)
        line = fgetl(fid);
        % Split line into columns
        columns = str2num(line);
        
        % Update values
        vals(t, 1) = columns(1) / 256;
        vals(t, 2) = columns(2) / 256;
        vals(t, 3) = columns(3) / 256;
        
        t = t + 1;
    end

    fclose(fid);
    
%     cmapeke = colormap(vals);
end
