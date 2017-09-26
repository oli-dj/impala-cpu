%% Script to test out IMPALA implementation on CPU

print = 1;
plots = 1;
pathtypes = [0, 0, 0]; %0 for raster, 1 for random, 2: TBD.
template_lengths = [12, 4, 3].^2;
multigrids = 3;
multistep = 2;
threshold = 2;

%Simulation grid size
sg_x = 64;
sg_y = 64;
sg_z = 64;

%% Training Image

%2D
%TI = channels;

%3D
%TI = read_eas_matrix('ti_fluvsim_250_250_100.dat');
TI = read_eas_matrix('ti_cb_6x6_40_40_40.dat');

dim = length(size(TI));
cat = unique(TI(:))';

%% Multigrid
if multigrids > 1
    switch dim
        case 2 %2D
            sg_x = [sg_x; zeros(multigrids-1,1)];
            sg_y = [sg_y; zeros(multigrids-1,1)];
            for i = 2:multigrids
                sg_x(i) = floor(sg_x(i-1)/multistep);
                sg_y(i) = floor(sg_y(i-1)/multistep);
            end
            
            SG = NaN(sg_x(end),sg_y(end));
            sg_x = flipud(sg_x);
            sg_y = flipud(sg_y);
        case 3 %3D
            sg_x = [sg_x; zeros(multigrids-1,1)];
            sg_y = [sg_y; zeros(multigrids-1,1)];
            sg_z = [sg_z; zeros(multigrids-1,1)];
            for i = 2:multigrids
                sg_x(i) = floor(sg_x(i-1)/multistep);
                sg_y(i) = floor(sg_y(i-1)/multistep);
                sg_z(i) = floor(sg_z(i-1)/multistep);
            end
            
            SG = NaN(sg_x(end),sg_y(end),sg_z(end));
            sg_x = flipud(sg_x);
            sg_y = flipud(sg_y);
            sg_z = flipud(sg_z);
            
    end
else
    switch dim
        case 2
            SG = NaN(sg_x,sg_y);
        case 3
            SG = NaN(sg_x,sg_y,sg_z);
    end
end

if plots
    figure;
end

%% Simulation

tic;
plotcount = 1;

for i = 1:multigrids
    % Template
    tau = mps_template(template_lengths(i),dim);
    
    %plots
    if plots
        switch dim
            case 2
                subplot(multigrids,3,plotcount)
                imagesc(TI(1:round(sg_x(end)/sg_x(i)):end,...
                    1:round(sg_y(end)/sg_y(i)):end));
                plotcount = plotcount + 1;
            case 3
                subplot(multigrids,3,plotcount)
                slice(TI(1:round(sg_x(end)/sg_x(i)):end,...
                    1:round(sg_y(end)/sg_y(i)):end,...
                    1:round(sg_z(end)/sg_z(i)):end),[10 25], [10 25], [5 12]);
                
                plotcount = plotcount + 1;
        end
    end
    switch length(size(TI))
        case 2
            list = populate_impala_list(...
                TI(1:round(sg_x(end)/sg_x(i)):end,...
                1:round(sg_y(end)/sg_y(i)):end), tau );
        case 3
            list = populate_impala_list(...
                TI(1:round(sg_x(end)/sg_x(i)):end,...
                1:round(sg_y(end)/sg_y(i)):end,...
                1:round(sg_z(end)/sg_z(i)):end), tau );
    end
    toc;
    time_elapsed = toc;
    if print
        fprintf('Time to populate list: %8.3f seconds.\n', time_elapsed);
    end
    
    %Display list if small enough
    if size(list,1) < 40
        print_impala_list( list );
    end
    
    %Generate new path
    switch  pathtypes(i)
        case 0
            [path, n_u] = raster_path(SG);
        case 1
            [path, n_u] = rand_path(SG);
    end
    
    %Pre-calculate random numbers
    rand_pre = rand(n_u,1);
    
    %Run IMPALA core
    [SG_temp, tauG] = impala_core(SG, list, path, tau, rand_pre, cat,...
        threshold, print);
    
    %Copy simulation grid points onto new, bigger grid
    %TODO: support non integer size factors.
    if i < multigrids
        switch dim
            case 2
                SG = NaN(sg_x(i+1),sg_y(i+1));
                SG(1:multistep:end,1:multistep:end) = SG_temp;
            case 3
                SG = NaN(sg_x(i+1),sg_y(i+1),sg_z(i+1));
                SG(1:multistep:end,1:multistep:end,1:multistep:end) = SG_temp;
                
        end
    else
        SG = SG_temp;
    end
    
    %%plots
    if plots
        switch dim
            case 2
                subplot(multigrids,3,plotcount)
                imagesc(SG_temp);
                plotcount = plotcount + 1;
                subplot(multigrids,3,plotcount)
                imagesc(tauG);
                colorbar;
                c_limits = [0 template_lengths(i)];
                caxis(c_limits);
                plotcount = plotcount + 1;
                drawnow;
            case 3
                subplot(multigrids,3,plotcount)
                slice(SG,[10 25], [10 25], [5 12]);
                plotcount = plotcount + 1;
                subplot(multigrids,3,plotcount)
                slice(tauG,[10 25], [10 25], [5 12]);
                plotcount = plotcount + 1;
                drawnow;
        end
    end
end