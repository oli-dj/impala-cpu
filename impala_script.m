%% Script to test out IMPALA implementation on CPU

print = 1;
plots = 1;
pathtypes = [1, 1, 1]; %0 for raster, 1 for random, 2: TBD.
template_lengths = [4, 3, 2].^2;
multigrids = 3;
multistep = 2;
threshold = 2;

%Simulation grid size
sg_x = 256;
sg_y = 256;

% Training Image
TI = channels;
cat = unique(TI(:))';

% Multigrid
if multigrids > 1
sg_x = [sg_x; zeros(multigrids-1,1)];
sg_y = [sg_y; zeros(multigrids-1,1)];
    for i = 2:multigrids
        sg_x(i) = floor(sg_x(i-1)/multistep);
        sg_y(i) = floor(sg_y(i-1)/multistep);
    end
end

SG = NaN(sg_x(end),sg_y(end));
sg_x = flipud(sg_x);
sg_y = flipud(sg_y);

if plots
    figure;
end

tic;
plotcount = 1;
for i = 1:multigrids
    % Template
    tau = mps_template(template_lengths(i));
    
    %plots
    if plots
        subplot(multigrids,3,plotcount)
        imagesc(TI(1:round(sg_x(end)/sg_x(i)):end,...
            1:round(sg_y(end)/sg_y(i)):end));
        plotcount = plotcount + 1;
    end
    
    list = populate_impala_list( TI(1:round(sg_x(end)/sg_x(i)):end,...
        1:round(sg_y(end)/sg_y(i)):end), tau );
    toc;
    time_elapsed = toc;
    %fprintf('Time to populate list: %8.3f seconds.\n', time_elapsed);
    
    %Display list if small enough
    if size(list,1) < 40
        print_impala_list( list );
    end
    
    %Generate new path
    switch  pathtypes(i)
        case 0
            [path, n_u] = raster_path_scn398(SG);
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
        SG = NaN(sg_x(i+1),sg_y(i+1));
        SG(1:multistep:end,1:multistep:end) = SG_temp;
    else
        SG = SG_temp;
    end
    
    %%plots
    if plots
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
    end
end