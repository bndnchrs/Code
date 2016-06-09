%% FD_timestep_plot
% Plotting routing for simple-thermo-advect run

npcolor = 50;
skip_contour = ceil(OPTS.nt/npcolor); % Normally, we don't sub-sample, but we will for making figures

if FSTD.i == 1
    
    addpath([OPTS.path_of_code 'Core/Plotting/']);
    
    close all
    figure;
    
    cmap = [255,247,236
        254,232,200
        253,212,158
        253,187,132
        252,141,89
        239,101,72
        215,48,31
        179,0,0
        127,0,0]/256;
    
    colormap(cmap);
    
    OPTS.ax_FSD = subplot(321);
    OPTS.ax_ITD = subplot(322);
    OPTS.ax_vol = subplot(3,4,9);
    OPTS.ax_mfs = subplot(3,4,12);
    OPTS.ax_conc = subplot(3,4,10);
    OPTS.ax_thick = subplot(3,4,11);
    
    
end

if FSTD.i == OPTS.nt
    
    
    
    nplots = 5;
    skip = ceil(OPTS.nt/nplots);
    % Plotting times
    inds = 1:skip:size(DIAG.FSTD.psi,3);
    
    
    cplots = [103,0,31
        178,24,43
        214,96,77
        244,165,130
        253,219,199
        255,255,255
        224,224,224
        186,186,186
        135,135,135
        77,77,77
        26,26,26]/256;
    
    div = size(cplots,1)/nplots;
    
    for i = 1:3
        colplots(:,i) = downsample(cplots(:,i),floor(div));
    end
    
    colplots(colplots > 1) = 1;
    colplots(colplots < 0) = 0;
    
end

%%

if mod(FSTD.i,skip_contour) == 1
    
    
    contour_plot_FSDITD;
    
    line_plot_FSDITD;
    
    
    make_single_variable_plots;
    
end

%%
if DIAG.DOPLOT
    
    drawnow
    
end

if FSTD.i == OPTS.nt
    
    drawnow
    
    %   legend(ax_vol,'Predicted - Advect','Predicted - Thermo','Predicted - Both','Calculated');
    
    legend(OPTS.ax_mfs,[OPTS.h5 OPTS.h6],'By Number','By Area');
    
    pos = [16 16];
    set(gcf,'windowstyle','normal','position',[0 0 pos],'paperposition',[0 0 pos],'papersize',pos,'units','inches');
    set(gcf,'windowstyle','normal','position',[0 0 pos],'paperposition',[0 0 pos],'papersize',pos,'units','inches');
    
    if isfield(OPTS,'figpath')
        
        saveas(gcf,[OPTS.figpath '.pdf'])
        saveas(gcf,[OPTS.figpath '.fig'])
        disp(['Saved figure to ' OPTS.figpath '.pdf'])
        
    end
    
end

