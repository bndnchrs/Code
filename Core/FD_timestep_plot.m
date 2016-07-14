%% FD_timestep_plot
% Plotting routing for simple-thermo-advect run

npcolor = 50;
skip_contour = ceil(OPTS.nt/npcolor); % Normally, we don't sub-sample, but we will for making figures

if FSTD.i == 1
    
    addpath([OPTS.path_of_code 'Core/Plotting/']);
    
    close all
    
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
    
    if DIAG.DO_PLOT_FSTD
        
        
        OPTS.plot_all = figure;
        OPTS.ax_FSD = subplot(321);
        OPTS.ax_ITD = subplot(322);
        OPTS.ax_vol = subplot(3,4,9);
        OPTS.ax_mfs = subplot(3,4,12);
        OPTS.ax_conc = subplot(3,4,10);
        OPTS.ax_thick = subplot(3,4,11);
        
    else
        
        if DIAG.DO_PLOT_OCEAN
            
            OPTS.plot_oc = figure;
            
        end
    end
    
end

if FSTD.i == OPTS.nt
    %%
    
    
    nplots = 7;
    skip = ceil(OPTS.nt/nplots);
    % Plotting times
    if ~isfield(OPTS,'plot_inds')
        inds = 1:skip:size(DIAG.FSTD.psi,3);
    else
        inds = OPTS.plot_inds;
    end
    
    cplots = [103,0,31
        247,251,255
        222,235,247
        198,219,239
        158,202,225
        107,174,214
        66,146,198
        33,113,181
        8,81,156
        8,48,107]/256;
    
    
    
    div = size(cplots,1)/nplots;
    
    for i = 1:3
        colplots(:,i) = downsample(cplots(:,i),floor(div));
    end
    
    
    
    colplots(colplots > 1) = 1;
    colplots(colplots < 0) = 0;
    
end

%%

if DIAG.DO_PLOT_FSTD && (mod(FSTD.i,skip_contour) == 1 || FSTD.i == OPTS.nt || skip_contour == 1)
    
    figure(OPTS.plot_all)
    
    FD_plot_contour_FSDITD;
    
    FD_plot_line_FSDITD;
    
    FD_plot_single_variables;
    
    if DIAG.plot_realtime
        
        drawnow
        
    end
    
end

if DIAG.DO_PLOT_OCEAN && (mod(FSTD.i,skip_contour) == 1 || FSTD.i == OPTS.nt || skip_contour == 1)
    
   figure(OPTS.plot_oc)

    
    FD_plot_ocean_stats;
    
    if DIAG.plot_realtime
        
        drawnow
        
    end
    
end

if FSTD.i == OPTS.nt
    
    if DIAG.DO_PLOT_FSTD
        
        figure(OPTS.plot_all)
        drawnow
        
        %   legend(ax_vol,'Predicted - Advect','Predicted - Thermo','Predicted - Both','Calculated');
        
        legend(OPTS.ax_mfs,[OPTS.h5 OPTS.h6],'By Number','By Area');
        legend(OPTS.ax_thick,[OPTS.h4 OPTS.h42],'Mean','Max');
        
        
        pos = [16 16];
        set(gcf,'windowstyle','normal','position',[0 0 pos],'paperposition',[0 0 pos],'papersize',pos,'units','inches');
        set(gcf,'windowstyle','normal','position',[0 0 pos],'paperposition',[0 0 pos],'papersize',pos,'units','inches');
        
        if isfield(OPTS,'figpath')
            
            saveas(gcf,[OPTS.figpath '.pdf'])
            saveas(gcf,[OPTS.figpath '.fig'])
            disp(['Saved figure to ' OPTS.figpath '.pdf'])
            
        end
        
    end
    
    if DIAG.DO_PLOT_OCEAN
        
        figure(OPTS.plot_oc)
        drawnow
        
        pos = [16 16];
        set(gcf,'windowstyle','normal','position',[0 0 pos],'paperposition',[0 0 pos],'papersize',pos,'units','inches');
        set(gcf,'windowstyle','normal','position',[0 0 pos],'paperposition',[0 0 pos],'papersize',pos,'units','inches');
        
        if isfield(OPTS,'figpath')
            
            saveas(gcf,[OPTS.figpath '-ocean.pdf'])
            saveas(gcf,[OPTS.figpath '-ocean.fig'])
            disp(['Saved figure to ' OPTS.figpath '-ocean.pdf'])
            
        end
        
    end
    
end

