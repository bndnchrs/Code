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
    
    
    if DIAG.PLOT_FSTD
        
        
        PLOTS.plot_all = figure;
        
        nfig = 1 + 2*DIAG.PLOT_OCEAN + DIAG.PLOT_OC_PROF; 
        pos = [1/nfig .5]; 
        paperpos = [8 4]; 

        set(gcf,'windowstyle','normal','position',[0 0 pos],'paperposition',[0 0 paperpos],'papersize',pos,'units','normalized','paperunits','inches');
        set(gcf,'windowstyle','normal','position',[0 0 pos],'paperposition',[0 0 paperpos],'papersize',pos,'units','normalized','paperunits','inches');
        
        
        PLOTS.ax_FSD = subplot(321);
        PLOTS.ax_ITD = subplot(322);
        PLOTS.ax_vol = subplot(3,4,9);
        PLOTS.ax_mfs = subplot(3,4,12);
        PLOTS.ax_conc = subplot(3,4,10);
        PLOTS.ax_thick = subplot(3,4,11);
        
        colormap(cmap);
        
    end
    
    if DIAG.PLOT_OCEAN
        
        PLOTS.plot_oc = figure;
        colormap(cmap);
        set(gcf,'windowstyle','normal','position',[pos(1) 0 pos],'paperposition',[0 0 paperpos],'papersize',pos,'units','normalized','paperunits','inches');
        set(gcf,'windowstyle','normal','position',[pos(1) 0 pos],'paperposition',[0 0 paperpos],'papersize',pos,'units','normalized','paperunits','inches');
        
        
        if DIAG.PLOT_OC_PROF
            PLOTS.oc_prof = figure;
            colormap(cmap);
            set(gcf,'windowstyle','normal','position',[2*pos(1) 0 pos],'paperposition',[0 0 paperpos],'papersize',pos,'units','normalized','paperunits','inches');
            set(gcf,'windowstyle','normal','position',[2*pos(1) 0 pos],'paperposition',[0 0 paperpos],'papersize',pos,'units','normalized','paperunits','inches');
        end
        
    end
    
    
end

if FSTD.i == OPTS.nt
    %%
    
    
    nplots = 7;
    skip = ceil(OPTS.nt/nplots);
    % Plotting times
    if ~isfield(PLOTS,'plot_inds')
        inds = 1:skip:size(DIAG.FSTD.psi,3);
    else
        inds = PLOTS.plot_inds;
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

if DIAG.PLOT_FSTD && ((DIAG.plot_realtime && mod(FSTD.i,skip_contour) == 1) || FSTD.i == OPTS.nt || skip_contour == 1)
    
    figure(PLOTS.plot_all)
    
    FD_plot_contour_FSDITD;
    
    FD_plot_line_FSDITD;
    
    FD_plot_single_variables;
    
end


%%
if DIAG.PLOT_OCEAN && ((DIAG.plot_realtime && mod(FSTD.i,skip_contour) == 1) || FSTD.i == OPTS.nt || skip_contour == 1)
    
    figure(PLOTS.plot_oc)
    
    
    FD_plot_ocean_stats;
    
    
    drawnow
end

if FSTD.i == OPTS.nt
    
    if DIAG.PLOT_FSTD
        
        figure(PLOTS.plot_all)
        drawnow
        
        %   legend(ax_vol,'Predicted - Advect','Predicted - Thermo','Predicted - Both','Calculated');
        
        legend(PLOTS.ax_mfs,[PLOTS.h5 PLOTS.h6],'By Number','By Area');
        legend(PLOTS.ax_thick,[PLOTS.h4 PLOTS.h42],'Mean','Max');
        
        
        pos = [16 12];
        set(gcf,'windowstyle','normal','position',[0 0 pos],'paperposition',[0 0 pos],'papersize',pos,'units','inches');
        set(gcf,'windowstyle','normal','position',[0 0 pos],'paperposition',[0 0 pos],'papersize',pos,'units','inches');
        
        if isfield(OPTS,'figpath') && OPTS.saveplots
            
            saveas(gcf,[OPTS.figpath '.pdf'])
            saveas(gcf,[OPTS.figpath '.fig'])
            disp(['Saved figure to ' OPTS.figpath '.pdf'])
            
        end
        
    end
    
    if DIAG.PLOT_OCEAN
        
        figure(PLOTS.plot_oc)
        drawnow
        
        pos = [16 12];
        set(gcf,'windowstyle','normal','position',[0 0 pos],'paperposition',[0 0 pos],'papersize',pos,'units','inches');
        set(gcf,'windowstyle','normal','position',[0 0 pos],'paperposition',[0 0 pos],'papersize',pos,'units','inches');
        
        if isfield(OPTS,'figpath') && OPTS.saveplots
            
            saveas(gcf,[OPTS.figpath '-ocean.pdf'])
            saveas(gcf,[OPTS.figpath '-ocean.fig'])
            disp(['Saved figure to ' OPTS.figpath '-ocean.pdf'])
            
        end
        
    end
    
end

