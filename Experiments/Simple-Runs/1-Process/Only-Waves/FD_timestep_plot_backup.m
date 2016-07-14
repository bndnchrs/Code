%% FD_timestep_plot
% Plotting routing for simple-thermo-advect run

if FSTD.i == 1
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
    
end

if mod(FSTD.i,1) == 0
    
    
    cols = [228,26,28
        55,126,184
        77,175,74
        152,78,163
        255,127,0
        155,155,155]/256;
    
    %% Plot FSD
    ax_FSD = subplot(321);
    
    if FSTD.i == 1
        
        
        
        OPTS.p1 = pcolor([0 FSTD.time/86400],FSTD.Rint,log10(squeeze(sum(bsxfun(@times,DIAG.FSTD.psi,DIAG.FSTD.dA),2)+eps)));
        hold on
        grid on
        box on
        ylabel('Floe Size')
        title('log_{10} of FSD(r)dr')
        xlabel('Time')
        set(gca,'ydir','normal','layer','top','fontname','helvetica','fontsize',14)
        
        ylim([FSTD.Rint(1) FSTD.Rint(end)]);
        
        
        set(gca,'clim',[-2 0])
        shading interp
        colorbar
        
    else
        
        
        
        delete(OPTS.p1)
        
        OPTS.p1 = pcolor([0 FSTD.time/86400],FSTD.Rint,log10(squeeze(sum(bsxfun(@times,DIAG.FSTD.psi,DIAG.FSTD.dA),2)+eps)));
        shading interp
        set(gca,'clim',[-4 0])
        
        if FSTD.i == OPTS.nt && WAVES.DO
            
            
            [a,b] = max(WAVES.spec);
            
            % This is R which is half the wavelength of the peak period
            R_longterm = FSTD.Rmid(b);
            plot(FSTD.time/86400,0*FSTD.time + R_longterm,'-k','linewidth',2)
            
        end
        
        
        if FSTD.i == OPTS.nt
            
            imcontour([0 FSTD.time/86400],FSTD.Rint,log10(squeeze(sum(bsxfun(@times,DIAG.FSTD.psi,DIAG.FSTD.dA),2)+eps)),[-4:.5:0],'--k','showtext','on');
            set(gca,'clim',[-4 0])
            
        end
        
        
        
    end
    
    %% Plot ITD
    ax_ITD = subplot(322);
    if length(FSTD.Hmid) > 1
        if FSTD.i == 1
            
            
            
            OPTS.p2 = pcolor([0 FSTD.time/86400],FSTD.Hmid,log10(squeeze(sum(bsxfun(@times,DIAG.FSTD.psi,DIAG.FSTD.dA),1)+eps)));
            
            hold on
            grid on
            box on
            ylabel('Floe Size')
            title('log_{10} of ITD(h)dh')
            xlabel('Time')
            set(gca,'ydir','normal','layer','top','fontname','helvetica','fontsize',14)
            
            ylim([FSTD.Hmid(1) FSTD.Hmid(end)]);
            
            
            set(gca,'clim',[-2 0])
            shading interp
            colorbar
            
        else
            
            
            
            delete(OPTS.p2)
            
            OPTS.p2 = pcolor([0 FSTD.time/86400],FSTD.Hmid,log10(squeeze(sum(bsxfun(@times,DIAG.FSTD.psi,DIAG.FSTD.dA),1)+eps)));
            shading interp
            set(gca,'clim',[-4 0])
            
            
            if FSTD.i == OPTS.nt
                
                imcontour([0 FSTD.time/86400],FSTD.Hmid,log10(squeeze(sum(bsxfun(@times,DIAG.FSTD.psi,DIAG.FSTD.dA),1)+eps)),[-4:.5:0],'--k','showtext','on');
                set(gca,'clim',[-4 0])
                
            end
        end
    end
    
    %% Plot FSd as line plots
    
    ax_FSD = subplot(323);
    
    if FSTD.i == OPTS.nt
        
        
        nplots = 10;
        skip = floor(OPTS.nt/nplots);
        llim = floor(log10(1/length(FSTD.Rint)) - 1);
        
        OPTS.p1 = plot(FSTD.Rint,log10(squeeze(sum(bsxfun(@times,DIAG.FSTD.psi(:,:,1:skip:end),FSTD.dA),2)+eps)));
        
        hold on
        grid on
        box on
        xlabel('Floe Size')
        title('log_{10} of FSD(r)dr')
        ylabel('log10(m^2/m^2)')
        set(gca,'ydir','normal','layer','top','fontname','helvetica','fontsize',14)
        
        set(gca,'ylim',[llim 0])
        
        shading interp
        colorbar
        
        
        
        
        if WAVES.DO
            
            
            [a,b] = max(WAVES.spec);
            
            % This is R which is half the wavelength of the peak period
            R_longterm = FSTD.Rmid(b);
            scatter(R_longterm,.5,200,'k')
            
        end
        
    end
    
    %% Plot ITD as line plots
    
    ax_FSD = subplot(324);
    
    if FSTD.i == OPTS.nt
        
        
        nplots = 10;
        skip = floor(OPTS.nt/nplots);
        
        OPTS.p1 = plot(FSTD.Hmid,log10(squeeze(sum(bsxfun(@times,DIAG.FSTD.psi(:,:,1:skip:end),FSTD.dA),1)+eps)));
        llim = floor(log10(1/length(FSTD.Hmid)) - 1);
        
        hold on
        grid on
        box on
        xlabel('Floe Size')
        title('log_{10} of ITD(h)dh')
        ylabel('log10(m^2/m^2)')
        set(gca,'ydir','normal','layer','top','fontname','helvetica','fontsize',14)
        
        set(gca,'ylim',[llim 0])
        
        shading interp
        colorbar
        
        
    end
    
    %% Plot Ice Volume
    
    
    ax_vol = subplot(349);
    
    if FSTD.i == 1
        
        hold on

        OPTS.h1 = plot(DIAG.FSTD.time(1:FSTD.i)/86400,DIAG.FSTD.Vtot(1:FSTD.i),'--','color',cols(5,:),'linewidth',2);
        OPTS.h2 = scatter(DIAG.FSTD.time(FSTD.i)/86400,DIAG.FSTD.Vtot(FSTD.i),200,'filled','markerfacecolor',cols(5,:));
        
        grid on
        box on
        xlabel('Time (days)')
        title('Ice Volume')
        ylabel('m^3/m^2')
        set(gca,'ydir','normal','layer','top','fontname','helvetica','fontsize',14)
        xlim([0 FSTD.time(end)/86400]);
        a = get(gca,'ylim');
        a(1) = 0;
        set(gca,'ylim',a);
    else
        
        delete(OPTS.h1)
        delete(OPTS.h2)
        OPTS.h1 = plot(DIAG.FSTD.time(1:FSTD.i)/86400,DIAG.FSTD.Vtot(1:FSTD.i),'--','color',cols(5,:),'linewidth',2);
        OPTS.h2 = scatter(DIAG.FSTD.time(FSTD.i)/86400,DIAG.FSTD.Vtot(FSTD.i),200,'filled','markerfacecolor',cols(5,:));
        
    end
    
    %% Plot Ice Concentration
    
    
    ax_conc = subplot(3,4,10);
    
    if FSTD.i == 1
        
        hold on
        OPTS.h3 = plot(DIAG.FSTD.time(1:FSTD.i)/86400,DIAG.FSTD.conc(1:FSTD.i),'-','color',cols(5,:),'linewidth',2);
        grid on
        box on
        xlabel('Time (days)')
        title('Ice Concentration')
        ylabel('m^2/m^2')
        set(gca,'ydir','normal','layer','top','fontname','helvetica','fontsize',14)
        xlim([0 FSTD.time(end)/86400]);
        ylim([0 1])
        
    else
        
        delete(OPTS.h3)
        OPTS.h3 = plot(DIAG.FSTD.time(1:FSTD.i)/86400,DIAG.FSTD.conc(1:FSTD.i),'-','color',cols(5,:),'linewidth',2);
        
    end
    %% Plot Ice Thickness
    
    
    ax_thick = subplot(3,4,11);
    
    if FSTD.i == 1
        
        hold on
        OPTS.h4 = plot(DIAG.FSTD.time(1:FSTD.i)/86400,DIAG.FSTD.Hmean(1:FSTD.i),'-','color',cols(5,:),'linewidth',2);
        
        grid on
        box on
        xlabel('Time (days)')
        title('Mean Ice Thickness')
        ylabel('m')
        set(gca,'ydir','normal','layer','top','fontname','helvetica','fontsize',14)
        xlim([0 FSTD.time(end)/86400]);
        a = get(gca,'ylim');
        a(1) = 0;
        set(gca,'ylim',a);
        
    else
        
        delete(OPTS.h4)
        OPTS.h4 = plot(DIAG.FSTD.time(1:FSTD.i)/86400,DIAG.FSTD.Hmean(1:FSTD.i),'-','color',cols(5,:),'linewidth',2);
        
    end
    
    %% Plot Mean Floe Sizes
    
    ax_mfs = subplot(3,4,12);
    
    if FSTD.i == 1
        
        
        OPTS.h5 = plot(DIAG.FSTD.time(1:FSTD.i)/86400,DIAG.FSTD.Rmeannum(1:FSTD.i),'-','color',cols(3,:),'linewidth',2);
        hold on
        OPTS.h6 = plot(DIAG.FSTD.time(1:FSTD.i)/86400,DIAG.FSTD.Rmeanarea(1:FSTD.i),'-','color',cols(5,:),'linewidth',2);
        
        legend('By Number','By Area')
            
        grid on
        box on
        xlabel('Time (days)')
        title('Mean Floe Size')
        ylabel('m')
        set(gca,'ydir','normal','layer','top','fontname','helvetica','fontsize',14)
        legend('By Number','By Area')
        
        xlim([0 FSTD.time(end)/86400]);
        a = get(gca,'ylim');
        a(1) = 0;
        set(gca,'ylim',a);
    else
        
        delete(OPTS.h5)
        delete(OPTS.h6)
        
        OPTS.h5 = plot(DIAG.FSTD.time(1:FSTD.i)/86400,DIAG.FSTD.Rmeannum(1:FSTD.i),'-','color',cols(3,:),'linewidth',2);
        OPTS.h6 = plot(DIAG.FSTD.time(1:FSTD.i)/86400,DIAG.FSTD.Rmeanarea(1:FSTD.i),'-','color',cols(5,:),'linewidth',2);
        
    end
    
    
    
end

drawnow

if FSTD.i == OPTS.nt
       
    pos = [16 16];
    set(gcf,'windowstyle','normal','position',[0 0 pos],'paperposition',[0 0 pos],'papersize',pos,'units','inches');
    set(gcf,'windowstyle','normal','position',[0 0 pos],'paperposition',[0 0 pos],'papersize',pos,'units','inches');
    saveas(gcf,'~/Dropbox/FSTD/Manuscripts/FSTD-Code-Feedbacks/Figures/Melt-Advect/Melt-Advect.fig')
    saveas(gcf,'~/Dropbox/FSTD/Manuscripts/FSTD-Code-Feedbacks/Figures/Melt-Advect/Melt-Advect.pdf')
    
    %     figure
    %     plot(DIAG.FSTD.time(1:end-1)/86400,DIAG.THERMO.dc_tot(2:end),'-','color',cols(3,:),'linewidth',2)
    %     hold on
    %     plot(DIAG.FSTD.time(1:end-1)/86400,DIAG.ADVECT.dc_adv(2:end),'-','color',cols(2,:),'linewidth',2)
    %     grid on
    %     box on
    %     xlabel('Time (days)')
    %     a = get(gca,'ylim');
    %     set(gca,'ylim',[-1 1] * max(abs(a)))
    %     legend('Thermodynamics','Advection')
    %     ylabel('1/s')
    %     xlabel('Time (days)')
    %     set(gca,'ydir','normal','layer','top','fontname','helvetica','fontsize',14)
    
end