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
    
    
    subplot(211)
    
    if FSTD.i == 1
        
        
        
        OPTS.p1 = pcolor([0 FSTD.time/86400],FSTD.Rint,log10(squeeze(sum(bsxfun(@times,DIAG.FSTD.psi,FSTD.dA),2)+eps)));
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
        
        OPTS.p1 = pcolor([0 FSTD.time/86400],FSTD.Rint,log10(squeeze(sum(bsxfun(@times,DIAG.FSTD.psi,FSTD.dA),2)+eps)));
        shading interp
        set(gca,'clim',[-2 0])
        
        if FSTD.i == OPTS.nt
            
            
            [a,b] = max(WAVES.spec);
            
            % This is R which is half the wavelength of the peak period
            R_longterm = FSTD.Rmid(b);
            plot(FSTD.time/86400,0*FSTD.time + R_longterm,'-k','linewidth',2)
            
        end
        
    end
    
    subplot(223)
    
    if FSTD.i == 1
        
        vest = DIAG.FSTD.Vtot(1);
        plot(FSTD.time/86400,0*FSTD.time + vest,'-','color',cols(6,:),'linewidth',2);
        hold on
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
        
        OPTS.h1 = plot(DIAG.FSTD.time(1:FSTD.i)/86400,DIAG.FSTD.Vtot(1:FSTD.i),'color',cols(3,:),'linewidth',2);
        
        legend('Predicted')
        
    else
        
        
        
        delete(OPTS.h1)
        
        OPTS.h1 = plot(DIAG.FSTD.time(1:FSTD.i)/86400,DIAG.FSTD.Vtot(1:FSTD.i),'--','color',cols(5,:),'linewidth',2);
        
    end
    
    
    
    subplot(224)
    
    if FSTD.i == 1
        
        [a,b] = max(WAVES.spec);
        
        % This is R which is half the wavelength of the peak period
        R_longterm = FSTD.Rmid(b);
        % This is the initial mean floe size
        R_init = DIAG.FSTD.Rmeanarea(1);
        
        R_pred = R_longterm + exp(1).^(-FSTD.time/WAVES.tau)*(R_init - R_longterm);
        plot(FSTD.time/86400,R_pred,'-','color',cols(6,:),'linewidth',2);
        
        hold on
        
        OPTS.h5 = plot(DIAG.FSTD.time(1:FSTD.i)/86400,DIAG.FSTD.Rmeanarea(1:FSTD.i),'-','color',cols(3,:),'linewidth',2);
        
        OPTS.h6 = plot(DIAG.FSTD.time(1:FSTD.i)/86400,DIAG.FSTD.Rmeannum(1:FSTD.i),'-','color',cols(3,:),'linewidth',2);
        
        grid on
        box on
        xlabel('Time (days)')
        title('Mean Floe Size')
        ylabel('m')
        set(gca,'ydir','normal','layer','top','fontname','helvetica','fontsize',14)
        xlim([0 FSTD.time(end)/86400]);
        legend('Decay to \lambda_max/2','By Area','By Number')
        a = get(gca,'ylim');
        a(1) = 0;
        set(gca,'ylim',a);
        
    else
        
        delete(OPTS.h5)
        delete(OPTS.h6)
        
        OPTS.h5 = plot(DIAG.FSTD.time(1:FSTD.i)/86400,DIAG.FSTD.Rmeanarea(1:FSTD.i),'-','color',cols(5,:),'linewidth',2);
        OPTS.h6 = plot(DIAG.FSTD.time(1:FSTD.i)/86400,DIAG.FSTD.Rmeannum(1:FSTD.i),'-','color',cols(3,:),'linewidth',2);
        
    end
    
end

% drawnow

if FSTD.i == OPTS.nt
    
    pos = [8 4];
    set(gcf,'windowstyle','normal','position',[0 0 pos],'paperposition',[0 0 pos],'papersize',pos,'units','inches');
    set(gcf,'windowstyle','normal','position',[0 0 pos],'paperposition',[0 0 pos],'papersize',pos,'units','inches');
    saveas(gcf,'~/Dropbox/FSTD/Manuscripts/FSTD-Code-Feedbacks/Figures/Waves-Only/waves-only.fig')
    saveas(gcf,'~/Dropbox/FSTD/Manuscripts/FSTD-Code-Feedbacks/Figures/Waves-Only/waves-only.pdf')
    
end

drawnow