if FSTD.i == 1
    close all
    figure; 
end


if mod(FSTD.i,1) == 0
    
    cols = [228,26,28
        55,126,184
        77,175,74
        152,78,163
        255,127,0
        155,155,155]/256;
    
    if FSTD.i == 1
        
        subplot(121)
        tau = EXFORC.Q_oc / (OPTS.L_f * OPTS.rho_ice);
        analytic = DIAG.FSTD.Vtot(1) - tau * FSTD.time;
        analytic(analytic < 0) = 0;
        
        plot(FSTD.time/86400,analytic,'color',cols(6,:),'linewidth',2)
        hold on
        grid on
        box on
        xlabel('Time (days)')
        title('Ice Volume')
        ylabel('m^3/m^2')
        set(gca,'ydir','normal','layer','top','fontname','helvetica','fontsize',14)
        pos = [16 8];
        set(gcf,'windowstyle','normal','position',[0 0 pos],'paperposition',[0 0 pos],'papersize',pos,'units','inches');
        set(gcf,'windowstyle','normal','position',[0 0 pos],'paperposition',[0 0 pos],'papersize',pos,'units','inches');
        xlim([0 FSTD.time(end)/86400]);
        a = get(gca,'ylim');
        a(1) = 0;
        set(gca,'ylim',a);
        OPTS.h1 = plot(DIAG.FSTD.time(1:FSTD.i)/86400,DIAG.FSTD.Vtot(1:FSTD.i),'color',cols(3,:),'linewidth',2);
        
        legend('Predicted','Computed')
        
    else
        
        subplot(121)
        
        
        delete(OPTS.h1)
        
        OPTS.h1 = plot(DIAG.FSTD.time(1:FSTD.i)/86400,DIAG.FSTD.Vtot(1:FSTD.i),'--','color',cols(5,:),'linewidth',2);
        
    end
    
    subplot(222)
    
    
    if FSTD.i == 1
        
        
        hold on
        OPTS.h3 = plot(DIAG.FSTD.time(1:FSTD.i)/86400,DIAG.FSTD.conc(1:FSTD.i),'-','color',cols(5,:),'linewidth',2);
        grid on
        box on
   %     xlabel('Time (days)')
        title('Ice Concentration')
        ylabel('m^2/m^2')
        set(gca,'ydir','normal','layer','top','fontname','helvetica','fontsize',14)
        
        xlim([0 FSTD.time(end)/86400]);
        ylim([0 1])
        
        
    else
        
        
        delete(OPTS.h3)
        OPTS.h3 = plot(DIAG.FSTD.time(1:FSTD.i)/86400,DIAG.FSTD.conc(1:FSTD.i),'-','color',cols(5,:),'linewidth',2);
        
    end
    
    subplot(224)
    
    if FSTD.i == 1
        
        
        
        OPTS.h4 = plot(DIAG.FSTD.time(1:FSTD.i)/86400,DIAG.FSTD.Hmean(1:FSTD.i),'-','color',cols(3,:),'linewidth',2);
        hold on
        
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
    
end

% drawnow

if FSTD.i == OPTS.nt
    
    pos = [8 4];
    set(gcf,'windowstyle','normal','position',[0 0 pos],'paperposition',[0 0 pos],'papersize',pos,'units','inches');
    set(gcf,'windowstyle','normal','position',[0 0 pos],'paperposition',[0 0 pos],'papersize',pos,'units','inches');
    saveas(gcf,'~/Dropbox/FSTD/Manuscripts/FSTD-Code-Feedbacks/Figures/Thermo-Only/thermo-freeze.fig')
    saveas(gcf,'~/Dropbox/FSTD/Manuscripts/FSTD-Code-Feedbacks/Figures/Thermo-Only/thermo-freeze.pdf')
    
end