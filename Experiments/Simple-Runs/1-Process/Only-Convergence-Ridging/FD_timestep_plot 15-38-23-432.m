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
    
    subplot(221)
    
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
    
    subplot(222)
    
    
    if FSTD.i == 1
        
        concest = DIAG.FSTD.conc(1) + cumsum(EXFORC.nu(:,1) * OPTS.dt);
        % concest = [DIAG.FSTD.conc(1); concest];
        
        plot(FSTD.time/86400,concest,'-','color',cols(6,:),'linewidth',2);
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
        OPTS.h3 = plot(DIAG.FSTD.time(1:FSTD.i)/86400,DIAG.FSTD.conc(1:FSTD.i),'--','color',cols(5,:),'linewidth',2);
        
    end
    
    subplot(223)
    
    if FSTD.i == 1
        
        hest = DIAG.FSTD.Hmean(1) ./ concest;
        plot(FSTD.time/86400,hest,'-','color',cols(6,:),'linewidth',2);
        
        
        hold on
        OPTS.h4 = plot(DIAG.FSTD.time(1:FSTD.i)/86400,DIAG.FSTD.Hmean(1:FSTD.i),'-','color',cols(3,:),'linewidth',2);
        
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
        
        OPTS.h4 = plot(DIAG.FSTD.time(1:FSTD.i)/86400,DIAG.FSTD.Hmean(1:FSTD.i),'--','color',cols(5,:),'linewidth',2);
        
    end
    
    subplot(224)
    
    if FSTD.i == 1
        
        
        OPTS.h5 = plot(DIAG.FSTD.time(1:FSTD.i)/86400,DIAG.FSTD.Rmeanarea(1:FSTD.i),'-','color',cols(3,:),'linewidth',2);
        hold on
        
        OPTS.h6 = plot(DIAG.FSTD.time(1:FSTD.i)/86400,DIAG.FSTD.Rmeannum(1:FSTD.i),'-','color',cols(3,:),'linewidth',2);
        
        grid on
        box on
        xlabel('Time (days)')
        title('Mean Floe Size')
        ylabel('m')
        set(gca,'ydir','normal','layer','top','fontname','helvetica','fontsize',14)
        xlim([0 FSTD.time(end)/86400]);
        legend('By Area','By Number')
        
        
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
    saveas(gcf,'~/Dropbox/FSTD/Manuscripts/FSTD-Code-Feedbacks/Figures/Mech-Only/mech-only.fig')
    saveas(gcf,'~/Dropbox/FSTD/Manuscripts/FSTD-Code-Feedbacks/Figures/Mech-Only/mech-only.pdf')
    
end

drawnow