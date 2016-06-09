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

if FSTD.i == OPTS.nt
    
            nplots = 5;

    
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
        
        nplots = 5;
        skip = floor(OPTS.nt/nplots);
        hold on
        inds = 1:skip:size(DIAG.FSTD.psi,3); 
        
        for i = 1:length(inds)
        
            plot(FSTD.Rint,log10(squeeze(sum(bsxfun(@times,DIAG.FSTD.psi(:,:,inds(i)),FSTD.dA),2)+eps)),'color',colplots(i,:));
        
        end
        
        
        nplots = 10;
        skip = floor(OPTS.nt/nplots);
        llim = floor(log10(1/length(FSTD.Rint)) - 1);
                
        grid on
        box on
        xlabel('Floe Size')
        title('log_{10} of FSD(r)dr')
        ylabel('log10(m^2/m^2)')
        set(gca,'ydir','normal','layer','top','fontname','helvetica','fontsize',14)
        
        set(gca,'ylim',[llim 0])
        
        shading interp
        
        
        
        
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
        
        
        skip = floor(OPTS.nt/nplots);
        hold on
        inds = 1:skip:size(DIAG.FSTD.psi,3); 
        
        for i = 1:length(inds)
        
            plot(FSTD.Hmid,log10(squeeze(sum(bsxfun(@times,DIAG.FSTD.psi(:,:,inds(i)),FSTD.dA),1)+eps)),'color',colplots(i,:));
        
        end
        
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
        
        
    end
    
    %% Plot Ice Volume
    
    
    ax_vol = subplot(349);
    
    if FSTD.i == 1
        
        hold on
        
        if ADVECT.DO
            
            V_init = DIAG.FSTD.Vtot(1);
            
            tauadvect = OPTS.Domainwidth / ADVECT.v1;
            
            V_in = integrate_FSTD(ADVECT.FSTD_in,FSTD.Hmid,FSTD.dA,0); % The total ice volume (area-weighted)
            
            
            % Compute the advective-only solution
            V_adv = V_in + exp(1).^(-FSTD.time/tauadvect)*(V_init - V_in);
            
            plot(FSTD.time/86400,V_adv,'--','color',cols(1,:),'linewidth',2)
            
            
        end
        
        if THERMO.DO
            
            tau_0 = THERMO.fixed_Q(1) / (OPTS.L_f * OPTS.rho_ice);
            
            % Compute the thermo solution
            V_thermo = DIAG.FSTD.Vtot(1) - tau_0 * FSTD.time;
            V_thermo(V_thermo < 0) = 0;
            plot(FSTD.time/86400,V_thermo,'--','color',cols(2,:),'linewidth',2)
            
        end
        
        
        if THERMO.DO && ADVECT.DO
            
            
            % Compute the full solution
            delta = V_in / tauadvect - tau_0;
            V_both = (V_init - tauadvect * delta) * exp(1).^(-FSTD.time/tauadvect) + tauadvect * delta;
            
            plot(FSTD.time/86400,V_both,'-','color',cols(3,:),'linewidth',2)
            
        end
        
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
        
        if ADVECT.DO
            
            c0 = integrate_FSTD(ADVECT.FSTD_in,1,FSTD.dA,0);
            c1 = integrate_FSTD(DIAG.FSTD.psi(:,:,1),1,FSTD.dA,0);
            
            C_adv = c0 + (c1 - c0) * exp(1).^(-FSTD.time/tauadvect);
            
            plot(FSTD.time/86400,C_adv,'--','color',cols(1,:),'linewidth',2)
            c0 = integrate_FSTD(ADVECT.FSTD_in,1,FSTD.dA,1);
            plot(FSTD.time/86400,0*FSTD.time + c0,'--k')
            
        end
        
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
        
        
        if ADVECT.DO
            
            h0 = integrate_FSTD(ADVECT.FSTD_in,FSTD.Hmid,FSTD.dA,0);
            h1 = integrate_FSTD(DIAG.FSTD.psi(:,:,1),FSTD.Hmid,FSTD.dA,0);
            
            
            plot(FSTD.time/86400,0*FSTD.time + h0,'--k')
            
            H_adv = h0 + (h1 - h0) * exp(1).^(-FSTD.time/tauadvect);
            H_adv = H_adv ./ C_adv;
            
            plot(FSTD.time/86400,H_adv,'--','color',cols(1,:),'linewidth',2)
            
        end
        
        
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
        
        
        if ADVECT.DO
            
            n0 = integrate_FSTD(ADVECT.FSTD_in./ (pi * FSTD.meshRmid.^2),1,FSTD.dA,0);
            n1 = integrate_FSTD(DIAG.FSTD.psi(:,:,1)./ (pi * FSTD.meshRmid.^2),1,FSTD.dA,0);
            
            N_adv = n0 + (n1 - n0) * exp(1).^(-FSTD.time/tauadvect);
            
            
            r0num0 = integrate_FSTD(ADVECT.FSTD_in,FSTD.meshRmid,FSTD.dA,0);
            r0num1 = integrate_FSTD(DIAG.FSTD.psi(:,:,1),FSTD.meshRmid,FSTD.dA,0);
            
            r0_adv = r0num0 + (r0num1 - r0num0) * exp(1).^(-FSTD.time/tauadvect);
            r0_adv = r0_adv ./ C_adv;
            
            
            r1num0 = integrate_FSTD(ADVECT.FSTD_in ./ (pi * FSTD.meshRmid.^2),FSTD.meshRmid,FSTD.dA,0);
            r1num1 = integrate_FSTD(DIAG.FSTD.psi(:,:,1) ./ (pi * FSTD.meshRmid.^2),FSTD.meshRmid,FSTD.dA,0);
            
            r1_adv = r1num0 + (r1num1 - r1num0) * exp(1).^(-FSTD.time/tauadvect);
            r1_adv = r1_adv ./ N_adv;
            
            OPTS.r0_pred = plot(FSTD.time/86400,r0_adv,'--','color',cols(1,:),'linewidth',2);
            OPTS.r1_pred = plot(FSTD.time/86400,r1_adv,'--','color',cols(1,:),'linewidth',2);
            
            r0 = ADVECT.FSTD_in ./ (pi * FSTD.meshRmid.^2);
            
            r0 = integrate_FSTD(r0,FSTD.Rmid',FSTD.dA,1);
            
            
            plot(FSTD.time/86400,0*FSTD.time + r0,'--k')
            
        end
        
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
    
    legend(ax_vol,'Predicted - Advect','Predicted - Thermo','Predicted - Both','Calculated');
    legend(ax_mfs,[OPTS.h5 OPTS.h6 OPTS.r0_pred],'By Number','By Area','Predicted  - Advect');
    
    pos = [16 16];
    set(gcf,'windowstyle','normal','position',[0 0 pos],'paperposition',[0 0 pos],'papersize',pos,'units','inches');
    set(gcf,'windowstyle','normal','position',[0 0 pos],'paperposition',[0 0 pos],'papersize',pos,'units','inches');
    saveas(gcf,'~/Dropbox/FSTD/Manuscripts/FSTD-Code-Feedbacks/Figures/Advect-Only/Advect-Only.fig')
    saveas(gcf,'~/Dropbox/FSTD/Manuscripts/FSTD-Code-Feedbacks/Figures/Advect-Only/Advect-Only.pdf')
    
    
end

