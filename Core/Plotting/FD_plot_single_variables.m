cols = [228,26,28
    55,126,184
    77,175,74
    152,78,163
    255,127,0
    155,155,155]/256;


%% Plot Ice Volume


set(gcf,'currentaxes',PLOTS.ax_vol);

if FSTD.i == 1
    
    
    
    hold on
    
    if ADVECT.DO
        
        V_init = DIAG.FSTD.Vtot(1);
        
        tauadvect = OPTS.Domainwidth / ADVECT.v1;
        
        V_in = integrate_FSTD(ADVECT.FSTD_in,FSTD.Hmid_i,FSTD.dA_i,0); % The total ice volume (area-weighted)
        
        
        % Compute the advective-only solution
        V_adv = V_in + exp(1).^(-FSTD.time/tauadvect)*(V_init - V_in);
        
        plot(FSTD.time/86400,V_adv,'-','color',cols(1,:),'linewidth',2)
        
        
    end
    
    if THERMO.DO && THERMO.fixQ
        
        tau_0 = THERMO.fixed_Q(1) / (OPTS.L_f * OPTS.rho_ice);
        
        % Compute the thermo solution
        V_thermo = DIAG.FSTD.Vtot(1) - tau_0 * FSTD.time;
        V_thermo(V_thermo < 0) = 0;
        plot(FSTD.time/86400,V_thermo,'-','color',cols(2,:),'linewidth',2)
        
    end
    
    if MECH.DO
        
        vest = DIAG.FSTD.Vtot(1);
        
        plot(FSTD.time/86400,0*FSTD.time + vest,'-','color',cols(4,:),'linewidth',2)
        
    end
    
    if THERMO.DO && ADVECT.DO && THERMO.fixQ
        
        
        % Compute the full solution
        delta = V_in / tauadvect - tau_0;
        V_both = (V_init - tauadvect * delta) * exp(1).^(-FSTD.time/tauadvect) + tauadvect * delta;
        
        plot(FSTD.time/86400,V_both,'-','color',cols(3,:),'linewidth',2)
        
    end
    
    PLOTS.h1 = plot(DIAG.FSTD.time(1:FSTD.i)/86400,DIAG.FSTD.Vtot(1:FSTD.i),'--','color',cols(5,:),'linewidth',2);
    PLOTS.h2 = scatter(DIAG.FSTD.time(FSTD.i)/86400,DIAG.FSTD.Vtot(FSTD.i),200,'filled','markerfacecolor',cols(5,:));
    
    
else
    
    if isfield(PLOTS,'h1')
        
        delete(PLOTS.h1)
    end
    
    if isfield(PLOTS,'h2')
        delete(PLOTS.h2)
    end
    
    PLOTS.h1 = plot(DIAG.FSTD.time(1:FSTD.i)/86400,DIAG.FSTD.Vtot(1:FSTD.i),'--','color',cols(5,:),'linewidth',2);
    hold on
    PLOTS.h2 = scatter(DIAG.FSTD.time(FSTD.i)/86400,DIAG.FSTD.Vtot(FSTD.i),200,'filled','markerfacecolor',cols(5,:));
    hold off
    
end

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


%% Plot Ice Concentration
set(gcf,'currentaxes',PLOTS.ax_conc);



if FSTD.i == 1
    
    if ADVECT.DO
        
        c0 = integrate_FSTD(ADVECT.FSTD_in,1,FSTD.dA,0);
        c1 = integrate_FSTD(DIAG.FSTD.psi(:,:,1),1,FSTD.dA,0);
        
        C_adv = c0 + (c1 - c0) * exp(1).^(-FSTD.time/tauadvect);
        
        plot(FSTD.time/86400,C_adv,'-','color',cols(1,:),'linewidth',2)
        c0 = integrate_FSTD(ADVECT.FSTD_in,1,FSTD.dA,1);
        plot(FSTD.time/86400,0*FSTD.time + c0,'--k')
        
    end
    
    if MECH.DO
        
        concest = DIAG.FSTD.conc(1) + cumsum(EXFORC.nu(:,1) * OPTS.dt);
        
        plot(FSTD.time/86400,concest,'-','color',cols(6,:),'linewidth',2)
        
    end
    
    hold on
    PLOTS.h3 = plot(DIAG.FSTD.time(1:FSTD.i)/86400,DIAG.FSTD.conc(1:FSTD.i),'-','color',cols(5,:),'linewidth',2);
    
    
else
    if isfield(PLOTS,'h3')
        
        delete(PLOTS.h3)
        
    end
    
    PLOTS.h3 = plot(DIAG.FSTD.time(1:FSTD.i)/86400,DIAG.FSTD.conc(1:FSTD.i),'--','color',cols(5,:),'linewidth',2);
    
    
end

grid on
box on
xlabel('Time (days)')
title('Ice Concentration')
ylabel('m^2/m^2')
set(gca,'ydir','normal','layer','top','fontname','helvetica','fontsize',14)
xlim([0 FSTD.time(end)/86400]);
ylim([0 1])

%% Plot Ice Thickness


set(gcf,'currentaxes',PLOTS.ax_thick)

if FSTD.i == 1
    
    
    if ADVECT.DO
        
        h0 = integrate_FSTD(ADVECT.FSTD_in,FSTD.Hmid,FSTD.dA,0);
        h1 = integrate_FSTD(DIAG.FSTD.psi(:,:,1),FSTD.Hmid,FSTD.dA,0);
        
        
        plot(FSTD.time/86400,0*FSTD.time + h0,'--k')
        
        H_adv = h0 + (h1 - h0) * exp(1).^(-FSTD.time/tauadvect);
        H_adv = H_adv ./ C_adv;
        
        plot(FSTD.time/86400,H_adv,'-','color',cols(1,:),'linewidth',2)
        
    end
    
    if MECH.DO
        
        hest = DIAG.FSTD.Vtot(1) ./ concest;
        plot(FSTD.time/86400,hest,'-','color',cols(6,:),'linewidth',2);
        
    end
    
    
    hold on
    PLOTS.h4 = plot(DIAG.FSTD.time(1:FSTD.i)/86400,DIAG.FSTD.Hmean(1:FSTD.i),'-','color',cols(5,:),'linewidth',2);
    PLOTS.h42 = plot(DIAG.FSTD.time(1:FSTD.i)/86400,DIAG.FSTD.Hmax(1:FSTD.i),'-','color',cols(4,:),'linewidth',2);
    
    
else
    
    if isfield(PLOTS,'h4')
        
        
        delete(PLOTS.h4)
        
    end
    
    if isfield(PLOTS,'h42')
        
        delete(PLOTS.h42)
        
    end
    
    PLOTS.h4 = plot(DIAG.FSTD.time(1:FSTD.i)/86400,DIAG.FSTD.Hmean(1:FSTD.i),'--','color',cols(5,:),'linewidth',2);
    hold on
    PLOTS.h42 = plot(DIAG.FSTD.time(1:FSTD.i)/86400,DIAG.FSTD.Hmax(1:FSTD.i),'--','color',cols(4,:),'linewidth',2);
    hold off
    
    
    
    
end

ylim auto
a = get(gca,'ylim');
a(1) = 0;
set(gca,'ylim',a);
grid on
box on
xlabel('Time (days)')
title('Ice Thickness')
ylabel('m')
set(gca,'ydir','normal','layer','top','fontname','helvetica','fontsize',14)
xlim([0 FSTD.time(end)/86400]);

%% Plot Mean Floe Sizes

set(gcf,'currentaxes',PLOTS.ax_mfs)

if FSTD.i == 1
    
    
    
    
    
    PLOTS.h5 = plot(DIAG.FSTD.time(1:FSTD.i)/86400,DIAG.FSTD.Rmeannum(1:FSTD.i),'-','color',cols(3,:),'linewidth',2);
    hold on
    PLOTS.h6 = plot(DIAG.FSTD.time(1:FSTD.i)/86400,DIAG.FSTD.Rmeanarea(1:FSTD.i),'-','color',cols(5,:),'linewidth',2);
    
    legend('By Number','By Area')
    
    if WAVES.DO
        
        [a,b] = max(WAVES.spec);
        
        % This is R which is half the wavelength of the peak period
        R_longterm = FSTD.Rmid(b);
        plot(FSTD.time/86400,0*FSTD.time + R_longterm,'-k','linewidth',2)
        
        r_0 = DIAG.FSTD.Rmeanarea(1);
        
        tau_waves = WAVES.tau;
        
        rbar = R_longterm + exp(1).^(-FSTD.time/tau_waves)*(r_0 - R_longterm);
        
        % plot(FSTD.time/86400,rbar,'-','color',cols(5,:),'linewidth',2)
        
    end
    
    
    
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
        
        PLOTS.r0_pred = plot(FSTD.time/86400,r0_adv,'-','color',cols(1,:),'linewidth',2);
        PLOTS.r1_pred = plot(FSTD.time/86400,r1_adv,'-','color',cols(1,:),'linewidth',2);
        
        r0 = ADVECT.FSTD_in ./ (pi * FSTD.meshRmid.^2);
        
        r0 = integrate_FSTD(r0,FSTD.Rmid',FSTD.dA,1);
        
        
        plot(FSTD.time/86400,0*FSTD.time + r0,'--k')
        
    end
    
    
else
    
    if isfield(PLOTS,'h5')
        
        delete(PLOTS.h5)
        
    end
    
    if isfield(PLOTS,'h6')
        
        delete(PLOTS.h6)
        
    end
    
    PLOTS.h5 = plot(DIAG.FSTD.time(1:FSTD.i)/86400,DIAG.FSTD.Rmeannum(1:FSTD.i),'--','color',cols(3,:),'linewidth',2);
    hold on
    PLOTS.h6 = plot(DIAG.FSTD.time(1:FSTD.i)/86400,DIAG.FSTD.Rmeanarea(1:FSTD.i),'--','color',cols(5,:),'linewidth',2);
    hold off
    
    
    
    
    
end
ylim auto
a = get(gca,'ylim');
a(1) = 0;
set(gca,'ylim',a);
grid on
box on
xlabel('Time (days)')
title('Mean Floe Size')
ylabel('m')
set(gca,'ydir','normal','layer','top','fontname','helvetica','fontsize',14)
xlim([0 FSTD.time(end)/86400]);


