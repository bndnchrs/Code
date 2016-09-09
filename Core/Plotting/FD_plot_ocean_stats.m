cols = [228,26,28
    55,126,184
    77,175,74
    152,78,163
    255,127,0
    155,155,155]/256;


%% Plot Temperatures
subplot(3,2,1)
cla

plot(DIAG.FSTD.time(2:FSTD.i)/86400,DIAG.OCEAN.T(2:FSTD.i),'-','color',cols(5,:),'linewidth',2);
hold on
plot(DIAG.FSTD.time(2:FSTD.i)/86400,DIAG.OCEAN.T_s(2:FSTD.i),'-','color',cols(4,:),'linewidth',2);
% plot(DIAG.FSTD.time(2:FSTD.i)/86400,DIAG.OCEAN.T_a(2:FSTD.i),'--','color',cols(4,:),'linewidth',2);

if FSTD.i == OPTS.nt

    legend({'Mixed Layer Temp','Surface Layer Temp'})

end
hold off

grid on
box on
xlabel('Time (days)')
title('Temps')
ylabel('^\circ C')
set(gca,'ydir','normal','layer','top','fontname','helvetica','fontsize',14)
xlim([0 FSTD.time(end)/86400]);

%% Plot Salinity
subplot(3,2,2)
cla

plot(DIAG.FSTD.time(2:FSTD.i)/86400,DIAG.OCEAN.S(2:FSTD.i),'-','color',cols(5,:),'linewidth',2);

grid on
box on
xlabel('Time (days)')
title('Mixed Layer Salinity')
ylabel('psu')
set(gca,'ydir','normal','layer','top','fontname','helvetica','fontsize',14)
xlim([0 FSTD.time(end)/86400]);


%% Plot Mixed Layer Depth
subplot(3,2,3)
cla

plot(DIAG.FSTD.time(2:FSTD.i)/86400,DIAG.OCEAN.H_ml(2:FSTD.i),'-','color',cols(5,:),'linewidth',2);

grid on
box on
xlabel('Time (days)')
title('Mixed Layer Depth')
ylabel('m')
set(gca,'ydir','normal','layer','top','fontname','helvetica','fontsize',14)
xlim([0 FSTD.time(end)/86400]);

%% Plot Radiative Fluxes
subplot(3,2,4)
cla

plot(DIAG.FSTD.time(2:FSTD.i)/86400,DIAG.OCEAN.QLW_out(2:FSTD.i),'-','color',cols(3,:),'linewidth',2);
hold on
plot(DIAG.FSTD.time(2:FSTD.i)/86400,DIAG.OCEAN.QLW_in(2:FSTD.i),'-','color',cols(4,:),'linewidth',2);
plot(DIAG.FSTD.time(2:FSTD.i)/86400,DIAG.OCEAN.QSW_in(2:FSTD.i),'--','color',cols(5,:),'linewidth',2);

strs = {'LW_out','LW_in','SW_in'}; 

if OCEAN.do_LH
plot(DIAG.FSTD.time(2:FSTD.i)/86400,DIAG.OCEAN.QLH(2:FSTD.i),'-','color',cols(1,:),'linewidth',2);
strs{end+1} = 'LH'; 
end

if OCEAN.do_SH
plot(DIAG.FSTD.time(2:FSTD.i)/86400,DIAG.OCEAN.QSH(2:FSTD.i),'-','color',cols(2,:),'linewidth',2);
strs{end+1} = 'SH';
end

hold off

if FSTD.i == OPTS.nt

    legend(strs)

end

grid on
box on
xlabel('Time (days)')
title('Ocean Surface Heat Fluxes')
ylabel('W/m^2')
set(gca,'ydir','normal','layer','top','fontname','helvetica','fontsize',14)
xlim([0 FSTD.time(end)/86400]);

%% Plot Ocean Fluxes
subplot(3,2,5)
cla

plot(DIAG.FSTD.time(2:FSTD.i)/86400,DIAG.OCEAN.Q_ml_out(2:FSTD.i),'-','color',cols(2,:),'linewidth',2);
hold on

plot(DIAG.FSTD.time(2:FSTD.i)/86400,DIAG.OCEAN.Q_base_mix(2:FSTD.i),'-','color',cols(5,:),'linewidth',2);
plot(DIAG.FSTD.time(2:FSTD.i)/86400,DIAG.OCEAN.Qsurf_ml(2:FSTD.i),'--','color',cols(3,:),'linewidth',2);
plot(DIAG.FSTD.time(2:FSTD.i)/86400,DIAG.OCEAN.Q_mi(2:FSTD.i),'--','color',cols(4,:),'linewidth',2);
plot(DIAG.FSTD.time(2:FSTD.i)/86400,-DIAG.OCEAN.QSW_in(2:FSTD.i),'--','color',cols(1,:),'linewidth',2);

Q_net = DIAG.OCEAN.dV_thermo * (OPTS.L_f * OPTS.rho_ice); 

plot(DIAG.FSTD.time(2:FSTD.i)/86400,Q_net(2:FSTD.i),'-k','linewidth',2);

if FSTD.i == OPTS.nt

legend({'Heat exchange in ML','Heat from below','Heat to Surface','Heat to Ice','Heat from SW','Net from Vol Loss'});

end

grid on
box on
xlabel('Time (days)')
title('Ocean Heat Fluxes')
ylabel('W/m^2')
set(gca,'ydir','normal','layer','top','fontname','helvetica','fontsize',14)
xlim([0 FSTD.time(end)/86400]);

subplot(3,2,6)
cla

plot(DIAG.FSTD.time(2:FSTD.i)/86400,DIAG.OCEAN.S_ml_out(2:FSTD.i),'-','color',cols(1,:),'linewidth',2);
hold on
plot(DIAG.FSTD.time(2:FSTD.i)/86400,DIAG.OCEAN.S_base_mix(2:FSTD.i),'-','color',cols(2,:),'linewidth',2);

if FSTD.i == OPTS.nt
    
    legend({'S change from ice','S mixing at ML base'})

end


grid on
box on
xlabel('Time (days)')
title('Ocean Salinity Fluxes')
ylabel('W/m^2')
set(gca,'ydir','normal','layer','top','fontname','helvetica','fontsize',14)
xlim([0 FSTD.time(end)/86400]);



%%
if FSTD.i == OPTS.nt && DIAG.PLOT_OC_PROF
    
    figure(PLOTS.oc_prof)
    
    Tprof = zeros(OPTS.nt,length(OCEAN.Z));
    Sprof = Tprof;
    
    cmap = [103,0,31
        178,24,43
        214,96,77
        244,165,130
        253,219,199
        247,247,247
        209,229,240
        146,197,222
        67,147,195
        33,102,172
        5,48,97]/256;
    
    clear cmap2
    for i = 1:size(cmap,2)
        cmap2(:,i) = interp(cmap(:,i),1);
    end
    
    cmap2(cmap2 < 0) = 0;
    cmap2(cmap2 > 1) = 1;
    
    for i = 1:OPTS.nt
        
        [~,zval] = find(OCEAN.Z < DIAG.OCEAN.H_ml(i),1,'last');
        
        
        Tprof(i,:) = OCEAN.T_b(OCEAN.Z);
        Tprof(i,1:zval) = DIAG.OCEAN.T(i);
        
        Sprof(i,:) = OCEAN.S_b(OCEAN.Z);
        Sprof(i,1:zval) = DIAG.OCEAN.S(i);
        
        
        
    end
    
    
    subplot(121)
    
    imagesc(FSTD.time/86400,-OCEAN.Z,Tprof');
    shading interp
    grid on
    box on
    set(gca,'ydir','normal','layer','top','fontname','helvetica','fontsize',12)
    
    subplot(122)
    
    imagesc(FSTD.time/86400,-OCEAN.Z,Sprof');
    shading interp
    grid on
    box on
    set(gca,'ydir','normal','layer','top','fontname','helvetica','fontsize',12)
    colormap(flipud(cmap2))
    
end
