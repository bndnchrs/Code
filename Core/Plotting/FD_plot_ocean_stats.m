cols = [228,26,28
    55,126,184
    77,175,74
    152,78,163
    255,127,0
    155,155,155]/256;


%% Plot Ice Volume


subplot(2,2,1)

plot(DIAG.FSTD.time(2:FSTD.i)/86400,DIAG.OCEAN.T(2:FSTD.i),'--','color',cols(5,:),'linewidth',2);
hold on
plot(DIAG.FSTD.time(2:FSTD.i)/86400,DIAG.OCEAN.T_s(2:FSTD.i),'--','color',cols(4,:),'linewidth',2);
hold off

grid on
box on
xlabel('Time (days)')
title('Ocean ML Temp')
ylabel('^\circ C')
set(gca,'ydir','normal','layer','top','fontname','helvetica','fontsize',14)
xlim([0 FSTD.time(end)/86400]);
    
subplot(2,2,2)

plot(DIAG.FSTD.time(2:FSTD.i)/86400,DIAG.OCEAN.S(2:FSTD.i),'--','color',cols(5,:),'linewidth',2);

grid on
box on
xlabel('Time (days)')
title('Mixed Layer Salinity')
ylabel('psu')
set(gca,'ydir','normal','layer','top','fontname','helvetica','fontsize',14)
xlim([0 FSTD.time(end)/86400]);

subplot(2,2,3)

plot(DIAG.FSTD.time(2:FSTD.i)/86400,DIAG.OCEAN.H_ml(2:FSTD.i),'--','color',cols(5,:),'linewidth',2);

grid on
box on
xlabel('Time (days)')
title('Mixed Layer Depth')
ylabel('m')
set(gca,'ydir','normal','layer','top','fontname','helvetica','fontsize',14)
xlim([0 FSTD.time(end)/86400]);

subplot(2,2,4)

plot(DIAG.FSTD.time(2:FSTD.i)/86400,DIAG.OCEAN.QLH(2:FSTD.i),'--','color',cols(1,:),'linewidth',2);
hold on
plot(DIAG.FSTD.time(2:FSTD.i)/86400,DIAG.OCEAN.QSH(2:FSTD.i),'--','color',cols(2,:),'linewidth',2);
plot(DIAG.FSTD.time(2:FSTD.i)/86400,DIAG.OCEAN.QLW_out(2:FSTD.i),'--','color',cols(3,:),'linewidth',2);
plot(DIAG.FSTD.time(2:FSTD.i)/86400,DIAG.OCEAN.QLW_in(2:FSTD.i),'--','color',cols(4,:),'linewidth',2);
plot(DIAG.FSTD.time(2:FSTD.i)/86400,DIAG.OCEAN.QSW_in(2:FSTD.i),'--','color',cols(5,:),'linewidth',2);
hold off

legend('LH','SH','LW_out','LW_in','SW_in')
grid on
box on
xlabel('Time (days)')
title('Ocean Surface Heat Fluxes')
ylabel('W/m^2')
set(gca,'ydir','normal','layer','top','fontname','helvetica','fontsize',14)
xlim([0 FSTD.time(end)/86400]);
