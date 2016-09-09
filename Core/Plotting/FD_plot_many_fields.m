nplots = 4;
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

FD_plot_line_FSDITD;


axis(PLOTS.ax_ITD_line)
figure; subplot(233); get(gca,'position'); close; set(PLOTS.ax_ITD_line,'position',ans)
ylim(PLOTS.ax_ITD_line,[1e-3 1])
xlim(PLOTS.ax_ITD_line,[.1 3.5])
set(PLOTS.ax_ITD_line,'ydir','normal','layer','top','fontname','helvetica','fontsize',12)

axis(PLOTS.ax_FSD_line)
figure; subplot(232); get(gca,'position'); close; set(PLOTS.ax_FSD_line,'position',ans)
set(PLOTS.ax_FSD_line,'xtick',[1 10 100 1000])
ylim(PLOTS.ax_FSD_line,[3e-4 .1])
xlim(PLOTS.ax_FSD_line,[.75 1600])
set(PLOTS.ax_FSD_line,'ydir','normal','layer','top','fontname','helvetica','fontsize',12)


subplot(231)
plot(FSTD.time/86400,DIAG.FSTD.Vtot(2:end)/DIAG.FSTD.Vtot(1),'linewidth',1)
hold on
plot(FSTD.time/86400,DIAG.FSTD.conc(2:end)/DIAG.FSTD.conc(1),'linewidth',1)
plot(FSTD.time/86400,DIAG.FSTD.Hmean(2:end)/DIAG.FSTD.Hmean(1),'linewidth',1)
hold off
grid on
box on
set(gca,'ydir','normal','layer','top','fontname','helvetica','fontsize',12)
xlim([FSTD.time(1) FSTD.time(end)]/86400)
legend('Volume','Area','Thickness')
ylabel('% of initial')
title('Ice Variables')


subplot(234)
[Ax,h1,h2] = plotyy(FSTD.time/86400,DIAG.OCEAN.T(2:end),FSTD.time/86400,DIAG.OCEAN.S(2:end));
hold on
h3 = plot(FSTD.time/86400,DIAG.OCEAN.T_s(2:end),'--k','linewidth',1)
legend('T','T_s','S')
set(h1,'color','k','linewidth',1)
set(h2,'linewidth',1)
xlim(Ax(1),[0 90])
xlim(Ax(2),[0 90])
title('ML S and T')
xlabel('Time (days)')

grid on
box on
set(gca,'ydir','normal','layer','top','fontname','helvetica','fontsize',12)
ylabel(Ax(1),'^\circ C')
set(Ax(1),'ycolor','k')
ylabel(Ax(2),'psu')

subplot(235)
plot(FSTD.time/86400,DIAG.OCEAN.Q_ml_SW(2:end),'linewidth',1)
hold on
plot(FSTD.time/86400,DIAG.OCEAN.Q_mi(2:end),'linewidth',1)
plot(FSTD.time/86400,DIAG.OCEAN.Qsurf_ml(2:end),'linewidth',1)
plot(FSTD.time/86400,-DIAG.OCEAN.Q_base_mix(2:end),'linewidth',1)
hold off
grid on
box on
set(gca,'ydir','normal','layer','top','fontname','helvetica','fontsize',12)
legend('SW heating','To Ice','To Surface','To Deep')
xlim([FSTD.time(1) FSTD.time(end)]/86400)
title('ML Heating')
ylabel('W/m^2')
xlabel('Time (days)')

subplot(236)
plot(FSTD.time/86400,DIAG.OCEAN.S_ml_out(2:end),'linewidth',1)
hold on
plot(FSTD.time/86400,DIAG.OCEAN.S_base_mix(2:end),'linewidth',1);
grid on
box on
set(gca,'ydir','normal','layer','top','fontname','helvetica','fontsize',12)
xlim([FSTD.time(1) FSTD.time(end)]/86400)
legend('From Ice','To Deep')
ylabel('psu m/s')
xlabel('Time (days)')
title('ML Salt Fluxes')
% for i = 1:6
%     
%     subplot(2,3,i)
%     grid on
%     box on
%     
%     if i ~=2 && i~=3
%     
%     xlim([FSTD.time(1) FSTD.time(end)]/86400)
%     
%     end
%     
% end
tightfig
pos = [10 6]
set(gcf,'windowstyle','normal','position',[0 0 pos],'paperposition',[0 0 pos],'papersize',pos,'units','inches','paperunits','inches');
set(gcf,'windowstyle','normal','position',[0 0 pos],'paperposition',[0 0 pos],'papersize',pos,'units','inches','paperunits','inches');
saveas(gcf,'~/Dropbox/FSTD/Manuscripts/FSTD-Response/Figures/Fig-4-Ocean-SS/Fig-4.pdf')
saveas(gcf,'~/Dropbox/FSTD/Manuscripts/FSTD-Response/Figures/Fig-4-Ocean-SS/Fig-4.fig')