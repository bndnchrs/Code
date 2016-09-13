% post_plot_data; 

cplots = [228,26,28
55,126,184
77,175,74
152,78,163
255,127,0
255,255,51
166,86,40]/256;


endyr = ceil(max(FSTD.time/OPTS.year)); 
ticks = (1:(4*endyr))/4  + 1/24 - .25; 
mos = {'Jan','Apr','Jul','Oct'}; 
mos = repmat(mos,[1 endyr]);

if isfield(PLOTS,'fig_pplot')
    
    try
    close(PLOTS.fig_pplot)
    catch err
    end
    
end

maxyr = ceil(FSTD.time(end)/OPTS.year); 
minyr = maxyr - 2; 

    xlimmer = [max(0,minyr) FSTD.time(end)/OPTS.year]; 



PLOTS.fig_pplot = figure;
% this is an optional routine run at the end that plots what I want for
% making output figures

Ax{1} = subplot(231);

plot(FSTD.time/OPTS.year,DIAG.FSTD.conc(2:end),'color',cplots(1,:))
hold on
plot(FSTD.time/OPTS.year,DIAG.FSTD.Vtot(2:end),'color',cplots(2,:))
plot(FSTD.time/OPTS.year,DIAG.FSTD.Hmean(2:end),'color',cplots(3,:))
hold off

ylabel('%,m')
title('Ice Variables')
legend('c','V','H','location','northwest')
set(gca,'xtick',ticks,'xticklabel',mos)

Ax{2} = subplot(232)
rho_oc = OCEAN.EOS(DIAG.OCEAN.T,DIAG.OCEAN.S); 
rho_b = OCEAN.EOS(EXFORC.T_b,OCEAN.S_b(1)); 

plot(FSTD.time/OPTS.year,rho_b - rho_oc(2:end),'color',cplots(1,:))
ylabel('kg/m^3')
set(gca,'xtick',ticks,'xticklabel',mos)
title('\Delta \rho')

Ax{3} = subplot(233);
[Ax1,h1,h2] = plotyy(FSTD.time/OPTS.year,DIAG.OCEAN.T(2:end),FSTD.time/OPTS.year,DIAG.OCEAN.S(2:end));
hold on
plot(FSTD.time/OPTS.year,DIAG.OCEAN.T_s(2:end),'--','color',cplots(1,:))
hold off

ylabel(Ax1(1),'^\circ C')
ylabel(Ax1(2),'psu')
set(Ax1(1),'ycolor',cplots(1,:))
set(Ax1(2),'ycolor',cplots(2,:))
set(h1,'color',cplots(1,:),'linewidth',1); 
set(h2,'color',cplots(2,:),'linewidth',2); 
title('ML T and S')
legend('T','T_s','S','location','northwest')
set(Ax1(1),'xtick',ticks,'xticklabel',mos)
set(Ax1(2),'xtick',ticks,'xticklabel',mos)

Ax{4} = subplot(234);
plot(FSTD.time/OPTS.year,DIAG.OCEAN.Qsurf_at(2:end),'color','k','linewidth',1)
hold on
plot(FSTD.time/OPTS.year,DIAG.OCEAN.Qsurf_ml(2:end),'color',cplots(1,:),'linewidth',1)
plot(FSTD.time/OPTS.year,DIAG.OCEAN.Qlead(2:end),'color',cplots(2,:),'linewidth',1)
% plot(FSTD.time/OPTS.year,DIAG.OCEAN.QSW_in(2:end),'color',cplots(3,:),'linewidth',1)
% plot(FSTD.time/OPTS.year,-DIAG.OCEAN.QSH(2:end),'color',cplots(4,:),'linewidth',1)
% plot(FSTD.time/OPTS.year,-DIAG.OCEAN.QLH(2:end),'color',cplots(5,:),'linewidth',1)
% plot(FSTD.time/OPTS.year,DIAG.OCEAN.QLW_in(2:end) - DIAG.OCEAN.QLW_out(2:end),'color',cplots(6,:),'linewidth',1)
hold off
title('Surface Fluxes')
% legend('SW','SH','LH','Net LW','To ML','To Ice','net')
legend('Surface Q','To ML','To Ice','location','northwest')
% legend('orientation','horizontal')
xlabel('Time (years)')
ylabel('W/m^2')
set(gca,'xtick',ticks,'xticklabel',mos)

Ax{5} = subplot(235);
% Plot the heat fluxes in the mixed layer
plot(FSTD.time/OPTS.year,DIAG.OCEAN.Qsurf_ml(2:end),'color',cplots(1,:),'linewidth',1)
hold on
plot(FSTD.time/OPTS.year,DIAG.OCEAN.Q_ml_SW(2:end),'color',cplots(2,:),'linewidth',1)
plot(FSTD.time/OPTS.year,DIAG.OCEAN.Q_base_mix(2:end),'color',cplots(3,:),'linewidth',1)
plot(FSTD.time/OPTS.year,OCEAN.rho * OCEAN.cw * DIAG.OCEAN.H_ml(2:end) .* DIAG.OCEAN.dTdt(2:end),'color',cplots(4,:),'linewidth',1)
hold off
title('ML T Fluxes')
legend('To surface','SW','Mixing at base','Net','location','northwest')
% legend('orientation','horizontal')
xlabel('Time (years)')
ylabel('W/m^2')
set(gca,'xtick',ticks,'xticklabel',mos)

Ax{6} = subplot(236);
% Plot the salinity fluxes
plot(FSTD.time/OPTS.year,DIAG.OCEAN.S_base_mix(2:end),'color',cplots(1,:),'linewidth',1);
hold on
plot(FSTD.time/OPTS.year,DIAG.OCEAN.S_ml_ice(2:end),'color',cplots(2,:),'linewidth',1);
plot(FSTD.time/OPTS.year,DIAG.OCEAN.S_ml_precip(2:end),'color',cplots(3,:),'linewidth',1);
plot(FSTD.time/OPTS.year,DIAG.OCEAN.S_ml_evap(2:end),'color',cplots(4,:),'linewidth',1);
hold off
title('ML S Fluxes')
xlabel('Time (years)')
ylabel('psu m/s')
legend('Mixing at base','Melting','Precip','Evap','location','northwest')
set(gca,'xtick',ticks,'xticklabel',mos)

%% Make all the axes to the right format
pos = [12 4];
set(gcf,'windowstyle','normal','position',[0 0 pos],'paperposition',[0 0 pos],'papersize',pos,'units','inches','paperunits','inches');
set(gcf,'windowstyle','normal','position',[0 0 pos],'paperposition',[0 0 pos],'papersize',pos,'units','inches','paperunits','inches');

drawnow

tightfig

letter = {'(a)','(b)','(c)','(d)','(e)','(f)','(g)','(h)'}; 

hAllAxes = findobj(gcf,'type','axes');
hLeg = findobj(hAllAxes,'tag','legend');
hAxes = setdiff(hAllAxes,hLeg); % All axes which are not

delete(findall(gcf,'Tag','legtag'))

for i = 1:length(Ax)
    posy = get(Ax{i},'position');
    annotation('textbox',[posy(1)-.04 posy(2)+posy(4)+.04 .025 .025], ...
        'String',letter{i},'LineStyle','none','FontName','Helvetica', ...
        'FontSize',14,'Tag','legtag');
    grid(Ax{i},'on')
    box(Ax{i},'on')
    set(Ax{i},'xlim',xlimmer)
    set(Ax{i},'ydir','normal','layer','top','fontname','helvetica','fontsize',12)
end

 grid(Ax1(2),'on')
    box(Ax1(2),'on')
    set(Ax1(2),'xlim',xlimmer)
    set(Ax1(2),'ydir','normal','layer','top','fontname','helvetica','fontsize',12)

    saveas(gcf,'~/Dropbox/FSTD/Manuscripts/FSTD-Response/Figures/Fig-6-seas-thermo/Fig-6.pdf')
    saveas(gcf,'~/Dropbox/FSTD/Manuscripts/FSTD-Response/Figures/Fig-6-seas-thermo/Fig-6.fig')
