% post_plot_data;
% 
% 
% if isfield(PLOTS,'fig_pplot')
%     
%     try
%         close(PLOTS.fig_pplot)
%     catch err
%     end
%     
% end

[~,b] = find(FSTD.Rint > 5,1);
[~,c] = find(FSTD.Rint >= 500,1);
[~,d] = find(DIAG.FSTD.conc(2:end) < .01,1);

PLOTS.fig_pplot = figure;
% this is an optional routine run at the end that plots what I want for
% making output figures


cplots = [228,26,28
    55,126,184
    77,175,74
    152,78,163
    255,127,0
    255,255,51
    166,86,40]/256;

cmap2 = [255,247,236
    254,232,200
    253,212,158
    253,187,132
    252,141,89
    239,101,72
    215,48,31
    179,0,0
    127,0,0]/256;




dayend = floor(FSTD.time(end)/86400);

dur = 30;

daybeg = 2;

xlimmer = [0 60];

%% Make the Contour plots of ITD/FSD
Ax{3} = subplot(133);

inds = 1:14:length(FSTD.time);
nplots = length(inds);
% Plot FSD as line plots
    lightplots = [103,0,31
        247,251,255
        222,235,247
        198,219,239
        158,202,225
        107,174,214
        66,146,198
        33,113,181
        8,81,156
        8,48,107]/256; 
         
    
    
    div = size(lightplots,1)/nplots;
    
    for i = 1:3
        colplots(:,i) = downsample(lightplots(:,i),floor(div));
    end
    
    
    
    colplots(colplots > 1) = 1;
    colplots(colplots < 0) = 0;



hold on

per = squeeze(sum(bsxfun(@times,DIAG.FSTD.psi(:,:,inds),DIAG.FSTD.dA(:,:,inds)),2)+eps);

per = bsxfun(@rdivide,per,DIAG.FSTD.conc(inds));

% per = log10(per);

for i = 1:length(inds)
    
    semilogy(DIAG.FSTD.R,per(:,i),'color',colplots(i,:),'linewidth',2);
    
    str{i} = ['Day: ' num2str(round(DIAG.FSTD.time(inds(i))/86400))];
    
end

%%
xlim([DIAG.FSTD.R(1) DIAG.FSTD.R(end)])


llim = floor(log10(1/length(DIAG.FSTD.R)) - 1);

grid on
box on
xlabel('Floe Size')
title('FSD(r) (normalized to 1)')
ylabel('log10(m^2/m^2)')
set(gca,'ydir','normal','layer','top','fontname','helvetica','fontsize',14)
% legend(str)
% set(gca,'ylim',[llim 0])

shading interp

set(gca,'xscale','log')

set(gca,'yscale','log')

axis xy
grid on
box on
ylabel('FSD')
title('log_{10} of FSD(r)')
xlabel('Floe Size (m)')
set(gca,'ydir','normal','layer','top','fontname','helvetica','fontsize',14)
ylim([1e-3 .1])
xlim([FSTD.Rmid(b) FSTD.Rmid(c)])
set(gca,'xtick',[FSTD.Rmid(b) 10 50 100 FSTD.Rmid(c) 1000],'xticklabel',{'5','10','50','10^{2}','500','10^{3}'}) 
shading interp
legend(str)
%% Plot ITD

% Ax{4} = subplot(224);
% 
% if length(FSTD.Hmid) > 1
%     
%     % determine interpolation grid: from min to max in equal steps
%     ax1 = [0 FSTD.time/86400];
%     ax1 = ax1(1:skip_contour:end);
%     
%     ax2 = FSTD.Hmid_i;
%     
%     Ax1v = linspace(ax1(1),ax1(end),numel(ax1));
%     Ax2v = linspace(ax2(1),ax2(end),numel(ax2));
%     
%     [ax1,ax2] = meshgrid(ax1,ax2);
%     [Ax1,Ax2] = meshgrid(Ax1v,Ax2v);
%     
%     
%     
%     plotter = log10(squeeze(sum(bsxfun(@times,DIAG.FSTD.psi(:,:,1:skip_contour:end),DIAG.FSTD.dA(:,:,1:skip_contour:end)),1)+eps));
%     
%     int = interp2(ax1,ax2,plotter,Ax1,Ax2);
%     
%     imagesc(Ax1v,Ax2v,interp2(int,2));
%     
%     
%     
%     hold on
%     imcontour(ax1,ax2,log10(squeeze(sum(bsxfun(@times,DIAG.FSTD.psi(:,:,1:skip_contour:end),DIAG.FSTD.dA(:,:,1:skip_contour:end)),1)+eps)),[-4:1:0],'--k','showtext','on');
%     hold off
%     
% end
% 
% grid on
% box on
% ylabel('Floe Size')
% title('log_{10} of ITD(h)dh')
% xlabel('Time')
% set(gca,'ydir','normal','layer','top','fontname','helvetica','fontsize',14)
% 
% if length(FSTD.Hmid > 1)
%     % ylim([FSTD.Hmid(1) FSTD.Hmid_i(end)]);
% end
% 
% set(gca,'clim',[-4 0])
% shading interp
% axis xy
% % colorbar
% 
% colormap(cmap2)

%% Plot the power-law behavior

[~,b] = find(FSTD.Rint > 5,1);
[~,c] = find(FSTD.Rint >= 500,1);
[~,d] = find(DIAG.FSTD.conc(2:end) < .01,1);

if isempty(d)
    d = length(DIAG.FSTD.conc);
end


Rinert = FSTD.Rint(b:c);
FSDinert = DIAG.FSTD.psi(b:c,:,1:d);
% dAinert = DIAG.FSTD.dA(b:c,:,1:d);
FSDinert = squeeze(sum(bsxfun(@times,FSDinert,FSTD.dH),2));
dRinert = FSTD.dR(b:c)';
% FSDinert = bsxfun(@rdivide,FSDinert, sum(FSDinert,1));
% mfs = sum(bsxfun(@times,Rinert',FSDinert));

FSDinertlog = log10(FSDinert + eps);

for i = 1:size(FSDinert,2)
    
    c_a(i) = sum(FSDinert(:,i).*dRinert,1);
    c_p(i) = 2*sum(FSDinert(:,i).*dRinert.*FSTD.Rint(b:c).^(-1)');
    p = polyfit(log10(Rinert),FSDinertlog(:,i)',1);
    coeff(i) = -p(1);
    
end

% Compare to the value of alpha suggested for the floe size distribution if
% we know it decays as f(r) = r^(-\alpha) from r1 to inf.

% This is the power-law decay of this function, beta
plaw = coeff(2:d);
% On the abcissa we plot the ice concentration
abciss = DIAG.FSTD.conc(2:d);
% Plot it!
% plot(abciss,beta,'linewidth',1,'color',cplots(1,:))
hold on



R = 1;
r1 = Rinert(1)-dRinert(1)/2;
D1 = 2*r1;

% Compare to the value of alpha suggested for the floe size distribution if
% we know it decays as f(r) = r^(-\alpha) from r1 to inf.

if min(abciss) == max(abciss)
    abciss = FSTD.time/86400;
end

horvat_alpha = (2* c_a ./ (2*c_a - r1*c_p));
beta = r1*c_p./(2*c_a);
beta_plaw = (plaw-1)./(plaw);

perovich_alpha = (8 * c_a * R - c_p * D1)./(4*c_a * R - c_p*D1) - 1;

% plot(abciss,perovich_alpha(2:d),'linewidth',1,'color',cplots(2,:))
% plot(abciss,horvat_alpha(2:d),'linewidth',1,'color',cplots(3,:))


%% Do goodness-of-fit-test

b = 1; 
c = size(DIAG.FSTD.psi,1); 
d = d;
VC_FSD = DIAG.FSTD.psi(b:c,:,1:d);
VC_FSD = squeeze(sum(bsxfun(@times,VC_FSD,FSTD.dH),2));
VC_dR = FSTD.dR(b:c)';
VC_R = FSTD.Rint(b:c);


[VC_alpha,VC_xmin,VC_L,UN_alpha,UN_xmin,UN_n,VC_p,VC_gof] = calc_power_law_stuff(VC_FSD,VC_dR,VC_R); 

%%
subplot(1,3,2)

Ax{2} = gca;

hold on
plot(c_a(2:d),plaw,'color',cplots(1,:))
plot(c_a(2:d),horvat_alpha(2:d),'--','color',cplots(1,:))

VC_y = smooth(VC_alpha(2:d)); 
VC_x = c_a(2:d); 
VC_err = UN_alpha(2:d); 

plot(c_a(2:d),VC_y,'-k')

patch([VC_x fliplr(VC_x)],[VC_y'+VC_err fliplr(VC_y'-VC_err)],[0.7 0.7 0.7],'facealpha',.5);
% plot(c_a(2:d),beta(2:d),'--','color',cplots(2,:))

%%
set(Ax{2},'xdir','reverse','ycolor','k','xcolor','k')
% set(Ax{5},'xdir','reverse','ycolor','r','xcolor','k')

grid on
box on
xlabel('Ice Concentration')
% ylabel('Power Law')
title('Power-Law Decay')
xlim([min(abciss) max(abciss)])
ylim([1 2.5])
legend('Fit \alpha','P+J \alpha')


%%


subplot(131);
hold off
[AX,H1,H2] = plotyy(FSTD.time/OPTS.day,DIAG.FSTD.Vtot(2:end),FSTD.time/OPTS.day,DIAG.FSTD.conc(2:end))

Ax{1} = AX(1); 
Ax{4} = AX(2);

% axis(Ax{1})
ylabel(Ax{4},'%')
% ,'color',cplots(1,:))
hold on
% plot(FSTD.time/OPTS.day,DIAG.FSTD.Vtot(2:end),'color',cplots(2,:))
axis(Ax{4})
plot(FSTD.time/OPTS.day,DIAG.FSTD.Hmean(2:end),'color',cplots(3,:))
hold off
ylabel('m')
xlabel('Time (days)')
ylim(Ax{4},[0 1])
ylim(Ax{1},[0 2])
title('Ice Variables')
legend('V','H','C','location','northwest')

% set(gca,'xtick',ticks,'xticklabel',mos)
%
% Ax{2} = subplot(232);
% rho_oc = OCEAN.EOS(DIAG.OCEAN.T,DIAG.OCEAN.S);
% rho_b = OCEAN.EOS(OCEAN.T_b(1),OCEAN.S_b(1));
%
% plot(FSTD.time/OPTS.day,rho_b - rho_oc(2:end),'color',cplots(1,:))
% ylabel('kg/m^3')
% % set(gca,'xtick',ticks,'xticklabel',mos)
% title('\Delta \rho')
%
% Ax{3} = subplot(233);
% [Ax1,h1,h2] = plotyy(FSTD.time/OPTS.day,DIAG.OCEAN.T(2:end),FSTD.time/OPTS.day,DIAG.OCEAN.S(2:end));
% hold on
% plot(FSTD.time/OPTS.day,DIAG.OCEAN.T_s(2:end),'--','color',cplots(1,:))
% hold off
%
% ylabel(Ax1(1),'^\circ C')
% ylabel(Ax1(2),'psu')
% set(Ax1(1),'ycolor',cplots(1,:))
% set(Ax1(2),'ycolor',cplots(2,:))
% set(h1,'color',cplots(1,:),'linewidth',1);
% set(h2,'color',cplots(2,:),'linewidth',2);
% title('ML T and S')
% legend('T','T_s','S','location','northwest')
% % set(Ax1(1),'xtick',ticks,'xticklabel',mos)
% % set(Ax1(2),'xtick',ticks,'xticklabel',mos)
%
% Ax{4} = subplot(234);
% plot(FSTD.time/OPTS.day,DIAG.OCEAN.Qsurf_at(2:end),'color','k','linewidth',1)
% hold on
% plot(FSTD.time/OPTS.day,DIAG.OCEAN.Qsurf_ml(2:end),'color',cplots(1,:),'linewidth',1)
% plot(FSTD.time/OPTS.day,DIAG.OCEAN.Qlead(2:end),'color',cplots(2,:),'linewidth',1)
% % plot(FSTD.time/OPTS.day,DIAG.OCEAN.QSW_in(2:end),'color',cplots(3,:),'linewidth',1)
% % plot(FSTD.time/OPTS.day,-DIAG.OCEAN.QSH(2:end),'color',cplots(4,:),'linewidth',1)
% % plot(FSTD.time/OPTS.day,-DIAG.OCEAN.QLH(2:end),'color',cplots(5,:),'linewidth',1)
% % plot(FSTD.time/OPTS.day,DIAG.OCEAN.QLW_in(2:end) - DIAG.OCEAN.QLW_out(2:end),'color',cplots(6,:),'linewidth',1)
% hold off
% title('Surface Fluxes')
% % legend('SW','SH','LH','Net LW','To ML','To Ice','net')
% legend('Surface Q','To ML','To Ice','location','northwest')
% % legend('orientation','horizontal')
% xlabel('Time (days)')
% ylabel('W/m^2')
% % set(gca,'xtick',ticks,'xticklabel',mos)
%
% Ax{5} = subplot(235);
% % Plot the heat fluxes in the mixed layer
% plot(FSTD.time/OPTS.day,DIAG.OCEAN.Qsurf_ml(2:end),'color',cplots(1,:),'linewidth',1)
% hold on
% plot(FSTD.time/OPTS.day,DIAG.OCEAN.Q_ml_SW(2:end),'color',cplots(2,:),'linewidth',1)
% plot(FSTD.time/OPTS.day,DIAG.OCEAN.Q_base_mix(2:end),'color',cplots(3,:),'linewidth',1)
% plot(FSTD.time/OPTS.day,OCEAN.rho * OCEAN.cw * DIAG.OCEAN.H_ml(2:end) .* DIAG.OCEAN.dTdt(2:end),'color',cplots(4,:),'linewidth',1)
% hold off
% title('ML T Fluxes')
% legend('To surface','SW','Mixing at base','Net','location','northwest')
% % legend('orientation','horizontal')
% xlabel('Time (days)')
% ylabel('W/m^2')
% % set(gca,'xtick',ticks,'xticklabel',mos)
%
% Ax{6} = subplot(236);
% % Plot the salinity fluxes
% plot(FSTD.time/OPTS.day,DIAG.OCEAN.S_base_mix(2:end),'color',cplots(1,:),'linewidth',1);
% hold on
% plot(FSTD.time/OPTS.day,DIAG.OCEAN.S_ml_ice(2:end),'color',cplots(2,:),'linewidth',1);
% plot(FSTD.time/OPTS.day,DIAG.OCEAN.S_ml_precip(2:end),'color',cplots(3,:),'linewidth',1);
% plot(FSTD.time/OPTS.day,DIAG.OCEAN.S_ml_evap(2:end),'color',cplots(4,:),'linewidth',1);
% hold off
% title('ML S Fluxes')
% xlabel('Time (days)')
% ylabel('psu m/s')
% legend('Mixing at base','Melting','Precip','Evap','location','northwest')
% % set(gca,'xtick',ticks,'xticklabel',mos)

%% Make all the axes to the right format
pos = [12 4];
set(gcf,'windowstyle','normal','position',[0 0 pos],'paperposition',[0 0 pos],'papersize',pos,'units','inches','paperunits','inches');
set(gcf,'windowstyle','normal','position',[0 0 pos],'paperposition',[0 0 pos],'papersize',pos,'units','inches','paperunits','inches');

drawnow


letter = {'(a)','(b)','(c)','(d)','(e)','(f)','(g)','(h)'};

hAllAxes = findobj(gcf,'type','axes');
hLeg = findobj(hAllAxes,'tag','legend');
hAxes = setdiff(hAllAxes,hLeg); % All axes which are not

delete(findall(gcf,'Tag','legtag'))

% tightfig


for i = 1:length(Ax)
    posy = get(Ax{i},'position');
    if i~=4
    annotation('textbox',[posy(1)-.04 posy(2)+posy(4)+.04 .025 .025], ...
        'String',letter{i},'LineStyle','none','FontName','Helvetica', ...
        'FontSize',14,'Tag','legtag');
    end
    grid(Ax{i},'on')
    box(Ax{i},'on')
    if i~= 3 && i~=2
        set(Ax{i},'xlim',xlimmer)
    end
    set(Ax{i},'ydir','normal','layer','top','fontname','helvetica','fontsize',12)
end


% saveas(gcf,'~/Dropbox/FSTD/Manuscripts/FSTD-Response/Figures/Fig-2/Fig-2.pdf')
% saveas(gcf,'~/Dropbox/FSTD/Manuscripts/FSTD-Response/Figures/Fig-2/Fig-2.fig')
