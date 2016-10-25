% post_plot_data;
if isfield(PLOTS,'fig_pplot')
    
    try
        close(PLOTS.fig_pplot)
    catch err
    end
    
end

[~,b] = find(FSTD.Rint > 5,1);
[~,c] = find(FSTD.Rint > 50,1);
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

xlimmer = [0 14];

%% Make the Contour plots of ITD/FSD
Ax{3} = subplot(133);

inds = 1:7*23:length(FSTD.time);
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

clear colplots

for i = 1:3
    colplots(:,i) = downsample(lightplots(:,i),floor(div));
end



colplots(colplots > 1) = 1;
colplots(colplots < 0) = 0;



hold on

per = squeeze(sum(bsxfun(@times,DIAG.FSTD.psi(:,:,inds),FSTD.dH),2)+eps);

per = bsxfun(@rdivide,per,DIAG.FSTD.conc(inds));

% per = log10(per);
str = cell(1);

for i = 1:length(inds)
    
    semilogy(DIAG.FSTD.R,per(:,i),'color',colplots(i,:),'linewidth',2);
    
    str{i} = ['Day: ' num2str(ceil(DIAG.FSTD.time(inds(i))/86400))];
    
end

plotter = squeeze(sum(bsxfun(@times,ADVECT.FSTD_in,FSTD.dH),2)) + eps;


plotter = plotter ./ integrate_FSTD(ADVECT.FSTD_in,1,FSTD.dA,0);

semilogy(DIAG.FSTD.R,plotter,'--','color',cplots(3,:),'linewidth',2)

str{end+1} = 'Pack Ice';
%%
xlim([DIAG.FSTD.R(1) DIAG.FSTD.R(end)])


llim = floor(log10(1/length(DIAG.FSTD.R)) - 1);

grid on
box on
xlabel('Floe Size')
title('FSD(r)dr (normalized to 1)')
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
ylabel('Floe Size')
title('log_{10} of FSD(r)dr')
xlabel('Time')
set(gca,'ydir','normal','layer','top','fontname','helvetica','fontsize',14)
ylim([1e-4 1])
xlim([FSTD.Rmid(b) 500])
set(gca,'xtick',[FSTD.Rmid(b) 10 50 100 500 1000],'xticklabel',{'5','10','50','10^{2}','500','10^{3}'})
shading interp
legend(str)

%% Plot the power-law behavior

Ax{2} = subplot(132);

cvals = [50 500];

for ind = 1
    
    [~,c] = find(FSTD.Rint > cvals(ind),1);
    
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
    
    % This is the power-law decay of this function, beta
    plaw = coeff(1:d);
    % On the abcissa we plot the ice concentration
    abciss = c_a;
    % Plot it!
    % plot(abciss,beta,'linewidth',1,'color',cplots(1,:))
    hold on
    
    
    
    R = 1;
    r1 = Rinert(1)-dRinert(1)/2;
    D1 = 2*r1;
    
    % Compare to the value of alpha suggested for the floe size distribution if
    % we know it decays as f(r) = r^(-\alpha) from r1 to inf.
    
    % if abs(min(abciss) - max(abciss)) <= .2
    abciss = [0 FSTD.time/86400];
    % abciss = DIAG.FSTD.conc;
    % end
    
    horvat_alpha = (2* c_a ./ (2*c_a - r1*c_p));
    beta = r1*c_p./(2*c_a);
    beta_plaw = (1+plaw)./plaw;
    
    perovich_alpha = (8 * c_a * R - c_p * D1)./(4*c_a * R - c_p*D1) - 1;
    % plot(abciss,perovich_alpha(2:d),'linewidth',1,'color',cplots(2,:))
    % plot(abciss,horvat_alpha(2:d),'linewidth',1,'color',cplots(3,:))
    
    if ind ==1
        plot(abciss,plaw,'color','k')
        hold on
        % plot(c_a(2:d),beta_plaw','color',cplots(3,:))
        plot(abciss,horvat_alpha,'--','color','k')
        % plot(c_a(2:d),beta(2:d),'--','color',cplots(2,:))
        good_plaw = plaw; 
    else
        plot(abciss,horvat_alpha,'--','color',cplots(3,:))
    end
    
    
    
    set(gca,'xdir','normal')
    grid on
    box on
    xlabel('Time (days)')
    ylabel('\alpha')
    title('Power-Law Decay')
    
    % xlim([min(abciss) max(abciss)])
    ylim([1.5 2.5])

    
end

legend('Computed - 50 m','P+J - 50 m','P+J - 500 m','location','southwest')

%%


Ax{1} = subplot(131);

plot(FSTD.time/OPTS.day,DIAG.FSTD.Rmeanarea(2:end),'color',cplots(1,:))
hold on
%plot(FSTD.time/OPTS.day,DIAG.FSTD.Vtot(2:end),'color',cplots(2,:))
% plot(FSTD.time/OPTS.day,DIAG.FSTD.Hmean(2:end),'color',cplots(3,:))
hold off

ylabel('m')
title('Mean Floe Size')
% legend('Mean Floe Size','location','northwest')

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
pos = [12 3];
set(gcf,'windowstyle','normal','position',[0 0 pos],'paperposition',[0 0 pos],'papersize',pos,'units','inches','paperunits','inches');
set(gcf,'windowstyle','normal','position',[0 0 pos],'paperposition',[0 0 pos],'papersize',pos,'units','inches','paperunits','inches');

drawnow


letter = {'(a)','(b)','(c)','(d)','(e)','(f)','(g)','(h)'};

hAllAxes = findobj(gcf,'type','axes');
hLeg = findobj(hAllAxes,'tag','legend');
hAxes = setdiff(hAllAxes,hLeg); % All axes which are not

delete(findall(gcf,'Tag','legtag'))

tightfig


for i = 1:length(Ax)
    posy = get(Ax{i},'position');
    annotation('textbox',[posy(1)-.025 posy(2)+posy(4)+.1 .025 .025], ...
        'String',letter{i},'LineStyle','none','FontName','Helvetica', ...
        'FontSize',14,'Tag','legtag');
    grid(Ax{i},'on')
    box(Ax{i},'on')
    if i~= 2 && i~=3
        set(Ax{i},'xlim',xlimmer)
    end
    set(Ax{i},'ydir','normal','layer','top','fontname','helvetica','fontsize',12)
end


saveas(gcf,'~/Dropbox/FSTD/Manuscripts/FSTD-Response/Figures/Model-Output/For-Paper/Waves-Advect/Fig-6.pdf')
saveas(gcf,'~/Dropbox/FSTD/Manuscripts/FSTD-Response/Figures/Model-Output/For-Paper/Waves-Advect/Fig-6.fig')
