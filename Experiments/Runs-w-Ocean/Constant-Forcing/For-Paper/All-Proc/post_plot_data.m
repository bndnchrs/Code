% post_plot_data;
% if isfield(PLOTS,'fig_pplot')
%     
%     try
%         close(PLOTS.fig_pplot)
%     catch err
%     end
%     
% end

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

xlimmer = [0 90];

%% Make the Contour plots of ITD/FSD
Ax{2} = subplot(232);

inds = [1 6 6*7 6*60];
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

xlim([DIAG.FSTD.R(1) DIAG.FSTD.R(end)])


llim = floor(log10(1/length(DIAG.FSTD.R)) - 1);

grid on
box on
xlabel('Floe Size')
title('FSD(r) (normalized to 1)')
% ylabel('log10(m^2/m^2)')
set(gca,'ydir','normal','layer','top','fontname','helvetica','fontsize',14)
% legend(str)
% set(gca,'ylim',[llim 0])

shading interp

set(gca,'xscale','log')

set(gca,'yscale','log')

axis xy
grid on
box on
% ylabel('Floe Size')
title('log_{10} of FSD(r)')
xlabel('')
set(gca,'ydir','normal','layer','top','fontname','helvetica','fontsize',14)
ylim([1e-7 .2])
xlim([FSTD.Rmid(b) 1500])
set(gca,'xtick',[FSTD.Rmid(b) 10 50 100 500 1000],'xticklabel',{'5','10','50','10^{2}','500','10^{3}'})
shading interp
legend(str)


xvals = 50 + 0*(1:100); 
ylimmer = get(gca,'ylim');
yvals = linspace(ylimmer(1),ylimmer(2),100); 
plot(xvals,yvals,'--k')

xvals = 150 + 0*(1:100); 
ylimmer = get(gca,'ylim');
yvals = linspace(ylimmer(1),ylimmer(2),100); 
plot(xvals,yvals,'--k')


%% Plot the power-law behavior

Ax{4} = subplot(234);

intervalx = [5 50 150];
intervalend = [50 150 1500];

for ind = 1
    
    
    [~,b] = find(FSTD.Rint > intervalx(ind),1);
    [~,c] = find(FSTD.Rint > intervalend(ind),1);
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
   
    
    FSDinertlog = log10(FSDinert + eps);
    
    Ninert = bsxfun(@rdivide,FSDinert,(pi * (Rinert').^2)); 
    Ntot = sum(bsxfun(@times,Ninert,dRinert),1); 
    
    for i = 1:size(FSDinert,2)
        
        c_a(i) = sum(FSDinert(:,i).*dRinert,1);
        c_p(i) = 2*sum(FSDinert(:,i).*dRinert.*FSTD.Rint(b:c).^(-1)');
        [p,S] = polyfit(log10(Rinert),FSDinertlog(:,i)',1);
    coeff(i) = -p(1);
    coeff2(i) = -p(2); 
    normr{ind}(i) = S.normr; 
    cc{ind}(i) = coeff(i); 
    
    FSD2 = Rinert .^(-coeff(i));
    N2 = FSD2 ./ (pi * (Rinert).^2);
    Ntot2 = sum(N2.*dRinert'); 
  %  FSD2 = FSD2 / sum(FSD2); 
    
    mfs{ind,1}(i) = sum(Ninert(:,i).*dRinert'.*2*pi*Rinert') / Ntot(i);
    mfs{ind,2}(i) = sum(2*FSD2.*dRinert./Rinert) / sum(FSD2.*dRinert ./ (pi*Rinert.^2));  
    
  %   sum(bsxfun(@times,2./Rinert',FSD2(:,i))) ./ sum(bsxfun(@times,1./(pi*Rinert.^2)',FSDinert));

    
    
    
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
    
    plot(abciss,plaw,'color',cplots(ind,:))
    hold on
    % plot(c_a(2:d),beta_plaw','color',cplots(3,:))
    
    plot(abciss,horvat_alpha,'--','color',cplots(ind,:))
    % plot(c_a(2:d),beta(2:d),'--','color',cplots(2,:))
    
    legstr{2*(ind-1) + 1} = ['Computed ' num2str(intervalx(ind)) '-' num2str(intervalend(ind))]; 
    legstr{2*(ind-1) + 2} = ['P+J ' num2str(intervalx(ind)) '-' num2str(intervalend(ind))]; 

    
    plaw(1)
end

set(gca,'xdir','normal')
grid on
box on
xlabel('Time (days)')
ylabel('\alpha')
title('Power-Law Decay')

xlim([0 90])
% ylim([1.5 2.5])
legend({'5-50','P+J','50-150','P+J','150-1500','P+J'})
%legend(legstr)
%%


Ax{1} = subplot(231);

[AX,H1,H2] = plotyy(FSTD.time/OPTS.day,DIAG.FSTD.Vtot(2:end),FSTD.time/OPTS.day,DIAG.FSTD.conc(2:end));
Ax{1} = AX(1); 
Ax{6} = AX(2);
ylim(Ax{1},[0 2])
ylim(Ax{6},[0 1])
set(Ax{1},'ytick',[0 .5 1 1.5 2])
set(Ax{6},'ytick',[0 .25 .5 .75 1])
% axis(Ax{1})
ylabel(Ax{6},'%')
% ,'color',cplots(1,:))
hold on
% plot(FSTD.time/OPTS.day,DIAG.FSTD.Vtot(2:end),'color',cplots(2,:))
axis(Ax{1})
plot(FSTD.time/OPTS.day,DIAG.FSTD.Hmean(2:end),'color',cplots(3,:))

hold off
ylabel('m')
xlabel('Time (days)')
title('Ice Variables')
legend('V','H','C','location','northwest')
set(Ax{1},'xlim',xlimmer)
set(Ax{6},'xlim',xlimmer)

%% Plot ITD

Ax{3} = subplot(233);

hold on

per = squeeze(sum(bsxfun(@times,DIAG.FSTD.psi(:,:,inds),FSTD.dR'),1)+eps);

per = bsxfun(@rdivide,per,DIAG.FSTD.conc(inds));

% per = log10(per);
str = cell(1);

for i = 1:length(inds)
    
    semilogy(DIAG.FSTD.H,per(:,i),'color',colplots(i,:),'linewidth',2);
    
    str{i} = ['Day: ' num2str(ceil(DIAG.FSTD.time(inds(i))/86400))];
    
end

plotter = squeeze(sum(bsxfun(@times,ADVECT.FSTD_in,FSTD.dR'),1)) + eps;


plotter = plotter ./ integrate_FSTD(ADVECT.FSTD_in,1,FSTD.dA,0);

semilogy(DIAG.FSTD.H,plotter,'--','color',colplots(3,:),'linewidth',2)

str{end+1} = 'Pack Ice';

xlim([DIAG.FSTD.H(1) DIAG.FSTD.H(end)])


llim = floor(log10(1/length(DIAG.FSTD.H)) - 1);

grid on
box on
% xlabel('Floe Size')
% title('FSD(r) (normalized to 1)')
% ylabel('log10 of FSD(r)')
set(gca,'ydir','normal','layer','top','fontname','helvetica','fontsize',14)
% legend(str)
% set(gca,'ylim',[llim 0])

shading interp
set(gca,'yscale','log')

axis xy
grid on
box on
% ylabel('Ice Thickness')
title('log_{10} of ITD(h)')
xlabel('')
set(gca,'ydir','normal','layer','top','fontname','helvetica','fontsize',14)
ylim([1e-3 2])
xlim([FSTD.H(1) FSTD.H(end)])
%set(gca,'xtick',[FSTD.Rmid(b) 10 50 100 500 1000],'xticklabel',{'5','10','50','10^{2}','500','10^{3}'})
shading interp
legend(str)

%% Plot the balances of the ITD
Ax{6} = subplot(236);
dur = 6*7; 


dFnet = DIAG.FSTD.diff_ITD; 
dF{1} = DIAG.ADVECT.diff_ITD; 
dF{2} = DIAG.THERMO.diff_ITD; 
dF{3} = DIAG.MECH.diff_ITD; 
dF{4} = DIAG.WAVES.diff_ITD;

hold on

for i = 1:4
    plotter = mean(dF{i}(:,end-dur:end),2);
    
    semilogx(FSTD.H,plotter,'color',cplots(i,:)); 
    
end

shading interp

axis xy
grid on
box on
ylabel('1/s')
title('dITD/dt')
xlabel('Ice Thickness')
set(gca,'ydir','normal','layer','top','fontname','helvetica','fontsize',14)
shading interp
xlim([FSTD.H(1) FSTD.H(end)])
legend('Advection','Thermo','Mech','Waves','location','southeast')



%% Plot the balances of the FSD
Ax{5} = subplot(235);
dur = 6*7; 

%%
dFnet = DIAG.FSTD.diff_FSD; 
dF{1} = DIAG.ADVECT.diff_FSD; 
dF{2} = DIAG.THERMO.diff_FSD; 
dF{3} = DIAG.MECH.diff_FSD; 
dF{4} = DIAG.WAVES.diff_FSD;

hold on

for i = 1:4
    plotter = mean(dF{i}(:,end-dur:end),2);
    
    semilogx(FSTD.Rint,plotter,'color',cplots(i,:),'linewidth',2); 
    
end

shading interp

[~,b] = find(FSTD.Rint > 5,1);
[~,c] = find(FSTD.Rint > 50,1);
[~,d] = find(DIAG.FSTD.conc(2:end) < .01,1);

set(gca,'xscale','log')
axis xy
grid on
box on
ylabel('1/s')
title('dFSD/dt')
xlabel('Floe Size')
set(gca,'ydir','normal','layer','top','fontname','helvetica','fontsize',14)
set(gca,'xtick',[FSTD.Rmid(b) 10 50 100 500 1000],'xticklabel',{'5','10','50','10^{2}','500','10^{3}'})
shading interp
xlim([FSTD.Rmid(b) 1500])

legend('Advection','Thermo','Mech','Waves')


xvals = 50 + 0*(1:100); 
ylimmer = get(gca,'ylim');
yvals = linspace(ylimmer(1),ylimmer(2),100); 
plot(xvals,yvals,'--k')

% annotation('textbox',[.435 .65 .025 .025], ...
%     'String','I','LineStyle','none','FontName','Times', ...
%     'FontSize',20);
% annotation('textbox',[.5 .65 .025 .025], ...
%     'String','II','LineStyle','none','FontName','Times', ...
%     'FontSize',20);
% annotation('textbox',[.59 .65 .025 .025], ...
%     'String','III','LineStyle','none','FontName','Times', ...
%     'FontSize',20);
% annotation('textbox',[.435 .15 .025 .025], ...
%     'String','I','LineStyle','none','FontName','Times', ...
%     'FontSize',20);
% annotation('textbox',[.5 .15 .025 .025], ...
%     'String','II','LineStyle','none','FontName','Times', ...
%     'FontSize',20);
% annotation('textbox',[.59 .15 .025 .025], ...
%     'String','III','LineStyle','none','FontName','Times', ...
%     'FontSize',20);

xvals = 150 + 0*(1:100);
ylimmer = get(gca,'ylim');
yvals = linspace(ylimmer(1),ylimmer(2),100);
plot(xvals,yvals,'--k')


%% Make all the axes to the right format
pos = [12 6];
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
    annotation('textbox',[posy(1)-.04 posy(2)+posy(4) .025 .025], ...
        'String',letter{i},'LineStyle','none','FontName','Helvetica', ...
        'FontSize',14,'Tag','legtag');
    grid(Ax{i},'on')
    box(Ax{i},'on')
    set(Ax{i},'ydir','normal','layer','top','fontname','helvetica','fontsize',12)
end


% saveas(gcf,'~/Dropbox/FSTD/Manuscripts/FSTD-Response/Figures/Fig-5/Fig-5.pdf')
% saveas(gcf,'~/Dropbox/FSTD/Manuscripts/FSTD-Response/Figures/Fig-5/Fig-5.fig')
