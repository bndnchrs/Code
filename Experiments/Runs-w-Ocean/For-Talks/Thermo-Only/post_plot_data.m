
load /Users/Horvat/Research/FSTD-Code/Output/For-Paper/Thermo-Only/Example.mat

close all
figure

skipper = 1; 

[~,b] = find(FSTD.Rint > 5,1);
[~,c] = find(FSTD.Rint >= 500,1);

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
% Ax{3} = subplot(133);

per = squeeze(sum(bsxfun(@times,DIAG.FSTD.psi(:,:,1:skipper:end),FSTD.dH),2)+eps);

per = bsxfun(@rdivide,per,DIAG.FSTD.conc(1:skipper:end));


loglog(DIAG.FSTD.R,per(:,1),'k','linewidth',2);
legend('Initial FSD')

xlim([6 1500])

hold on


axis xy
grid on
box on
ylabel('FSD')
title('FSD(r): Day 0')
xlabel('Floe Size (m)')
set(gca,'ydir','normal','layer','top','fontname','helvetica','fontsize',14)
ylim([1e-5 .2])
xlim([FSTD.Rmid(b) FSTD.Rmid(c)])
set(gca,'xtick',[FSTD.Rmid(b) 10 50 100 500 1000],'xticklabel',{'5','10','50','10^{2}','500','10^{3}'})
create_movie_FD('FSD',0,[0 0 6 4])

endt = 68; 

%%

timer = DIAG.FSTD.time(1:skipper:end);

for i = 1:endt-5
    
    
    pl = loglog(DIAG.FSTD.R,per(:,i),'color',[33,113,181]/256,'linewidth',2);
    
    
    xlim([FSTD.Rmid(b) FSTD.Rmid(c)])
    
    grid on
    box on
    xlabel('Floe Size')
    title('FSD(r)dr (normalized to 1)')
    ylabel('log10(m^2/m^2)')
    set(gca,'ydir','normal','layer','top','fontname','helvetica','fontsize',14)
    
    daystr = num2str(floor(timer(i)/86400));
    shading interp
    
    set(gca,'xscale','log')
    
    set(gca,'yscale','log')
    
    axis xy
    grid on
    box on
    ylabel('FSD')
    title(['Normalized FSD: Day ' daystr])
    xlabel('Floe Size (m)')
    set(gca,'ydir','normal','layer','top','fontname','helvetica','fontsize',14)
    ylim([1e-5 .2])
    xlim([FSTD.Rmid(b) 1500])
    set(gca,'xtick',[FSTD.Rmid(b) 10 50 100 500 1000],'xticklabel',{'5','10','50','10^{2}','500','10^{3}'})
    create_movie_FD('FSD',i,[0 0 6 4])
    
    
    delete(pl)
    
end

create_movie_FD('close',i,[0 0 8 3])

%% Plot the power-law behavior

% if isempty(d)
    d = length(DIAG.FSTD.conc);
%end


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
    [p,S] = polyfit(log10(Rinert),FSDinertlog(:,i)',1);
    [f,gof2] = fit(log10(Rinert)',FSDinertlog(:,i),'poly1');
    
    coeff(i) = -f.p1;
    errval(i) = gof2.rsquare;
    
end
    
% This is the power-law decay of this function, beta
plaw = coeff(2:d);
% On the abcissa we plot the ice concentration
abciss = DIAG.FSTD.conc(2:d);
% Plot it!
% plot(abciss,beta,'linewidth',1,'color',cplots(1,:))
hold on

figure
 

% plot(DIAG.FSTD.time/86400,errval,'--','color',cplots(1,:))


%%



R = 1;
r1 = Rinert(1)-dRinert(1)/2;
D1 = 2*r1;

% Compare to the value of alpha suggested for the floe size distribution if
% we know it decays as f(r) = r^(-\alpha) from r1 to inf.

%if min(abciss) == max(abciss)
    abciss = FSTD.time/86400;
%end

horvat_alpha = (2* c_a ./ (2*c_a - r1*c_p));
beta = r1*c_p./(2*c_a); 
beta_plaw = (1+plaw)./plaw; 

perovich_alpha = (8 * c_a * R - c_p * D1)./(4*c_a * R - c_p*D1) - 1;

%%
% plot(abciss,perovich_alpha(2:d),'linewidth',1,'color',cplots(2,:))
% plot(abciss,horvat_alpha(2:d),'linewidth',1,'color',cplots(3,:))

plot(DIAG.FSTD.time(1:end-1)/86400,plaw,'color',cplots(1,:))
hold on
% plot(c_a(2:d),beta_plaw','color',cplots(3,:))

plot(DIAG.FSTD.time/86400,horvat_alpha,'--','color',cplots(1,:))
% plot(c_a(2:d),beta(2:d),'--','color',cplots(2,:))

% set(gca,'xdir','reverse')
grid on
box on
xlabel('Time (days)')
ylabel('\alpha')
title('Power-Law Decay')

xlim([min(abciss) max(abciss)])
ylim([1 2])
legend('Computed','P+J')

%%

close all
figure

p = plaw;
alp = horvat_alpha(1:skipper:end); 

 plot(DIAG.FSTD.time(2:end)/86400,plaw,'color',cplots(1,:))
 hold on
 plot(DIAG.FSTD.time/86400,horvat_alpha,'color',cplots(2,:))

axis xy
grid on
box on
ylabel('\alpha')
title(['Power Laws'])
xlabel('Time (days)')
set(gca,'ydir','normal','layer','top','fontname','helvetica','fontsize',14)
ylim([1 2])
xlim([0 timer(endt-5)]/86400)

pos = [6 4]; 
set(gcf,'windowstyle','normal','position',[0 0 pos],'paperposition',[0 0 pos],'papersize',pos,'units','inches','paperunits','inches');
set(gcf,'windowstyle','normal','position',[0 0 pos],'paperposition',[0 0 pos],'papersize',pos,'units','inches','paperunits','inches');
legend('Fit','P+J')

saveas(gcf,'plaw_plot.pdf')

%%
create_movie_FD('PLAW',0,[0 0 4 3])

legend('Fit','P+J')



for i = 2:endt-5
    
    
    pl = scatter(timer(i)/86400,plaw(i-1),200,cplots(1,:),'filled');
    pl2 = scatter(timer(i)/86400,horvat_alpha(i),100,cplots(2,:),'filled');
    
    
    grid on
    box on
    set(gca,'ydir','normal','layer','top','fontname','helvetica','fontsize',14)
    
    daystr = num2str(floor(timer(i)/86400));
    shading interp
    
    axis xy
    grid on
    box on
    title(['Power Laws: Day ' daystr])
    create_movie_FD('PLAW',i,[0 0 4 3])
    
    
    delete(pl)
    delete(pl2)
    
end

create_movie_FD('close',i,[0 0 4 3])

title('Power-Law Decay')
xlim([0 DIAG.FSTD.time(endt-5)]/86400)
ylim([1 2.5])
%%

%% Make the Contour plots of ITD/FSD
% Ax{3} = subplot(133);

close all

per = squeeze(sum(bsxfun(@times,DIAG.FSTD.psi(:,:,:),FSTD.dH),2)+eps);

per = bsxfun(@rdivide,per,DIAG.FSTD.conc);
per = bsxfun(@rdivide,per,FSTD.Rint'.^2);
per = sum(bsxfun(@times,per,FSTD.dR'),1) - cumsum(bsxfun(@times,per,FSTD.dR'),1);


% per = bsxfun(@rdivide,per,FSTD.dR');

loglog(DIAG.FSTD.R,1*per(:,1),'k','linewidth',2);
legend('Initial CND')

xlim([6 1500])

hold on

%% 
axis xy
grid on
box on
ylabel('CND')
title('CND(r): Day 0')
xlabel('Floe Size (m)')
set(gca,'ydir','normal','layer','top','fontname','helvetica','fontsize',14)
ylim([1e-8 1e-2])
set(gca,'xtick',[FSTD.Rmid(b) 10 50 100 500 1000],'xticklabel',{'5','10','50','10^{2}','500','10^{3}'})
create_movie_FD('CND',0,[0 0 8 3])

%%

for i = 1:endt-5
    
    
    pl = loglog(DIAG.FSTD.R,per(:,i),'color',[33,113,181]/256,'linewidth',2);
    
    
    grid on
    box on
    xlabel('Floe Size')
    ylabel('log10(m^2/m^2)')
    set(gca,'ydir','normal','layer','top','fontname','helvetica','fontsize',14)
    
    daystr = num2str(floor(timer(i)/86400));
    shading interp
    
    set(gca,'xscale','log')
    
    set(gca,'yscale','log')
    
    axis xy
    grid on
    box on
    ylabel('CND')
    title(['CND(r): Day ' daystr])
    xlabel('Floe Size (m)')
    set(gca,'ydir','normal','layer','top','fontname','helvetica','fontsize',14)
    ylim([1e-8 1e-2])
    set(gca,'xtick',[FSTD.Rmid(b) 10 50 100 500 1000],'xticklabel',{'5','10','50','10^{2}','500','10^{3}'})
    create_movie_FD('CND',i,[0 0 8 3])
    
    
    delete(pl)
    
end

create_movie_FD('close',i,[0 0 8 3])
