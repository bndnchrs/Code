
close

load /Users/Horvat/Research/FSTD-Code/Output/For-Paper/Sens/one.mat

load plaws

plaws = [0 .5 1 1.5 2 2.5 3]; 


subplot(224)
loglog(FSTD.Rmid,saveper{1}(:,4))
xlabel('Floe Size (m)')
set(gca,'ytick',[1e-6 1e-4 1e-2]); 
set(gca,'xtick',[1 10 100 1000])
hold on
for i = 2:7
    loglog(FSTD.Rmid,saveper{i}(:,4))
end
legend('0','0.5','1','1.5','2','2.5','3','orientation','horizontal')
set(gco,'position',[.25    0.94    0.5    0.05])
ylim([1e-7 1])
xlim([6 1500])

subplot(223)
hold on
xlabel('Initial Power Law')
ylabel('Decay Coeff')
xlim([0 3])
ylim([-1 10])
for ind = 1:3
    for i = 1:7
        plawwer(i,ind) = mean(save_plaw{i,ind}(end-14:end));
    end
end
set(gca,'ytick',[-1 0 1 2 3 5 7 9])

m = {'d','s','o'}


intervalx = [5 50 200];
intervalend = [20 100 1500];

[~,b] = find(FSTD.Rint > intervalx(ind),1);
[~,c] = find(FSTD.Rint > intervalend(ind),1);
[~,d] = find(DIAG.FSTD.conc(2:end) < .01,1);

for ind = 1:3
    for i = 1:7
        
        
         load([OPTS.savepath OPTS.names{i}],'DIAG')
         
         
        [~,b] = find(FSTD.Rint > intervalx(ind),1);
        [~,c] = find(FSTD.Rint > intervalend(ind),1);
        [~,d] = find(DIAG.FSTD.conc(2:end) < .01,1);
        
        Rinert = FSTD.Rint(b:c);
        FSDinert = DIAG.FSTD.psi(b:c,:,end-14:end);
        FSDinert = squeeze(sum(bsxfun(@times,FSDinert,FSTD.dH),2));
        FSDinertlog = log10(mean(FSDinert,3)); 
  
        [p,S] = polyfit(log10(Rinert),FSDinertlog(:,i)',1);
        
        plawwer(i,ind) = -p(1); 
        
    end
    
end
        
    %%
    

for ind = 1:3
plot(plaws,plawwer(:,ind),'marker',m{ind},'linestyle','none'); 
end

legend('I','II','III')

subplot(221)
hold on
xlabel('Time (days)')
ylabel('Concentration')
ylim([0 1])
xlim([0 90])

subplot(222)
hold on
ylim([0 1])
xlim([0 90])
xlabel('Time (days)')
ylabel('Volume (m^3/m^2)')

for i = 1:7
    
    load([OPTS.savepath OPTS.names{i}],'DIAG')
    
    subplot(221)
    
    plot([0 FSTD.time]/86400,DIAG.FSTD.conc); 
    
    subplot(222)
        plot([0 FSTD.time]/86400,DIAG.FSTD.Vtot); 

    
end


for i = 1:4
    subplot(2,2,i)
    grid on
    box on
    set(gca,'ydir','normal','layer','top','fontname','helvetica','fontsize',18)
    
end

%%
pos = [12 8]; 

set(gcf,'windowstyle','normal','position',[0 0 pos],'paperposition',[0 0 pos],'papersize',pos,'units','inches','paperunits','inches');
set(gcf,'windowstyle','normal','position',[0 0 pos],'paperposition',[0 0 pos],'papersize',pos,'units','inches','paperunits','inches');

saveas(gcf,'~/Dropbox/FSTD/Manuscripts/FSTD-Response/Submission-2/Figures/Fig-S1/Fig-S1.pdf')
saveas(gcf,'~/Dropbox/FSTD/Manuscripts/FSTD-Response/Submission-2/Figures/Fig-S1/Fig-S1.fig')
