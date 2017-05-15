
close all

xlimmer = [0 12];

cplots = [228,26,28
    55,126,184
    77,175,74
    152,78,163
    255,127,0
    255,255,51
    166,86,40]/256;

ADvals = [integrate_FSTD(ADVECT.FSTD_in,1,FSTD.dA,0),  ...
    integrate_FSTD(ADVECT.FSTD_in,FSTD.Hmid,FSTD.dA,0), ...
    integrate_FSTD(ADVECT.FSTD_in,FSTD.Hmid,FSTD.dA,1), ...
    integrate_FSTD(ADVECT.FSTD_in./(pi * FSTD.meshR.^2),FSTD.Rmid',FSTD.dA,1), ...
    integrate_FSTD(ADVECT.FSTD_in,FSTD.Rmid',FSTD.dA,1)];

vals = {DIAG.FSTD.conc, ...
    DIAG.FSTD.Vtot, ...
    DIAG.FSTD.Hmean, ...
    DIAG.FSTD.Rmeannum, ...
    DIAG.FSTD.Rmeanarea};

legstr = {'Concentration','Volume','Mean Thickness','Mean Floe Size (N)','Mean Floe Size (f)'};

Ax{5} = subplot(122);
hold on

ylabs = {'%','m^3/m^3','m','m','m','%'};
titles = {'Ice Concentration','Ice Volume','Mean Thickness','Mean Floe Size','Approach to Advected Value'};

ploti = [1 2 3 4 4 5]; 
splot = [1 2 5 6 6 7];

abciss = DIAG.FSTD.time/86400;

subplot(122)
xlabel('Time (days)')

for i = 1:length(vals)
    
    subplot(122)
    
    if i ~= 2 && i~=5
        
        plot(abciss,1 - (vals{i} - ADvals(i))/(vals{i}(1) - ADvals(i)),'linewidth',1,'color',cplots(i,:));
        
    else
        
        plot(abciss,1 - (vals{i} - ADvals(i))/(vals{i}(1) - ADvals(i)),'--','linewidth',1,'color',cplots(i,:));
        
    end
    
    
    
    Ax{ploti(i)} = subplot(2,4,splot(i));
    
    if i~=5
    
    plot(abciss,vals{i},'linewidth',1,'color',cplots(i,:));
    
    hold on
    
    plot(abciss,0*(vals{i}) + ADvals(i),'--','linewidth',1,'color','k');
    
    else
        
        plot(abciss,vals{i},'--','linewidth',1,'color',cplots(i,:));
    
    end
    
    title(Ax{ploti(i)},titles{ploti(i)})
    ylabel(Ax{ploti(i)},ylabs{ploti(i)})
   
    if i == 3 || i == 4 || i == 5
        xlabel('Time (days)')
    else
        set(Ax{i},'xticklabel',{});
    end
    
    
end

subplot(122)

legend(legstr)

%%
tightfig

pos = [8 3.5];
set(gcf,'windowstyle','normal','position',[0 0 pos],'paperposition',[0 0 pos],'papersize',pos,'units','inches','paperunits','inches');
set(gcf,'windowstyle','normal','position',[0 0 pos],'paperposition',[0 0 pos],'papersize',pos,'units','inches','paperunits','inches');

xlabel(Ax{4},'Time (days)')

drawnow


letter = {'(a)','(b)','(c)','(d)','(e)','(f)','(g)','(h)'};

hAllAxes = findobj(gcf,'type','axes');
hLeg = findobj(hAllAxes,'tag','legend');
hAxes = setdiff(hAllAxes,hLeg); % All axes which are not

delete(findall(gcf,'Tag','legtag'))

% 


% tightfig

for i = 1:length(Ax)
    posy = get(Ax{i},'position');
    annotation('textbox',[posy(1)+.001 posy(2)+posy(4)-.025 .025 .025], ...
        'String',letter{i},'LineStyle','none','FontName','Helvetica', ...
        'FontSize',14,'Tag','legtag');
    grid(Ax{i},'on')
    box(Ax{i},'on')
    set(Ax{i},'xlim',xlimmer)

    set(Ax{i},'ydir','normal','layer','top','fontname','helvetica','fontsize',12)
end



saveas(gcf,'~/Dropbox/FSTD/Manuscripts/FSTD-Response/Figures/Fig-1/Fig-1.pdf')
saveas(gcf,'~/Dropbox/FSTD/Manuscripts/FSTD-Response/Figures/Fig-1/Fig-1.fig')

