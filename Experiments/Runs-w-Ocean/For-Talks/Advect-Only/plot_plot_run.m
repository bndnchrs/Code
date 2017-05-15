load /Users/Horvat/Research/FSTD-Code/Output/Simple_Advection/Example.mat


close all

%%

timer = DIAG.FSTD.time/86400;

xlimmer = [0 max(timer)];

xlimmer = [0 12];

cplots = [228,26,28
    55,126,184
    77,175,74
    152,78,163
    255,127,0
    255,255,51
    166,86,40]/256;

ADvals = [100*integrate_FSTD(ADVECT.FSTD_in,1,FSTD.dA,0),  ...
    integrate_FSTD(ADVECT.FSTD_in,FSTD.Hmid,FSTD.dA,0), ...
    integrate_FSTD(ADVECT.FSTD_in,FSTD.Hmid,FSTD.dA,1), ...
    integrate_FSTD(ADVECT.FSTD_in./(pi * FSTD.meshR.^2),FSTD.Rmid',FSTD.dA,1), ...
    integrate_FSTD(ADVECT.FSTD_in,FSTD.Rmid',FSTD.dA,1)];

vals = {DIAG.FSTD.conc*100, ...
    DIAG.FSTD.Vtot, ...
    DIAG.FSTD.Hmean, ...
    DIAG.FSTD.Rmeannum, ...
    DIAG.FSTD.Rmeanarea};

legstr = {'Ice Concentration','Ice Volume','Mean Thickness','Mean Floe Size (N)','Mean Floe Size (f)'};


ylabs = {'%','m^3/m^3','m','m','m','%'};
titles = {'Ice Concentration','Ice Volume','Mean Thickness','Mean Floe Size','Approach to Advected Value'};

names = {'c','V','H','r'}; 

skipper = 6; 

for i = 1:length(vals)
    
    close
    
    plot(timer,0*timer + ADvals(i),'--k','linewidth',2);
    
    create_movie_FD(names{i},0,[0 0 6 3])

    hold on
    
    grid on
    box on
    xlim(xlimmer)
    set(gca,'ydir','normal','layer','top','fontname','helvetica','fontsize',14)
    title(titles{i})
    xlabel('Time (days)')
    ylim([0 max(max(vals{i}),ADvals(i))])
    
    ylabel(ylabs{i})
    ind = 0; 
    
    for jj = 1:skipper:length(timer)
        
        ind = ind + 1; 
        
       pl =  plot(timer(1:jj),vals{i}(1:jj),'linewidth',2,'color',cplots(i,:));
       sc = scatter(timer(jj),vals{i}(jj),200,cplots(i,:),'filled'); 
                  create_movie_FD(names{i},ind,[0 0 6 3])

       delete(pl);
       delete(sc); 
       

    end
    
    create_movie_FD('close',jj,[0 0 6 3])

    
end
