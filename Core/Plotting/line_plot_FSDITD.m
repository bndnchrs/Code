    %% Plot FSD as line plots
    
    
    if FSTD.i == OPTS.nt
        
        OPTS.ax_FSD_line = subplot(323);
        
        
        hold on
        
        per = log10(squeeze(sum(bsxfun(@times,DIAG.FSTD.psi(:,:,inds),FSTD.dA),2)+eps));
        
        for i = 1:length(inds)
            
            plot(FSTD.Rint,per(:,i),'color',colplots(i,:));
            
            str{i} = ['Day: ' num2str(round(DIAG.FSTD.time(inds(i))/86400))];
            
        end
        xlim([FSTD.Rint(1) FSTD.Rint(end)])
        
        
        llim = floor(log10(1/length(FSTD.Rint)) - 1);
        
        grid on
        box on
        xlabel('Floe Size')
        title('log_{10} of FSD(r)dr')
        ylabel('log10(m^2/m^2)')
        set(gca,'ydir','normal','layer','top','fontname','helvetica','fontsize',14)
        % legend(str)
        set(gca,'ylim',[llim 0])
        
        shading interp
        
        
        
        
        if WAVES.DO
            
            
            [a,b] = max(WAVES.spec);
            
            % This is R which is half the wavelength of the peak period
            R_longterm = FSTD.Rmid(b);
            scatter(R_longterm,.5,200,'k')
            
        end
        
    end
    
    %% Plot ITD as line plots
    
    
    if FSTD.i == OPTS.nt
        
        OPTS.ax_ITD_line = subplot(324);
        
        hold on
        
        for i = 1:length(inds)
            
            plot(FSTD.Hmid_i,log10(squeeze(sum(bsxfun(@times,DIAG.FSTD.psi(:,:,inds(i)),FSTD.dA),1)+eps)),'color',colplots(i,:));
            
        end
        
        llim = floor(log10(1/length(FSTD.Hmid)) - 1);
        
        hold on
        grid on
        box on
        xlabel('Ice Thickness')
        title('log_{10} of ITD(h)dh')
        ylabel('log10(m^2/m^2)')
        set(gca,'ydir','normal','layer','top','fontname','helvetica','fontsize',14)
        legend(str)
        xlim([FSTD.Hmid_i(1) FSTD.Hmid_i(end)])
        set(gca,'ylim',[llim 0])
        
        shading interp
        
        
    end