    %% Plot FSD
    
    set(gcf,'currentaxes',OPTS.ax_FSD);
    
    
    if FSTD.i == 1
        
        
        OPTS.p1 = pcolor([0 FSTD.time/86400],FSTD.Rint,log10(squeeze(sum(bsxfun(@times,DIAG.FSTD.psi,DIAG.FSTD.dA),2)+eps)));
        hold on
        grid on
        box on
        ylabel('Floe Size')
        title('log_{10} of FSD(r)dr')
        xlabel('Time')
        set(gca,'ydir','normal','layer','top','fontname','helvetica','fontsize',14)
        
        ylim([FSTD.Rint(1) FSTD.Rint(end)]);
        
        
        set(gca,'clim',[-2 0])
        shading interp
        colorbar
        
    else
        
        
        
        delete(OPTS.p1)
        
        % determine interpolation grid: from min to max in equal steps
        ax1 = [0 FSTD.time/86400];
        ax1 = ax1(1:skip_contour:end);
        ax2 = FSTD.Rint;
        
        Ax1v = linspace(ax1(1),ax1(end),numel(ax1));
        Ax2v = linspace(ax2(1),ax2(end),numel(ax2));
        
        [ax1,ax2] = meshgrid(ax1,ax2);
        [Ax1,Ax2] = meshgrid(Ax1v,Ax2v);
        
                    
        plotter = log10(squeeze(sum(bsxfun(@times,DIAG.FSTD.psi(:,:,1:skip_contour:end),DIAG.FSTD.dA(:,:,1:skip_contour:end)),2)+eps));
        
        int = interp2(ax1,ax2,plotter,Ax1,Ax2);
        
        OPTS.p1 = pcolor(Ax1v,Ax2v,int);
        shading interp
        axis xy
        
        
        %  OPTS.p1 = pcolor([0 FSTD.time/86400],FSTD.Rint,log10(squeeze(sum(bsxfun(@times,DIAG.FSTD.psi,DIAG.FSTD.dA),2)+eps)));
        %  shading interp
        set(gca,'clim',[-4 0])
        
        if FSTD.i == OPTS.nt && WAVES.DO
            
            
            [a,b] = max(WAVES.spec);
            
            % This is R which is half the wavelength of the peak period
            R_longterm = FSTD.Rmid(b);
            plot(FSTD.time/86400,0*FSTD.time + R_longterm,'-k','linewidth',2)
            
        end
        
        
        if FSTD.i == OPTS.nt
            
            imcontour(ax1,FSTD.Rint,log10(squeeze(sum(bsxfun(@times,DIAG.FSTD.psi(:,:,1:skip_contour:end),DIAG.FSTD.dA(:,:,1:skip_contour:end)),2)+eps)),[-4:1:0],'--k','showtext','on');
            % set(gca,'clim',[-4 0])
            
        end
        
        
        
    end
    
    %% Plot ITD
    
    set(gcf,'currentaxes',OPTS.ax_ITD);
    
    
    if length(FSTD.Hmid) > 1
        if FSTD.i == 1
            
            
            
            
            
            OPTS.p2 = pcolor([0 FSTD.time/86400],FSTD.Hmid,log10(squeeze(sum(bsxfun(@times,DIAG.FSTD.psi,DIAG.FSTD.dA),1)+eps)));
            
            hold on
            grid on
            box on
            ylabel('Floe Size')
            title('log_{10} of ITD(h)dh')
            xlabel('Time')
            set(gca,'ydir','normal','layer','top','fontname','helvetica','fontsize',14)
            
            ylim([FSTD.Hmid(1) FSTD.Hmid_i(end)]);
            
            
            set(gca,'clim',[-2 0])
            %       shading interp
            colorbar
            
        else
            
            
            delete(OPTS.p2)
            
            % determine interpolation grid: from min to max in equal steps
            ax1 = [0 FSTD.time/86400];
            ax1 = ax1(1:skip_contour:end);
            
            ax2 = FSTD.Hmid_i;
            
            Ax1v = linspace(ax1(1),ax1(end),numel(ax1));
            Ax2v = linspace(ax2(1),ax2(end),numel(ax2));
            
            [ax1,ax2] = meshgrid(ax1,ax2);
            [Ax1,Ax2] = meshgrid(Ax1v,Ax2v);
            
            
            
            plotter = log10(squeeze(sum(bsxfun(@times,DIAG.FSTD.psi(:,:,1:skip_contour:end),DIAG.FSTD.dA(:,:,1:skip_contour:end)),1)+eps));
            
            int = interp2(ax1,ax2,plotter,Ax1,Ax2);
            
            OPTS.p2 = pcolor(Ax1v,Ax2v,int);
            shading interp
            
            axis xy
            
            
            set(gca,'clim',[-4 0])
            
            
            if FSTD.i == OPTS.nt
                
                imcontour(ax1,ax2,log10(squeeze(sum(bsxfun(@times,DIAG.FSTD.psi(:,:,1:skip_contour:end),DIAG.FSTD.dA(:,:,1:skip_contour:end)),1)+eps)),[-4:1:0],'--k','showtext','on');
                set(gca,'clim',[-4 0])
                
            end
        end
    end
    
