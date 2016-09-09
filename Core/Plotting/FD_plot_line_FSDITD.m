%% Plot FSD as line plots
cols = [228,26,28
    55,126,184
    77,175,74
    152,78,163
    255,127,0
    155,155,155]/256;

if FSTD.i == OPTS.nt
    
    PLOTS.ax_FSD_line = subplot(323);
    
    
    hold on
    
    per = squeeze(sum(bsxfun(@times,DIAG.FSTD.psi(:,:,inds),DIAG.FSTD.dA(:,:,inds)),2)+eps);
    
    per = bsxfun(@rdivide,per,DIAG.FSTD.conc(inds));
    
    % per = log10(per);
    
    for i = 1:length(inds)
        
        semilogy(DIAG.FSTD.R,per(:,i),'color',colplots(i,:),'linewidth',2);
        
        str{i} = ['Day: ' num2str(round(DIAG.FSTD.time(inds(i))/86400))];
        
    end
    
    xlim([DIAG.FSTD.R(1) DIAG.FSTD.R(end)])
    
    
    llim = floor(log10(1/length(DIAG.FSTD.R)) - 1);
    
    grid on
    box on
    xlabel('Floe Size')
    title('log_{10} of FSD(r)dr')
    ylabel('log10(m^2/m^2)')
    set(gca,'ydir','normal','layer','top','fontname','helvetica','fontsize',14)
    % legend(str)
    % set(gca,'ylim',[llim 0])
    
    shading interp
    
    plotter = sum(ADVECT.FSTD_in .* FSTD.dA,2) + eps;
    plotter = plotter ./ integrate_FSTD(ADVECT.FSTD_in,1,FSTD.dA,0);
    
    semilogy(DIAG.FSTD.R,plotter,'--','color',cols(3,:),'linewidth',2)
    
    set(gca,'xscale','log')
    
    set(gca,'yscale','log')
    set(gca,'ylim',[1e-6 1])
    
    
    if WAVES.DO
        
        
        [a,b] = max(WAVES.spec);
        
        % This is R which is half the wavelength of the peak period
        R_longterm = FSTD.Rmid(b);
        scatter(R_longterm,.5,200,'k')
        
    end
    
    
    
end

%% Plot ITD as line plots


if FSTD.i == OPTS.nt
    
    PLOTS.ax_ITD_line = subplot(324);
    
    hold on
    
    per = squeeze(sum(bsxfun(@times,DIAG.FSTD.psi(:,:,inds),DIAG.FSTD.dA(:,:,inds)),1)+eps);
    
    per = bsxfun(@rdivide,per,DIAG.FSTD.conc(inds));
    
    %  per = log10(per);
    
    for i = 1:length(inds)
        
        semilogy(DIAG.FSTD.Hmid(:,i),per(:,i),'color',colplots(i,:),'linewidth',2);
        
    end
    
    llim = floor(log10(1/length(DIAG.FSTD.H)) - 1);
    
    hold on
    grid on
    box on
    xlabel('Ice Thickness')
    title('log_{10} of ITD(h)dh')
    ylabel('log10(m^2/m^2)')
    set(gca,'ydir','normal','layer','top','fontname','helvetica','fontsize',14)
    legend(str)
    xlim([DIAG.FSTD.H(1) DIAG.FSTD.H(end)])
    % set(gca,'ylim',[llim 0])
    
    plotter = sum(ADVECT.FSTD_in .* FSTD.dA,1) + eps;
    plotter = plotter ./ integrate_FSTD(ADVECT.FSTD_in,1,FSTD.dA,0);
    
    semilogy(FSTD.Hmid,plotter,'--','color',cols(3,:),'linewidth',2)
    
    set(gca,'yscale','log')
    set(gca,'ylim',[1e-6 1])
    
    shading interp
    
    
end

