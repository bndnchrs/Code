

load('/Users/Horvat/Research/FSTD-Code/Output/For-Paper/All-Proc/Example.mat','FSTD','DIAG');

intervalx = [5 50 150];
intervalend = [50 150 1500];

intmult = linspace(.75,1.25,5); 

beg_int_vals = [5 7.5 10; 
    25 50 75; 
    100 150 300]; 

corr_beg = [1 2 2]; 

end_int_vals = [20 40 50; 
    100 150 200; 
    500 1000 1500]; 

%%
for ind = 1:3
    
    for beg_int_ind = 1:3
        
        for end_int_ind = 1:3
            
            clearvars -except FSTD DIAG *_ind ind save* int* *_int_*
            
            
            
            beg_int = beg_int_vals(ind,beg_int_ind);
            end_int = end_int_vals(ind,end_int_ind);
            
            
           
            [~,b] = find(FSTD.Rint > beg_int,1);
            [~,c] = find(FSTD.Rint > end_int,1);
            [~,d] = find(DIAG.FSTD.conc(2:end) < .01,1);
            
            if b < c
                
                if isempty(d)
                    d = length(DIAG.FSTD.conc);
                end
                %%
                
                
                Rinert = FSTD.Rint(b:c);
                FSDinert = DIAG.FSTD.psi(b:c,:,1:d);
                % dAinert = DIAG.FSTD.dA(b:c,:,1:d);
                FSDinert = squeeze(sum(bsxfun(@times,FSDinert,FSTD.dH),2));
                dRinert = FSTD.dR(b:c)';
                % FSDinert = bsxfun(@rdivide,FSDinert, sum(FSDinert,1));
                
                
                FSDinertlog = log10(FSDinert + eps);
                
                Ninert = bsxfun(@rdivide,FSDinert,(pi * (Rinert').^2));
                Ntot = sum(bsxfun(@times,Ninert,dRinert),1);
                
                coeff = zeros(1,size(FSDinert,2));
                coeff2 = coeff;
                
                
                for i = 1:size(FSDinert,2)
                    
                    c_a(i) = sum(FSDinert(:,i).*dRinert,1);
                    c_p(i) = 2*sum(FSDinert(:,i).*dRinert.*FSTD.Rint(b:c).^(-1)');
                    [p,S] = polyfit(log10(Rinert),FSDinertlog(:,i)',1);
                    coeff(i) = -p(1);
                    coeff2(i) = -p(2);
                    
                    FSD2 = Rinert .^(-coeff(i));
                    N2 = FSD2 ./ (pi * (Rinert).^2);
                    Ntot2 = sum(N2.*dRinert');
                    %  FSD2 = FSD2 / sum(FSD2);
                    
                    %mfs{ind,1}(i) = sum(Ninert(:,i).*dRinert'.*2*pi*Rinert') / Ntot(i);
                    % mfs{ind,2}(i) = sum(2*FSD2.*dRinert./Rinert) / sum(FSD2.*dRinert ./ (pi*Rinert.^2));
                    
                    %   sum(bsxfun(@times,2./Rinert',FSD2(:,i))) ./ sum(bsxfun(@times,1./(pi*Rinert.^2)',FSDinert));
                    
                    
                    
                    
                end
                
                % This is the power-law decay of this function, beta
                plaw = coeff(1:d);
                
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
                
                save_plaw{ind,beg_int_ind,end_int_ind,1} = plaw;
                save_plaw{ind,beg_int_ind,end_int_ind,2} = horvat_alpha;
                
                
            end
            
        end
        
    end
    
end
%%
clines = [
27,158,119
217,95,2
117,112,179]/256;


corr_beg = [1 2 2]; 
corr_end = [3 2 3]; 

for i = 1:3
    
    Ax{i} = subplot(2,3,i)
    hold on
    
    
    for j = 1:3
        
        plotter = save_plaw{i,corr_beg(i),j,1};
        
        
        if length(plotter) > 0
            
            
            plot([0 FSTD.time]/86400,plotter,'color',clines(j,:),'linewidth',2)
            
        end
        
    end
    
    Ax{3+i} = subplot(2,3,3+i)
    hold on
    
    
    for j = 1:3
        
        plotter = save_plaw{i,j,corr_end(i),1};
        
        
        if length(plotter) > 0
            
            
            plot([0 FSTD.time]/86400,plotter,'color',clines(j,:),'linewidth',2)
            
        end
        
    end
    
    
    
end

titles = {'Regime I','Regime II','Regime III'}; 

ylabs = {'Vary r_1','','','Vary r_2'}; 

for i = 1:6
    subplot(2,3,i)
    grid on
    box on
    set(gca,'ydir','normal','layer','top','fontname','helvetica','fontsize',18)
    if i > 3
        
        xlabel('Time (days)')
        
        legend(['r1 = ' num2str(beg_int_vals(i-3,corr_beg(i-3))) ', r2 = ' num2str(end_int_vals(i-3,1))],['r1 = ' num2str(beg_int_vals(i-3,corr_beg(i-3))) ', r2 = ' num2str(end_int_vals(i-3,2))],['r1 = ' num2str(beg_int_vals(i-3,corr_beg(i-3))) ', r2 = ' num2str(end_int_vals(i-3,3))]); 

        
    else
        
        legend(['r1 = ' num2str(beg_int_vals(i,1)) ', r2 = ' num2str(end_int_vals(i,corr_end(i)))], ['r1 = ' num2str(beg_int_vals(i,2)) ', r2 = ' num2str(end_int_vals(i,corr_end(i)))],['r1 = ' num2str(beg_int_vals(i,3)) ', r2 = ' num2str(end_int_vals(i,corr_end(i)))]); 

    end
    
    if i < 4
        title(titles{i})
    end
    
    if i == 1 || i == 4
        ylabel(ylabs{i}); 
    end
    
    xlim([0 90])
    
    ylim([0 8])
    
end

for i = 1:3
    
    subplot(2,3,i)
    hold on
    
    
    for j = 1:3
        
        plotter = save_plaw{i,corr_beg(i),j,2};
        
        
        if length(plotter) > 0
            
            
            plot([0 FSTD.time]/86400,plotter,'--','color',clines(j,:),'linewidth',2)
            
        end
        
    end
    
    subplot(2,3,3+i)
    hold on
    
    
    for j = 1:3
        
        plotter = save_plaw{i,j,corr_end(i),2};
        
        
        if length(plotter) > 0
            
            
            plot([0 FSTD.time]/86400,plotter,'--','color',clines(j,:),'linewidth',2)
            
        end
        
    end
    
    
    
end



pos = [12 8]; 
 set(gcf,'windowstyle','normal','position',[0 0 pos],'paperposition',[0 0 pos],'papersize',pos,'units','inches','paperunits','inches');
  set(gcf,'windowstyle','normal','position',[0 0 pos],'paperposition',[0 0 pos],'papersize',pos,'units','inches','paperunits','inches');
tightfig

%%
letter = {'(a)','(b)','(c)','(d)','(e)','(f)','(g)','(h)'};
delete(findall(gcf,'Tag','legtag'))
for i = 1:6
    posy = get(Ax{i},'position');
    annotation('textbox',[posy(1)-.06 posy(2)+posy(4)+.025 .025 .025], ...
        'String',letter{i},'LineStyle','none','FontName','Helvetica', ...
        'FontSize',14,'Tag','legtag');
end

%%

saveas(gcf,'Fig_S2.pdf');
saveas(gcf,'Fig_S2.fig'); 
  