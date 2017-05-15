clear


load('/Users/Horvat/Research/FSTD-Code/Output/For-Paper/Mech-Advect/Example.mat','DIAG','FSTD');

% The full FSTD
psi = bsxfun(@times,DIAG.FSTD.psi,FSTD.dH); % Convert to FSD
psi = squeeze(sum(psi,2)); % This is the FSD f(r)

conc = sum(bsxfun(@times,psi,FSTD.dR'),1);

nt = find(conc == 0,1)-1; % Number of nonzero concentrations
if isempty(nt)
    nt = length(conc); 
end
% The floe size categories
psi = psi(1:end-1,:); 
R = FSTD.Rmid(1:end-1);
dR = FSTD.dR(1:end-1);

% We will generate synthetic data to fit a power law to
Npts = 1e6; % number of floes to "measure"

KL_div = zeros(round(nt/24),5); 
KL_div_symm = KL_div; 
ET_distance = KL_div; 
CH_distance = KL_div; 


%% 
i = 0; 

for ind = 1:24:nt   
    
    i = i + 1; 
    % Take the distribution at each timestep. This is 
    % f(r), not the probability distribution, i.e. it is not normalized.
    % f(r) is the thing that decays as a "power law"
    FSD = psi(:,i);
    
    %% If the FSD is nonzero somewhere
    if sum(FSD) ~= 0
        %% We will create a random sample of floe sizes from the distribution 
        % f(r), and find the maximum likelihood estimate of the power law
        % that fits this sample 
        
        % Convert f(r) to f(r) dr - the binned probability distribution
        % Now this is a discrete probability distribution which sums to 1
        f_pdf = FSD.*dR' / sum(FSD.*dR');
        
        % Convert the distribution into one on the unit interval
        unit_f = cumsum(f_pdf)/sum(f_pdf);
        
        % Make many "observations"
        r = rand(Npts,1);
        % Bin these observations over the categories with their assigned [0,1]
        % values
        [p,VC_data] = histc(r,unit_f);
        % Convert this histogram into a synthetic dataset
        
        f_vals = (R(VC_data+1)); % our synthetic data
        
        
        %% Using the Virkar and Clauset MLE, compute the power law decay
        % coefficient alpha, and the beginning of the power law tail xmin
       
        [fit_alpha(i), xmin(i), L(i)]=plfit(f_vals);
        
        tail_beg = find(R > xmin(i),1)-1; % Beginning of tail

        
        %% Next we want to form the pdf from this best fit and compare
        % the tails
        
        % Only consider the tail of the FSD above x_min
        f_tail = FSD'; 
        f_tail(R<xmin(i)) = 0;
        % Normalize to 1
        f_tail = f_tail / sum(f_tail.*dR);
        
        % Form the distribution using the best fitting power law VC_alpha(i)
        g_tail = (R).^(-fit_alpha(i));
        % And ignore all points below the minimum size
        g_tail(R<xmin(i)) = 0;
        % Normalize the pdf
        g_tail = g_tail / sum(g_tail.*dR);
        
        %% Next we compare alternative hypotheses for the fitting and
        % compute the distances for these hypotheses as well to compare
        
        %% Fit using a least squares regression to the plot
    
        Rinert = R(tail_beg:end); 
        FSDinert = FSD(tail_beg:end)'; 
        p = polyfit(log10(Rinert),log10(FSDinert),1);
        coeff(i) = -p(1);
    
        gls_pdf = (R).^(-coeff(i));
        gls_tail = gls_pdf; 
        gls_tail(R<xmin(i)) = 0; 
        gls_tail = gls_tail/sum(gls_tail.*dR); 
    
        %% Fit f_vals to a generalized pareto distribution
        
        % take all measurements in the tail
        tail_vals = f_vals(f_vals >= xmin(i));
        
        % Fit the tail values, offset from zero
        gp_par = gpfit(tail_vals - min(tail_vals) + eps); 
        
        % Compute the pdf from these values, restored by the offset
        gp_pdf = gppdf(R,gp_par(1),gp_par(2),min(tail_vals));
       
        % Only concern ourselves
        gp_tail = gp_pdf; 
        gp_tail(f_tail == 0) = 0; 
        gp_tail = gp_tail / sum(gp_tail.*dR); 
         
        
        %% Compute an exponential fit
        
        ex_par = expfit(tail_vals-min(tail_vals)+eps); 
        ex_pdf = exppdf(R-min(tail_vals),ex_par); 
        ex_tail = ex_pdf; 
        ex_tail(f_tail == 0) = 0; 
        ex_tail = ex_tail / sum(ex_tail.*dR); 
        
        %% Compute the GLV fit
        fo = fitoptions('Method','Nonlinearleastsquares','StartPoint',1);
        ft = fittype('x^(-1-a)*exp((1-a)/x)','independent','x','options',fo);
        
        [curve2,gof2] = fit(R(tail_beg),f_tail(tail_beg)',ft);
         
        glv_pdf = curve2(R)./dR'; 
        glv_tail = glv_pdf; 
        glv_tail(f_tail == 0) = 0; 
        glv_tail = glv_tail' / sum(glv_tail.*dR'); 
        

        
        %% Now we will compute the difference between the two distributions 
        %% And the hypothesized distributions in a few ways
        
        %% The KL distance
        
        KL_div(i,1) = KLdiv(g_tail,f_tail);
        KL_div(i,2) = KLdiv(gp_tail,f_tail); 
        KL_div(i,3) = KLdiv(ex_tail,f_tail); 
        KL_div(i,4) = KLdiv(glv_tail,f_tail); 
        KL_div(i,5) = KLdiv(gls_tail,f_tail); 

        
        %% The Symmetric KL distance
        
        KL_div_symm(i,1) = .5*KLdiv(g_tail,f_tail) + .5*KLdiv(f_tail,g_tail);
        KL_div_symm(i,2) = .5*KLdiv(gp_tail,f_tail) + .5*KLdiv(f_tail,gp_tail);
        KL_div_symm(i,3) = .5*KLdiv(ex_tail,f_tail) + .5*KLdiv(f_tail,ex_tail);
        KL_div_symm(i,4) = .5*KLdiv(glv_tail,f_tail) + .5*KLdiv(f_tail,glv_tail);
        KL_div_symm(i,5) = .5*KLdiv(gls_tail,f_tail) + .5*KLdiv(f_tail,gls_tail);

        
        %% The ET distance
        
        ET_distance(i,1) = ET_dist(g_tail,f_tail); 
        ET_distance(i,2) = ET_dist(gp_tail,f_tail); 
        ET_distance(i,3) = ET_dist(ex_tail,f_tail); 
        ET_distance(i,4) = ET_dist(glv_tail,f_tail); 
        ET_distance(i,5) = ET_dist(gls_tail,f_tail); 

        
        
        %% The CH distance
        
        CH_distance(i,1) = CH_dist(g_tail,f_tail); 
        CH_distance(i,2) = CH_dist(gp_tail,f_tail); 
        CH_distance(i,3) = CH_dist(ex_tail,f_tail); 
        CH_distance(i,4) = CH_dist(glv_tail,f_tail); 
        CH_distance(i,5) = CH_dist(gls_tail,f_tail); 
        
        %% The HT distance
        
        
        
        
    end
    
    
    
    
end


%% Do some plotting
cla

subplot(221)

time = FSTD.time(1:i)/86400; 

plot(time,KL_div,'linewidth',1)
title('KL Divergence')


subplot(222)

plot(time,KL_div_symm,'linewidth',1)
title('Symmetric KL')

subplot(223)
plot(time,ET_distance,'linewidth',1)
title('ET Distance')

subplot(224)

plot(time,CH_distance,'linewidth',1)
title('CH Distance')

legend('Power Law','Pareto','Exponential','GLV','Naive Least Squares')

drawnow

disp('done')

for i = 1:4
    subplot(2,2,i)
    hold on
    xlim([time(1) time(nt)])
    ylim([0 max(max(get(gca,'ylim')),1)])
    grid on
    box on
    set(gca,'layer','top','fontname','helvetica','fontsize',18,'ydir','normal');
end

pos = [12 8]
set(gcf,'windowstyle','normal','position',[0 0 pos],'paperposition',[0 0 pos],'papersize',pos,'paperunits','inches','units','inches');
set(gcf,'windowstyle','normal','position',[0 0 pos],'paperposition',[0 0 pos],'papersize',pos,'paperunits','inches','units','inches');