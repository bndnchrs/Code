%% Update_local_variables

% This script updates all of the variables that change on the order of
% one internal timestep but are not reset, mostly counting variables



% Marginal Distributions
FSTD.FSD = sum(bsxfun(@times,FSTD.psi,FSTD.dH),2);
FSTD.ITD = sum(bsxfun(@times,FSTD.psi,FSTD.dR'),1);

% How much time left to go in the timestep
OPTS.dt_sub = OPTS.dt_sub - OPTS.dt_temp;

% Counter of sub-cycles, both per global timestep and in totality
OPTS.numSC = OPTS.numSC + 1;
OPTS.totnum = OPTS.totnum + 1;

%%

% How much volume is in the largest floe class.
FSTD.V_max = FSTD.V_max + OPTS.dt_temp*FSTD.dV_max;

% To ensure the code converges, if this value is too small, we just set it
% to zero. If it is that small, then the maximum thickness is just equal to
% the initial value we set for the maximum thickness. 
if abs(FSTD.V_max) < 1e-8 % && sum(FSTD.ITD(1:end-1).*FSTD.H) < eps
    FSTD.V_max = 0; 
    FSTD.H_max = FSTD.H_max_i;
    FSTD.A_max = 0; 
    if integrate_FSTD(FSTD.ITD,FSTD.H,FSTD.dH,0) < eps
    FSTD.psi = 0*FSTD.psi; 
    end
else
    FSTD.A_max = integrate_FSTD(FSTD.psi(:,end),1,FSTD.dA(:,end),0);
end

% if FSTD.H_max > 10 && FSTD.A_max < 1e-6
%     disp('H_max is too huge. This is likely just an I.C. issue. Deleting the ice area there.')
%     disp('If this message re-appears a bunch, need to examine the code.')
%     disp(['Timestep ' num2str(FSTD.i)])    
%     FSTD.psi(:,end) = 0; 
%     FSTD.A_max = 0; 
%     FSTD.V_max = 0; 
% end


% The mean thickness is equal to the sum of the 
FSTD.Hmean = integrate_FSTD(FSTD.psi,FSTD.Hmid,FSTD.dA,1);

if FSTD.Hmean == 0
    FSTD.Hmean = OPTS.h_p;
end

Ameps = 0;

if FSTD.A_max == 0
    
    Ameps = eps;
    FSTD.V_max = 0; 
    
end

% The top category thickness is simply the ratio of the volume to the
% concentration

% if FSTD.i == 42
%     disp('ok')
% end

FSTD.H_max = FSTD.V_max / (Ameps + FSTD.A_max);

if FSTD.H_max == 0
    FSTD.H_max = FSTD.H_max_i; 
end

FSTD.Hmid(end) = FSTD.H_max; 
%% We need to correct for changes to the thickest class of ice. 

% If the thickest ice category decreases in thickness so that it becomes
% thinner than the right-most boundary of the thickest non-variable ice
% thickness category, we have to adjust where the ice is. 

if size(FSTD.H,2) > 1 && FSTD.H_max < FSTD.H(end)
 %%   
    % This is the initial total ice volume and concentration
    V_init = integrate_FSTD(FSTD.psi,FSTD.Hmid,FSTD.dA,0); 
    c_init = integrate_FSTD(FSTD.psi,1,FSTD.dA,0); 
    
    % This is the total initial ice volume and concentration in each of the
    % thickest floe size categories. 
    V_i = FSTD.psi(:,end) .* FSTD.dA(:,end) * FSTD.H_max; 
    c_i = FSTD.psi(:,end) .* FSTD.dA(:,end); 
    
    
    % When the ice thickness decreases to be within the thickest fixed
    % category, i.e. when FSTD.H_max < FSTD.H(end), we have an issue, as dH
    % is defined to be negative. Therefore we need to put this ice into
    % that category. 
    
    % We do this by ensure that given an initial volume V_i of ice in the
    % thickest category, and an initial concentration c_i, this is
    % conserved, and re-distributed to a thinner ice thickness category.
    % After this happens, the grid is restored to its initial values, i.e.
    % FSTD.H_max = FSTD.H_max_i, and FSTD.dA = FSTD.dA_i. 

    % This can be done with a simple algebraic relation. If the ice is
    % moved to the thickness category (ind), we have
    % V_i = H(ind) * c(ind) + H(end) * c(end) 
    % c_i = c(ind) + c(end). 
    
    % This can be solved analytically. With c_1 = c(ind) and c_2 = c(end),
    % and with H_2 = H(end) = H_max_i, and H_1 = H(ind), 

    % c_1 = (c_i H_2 - V_i) / (H_2 - H_1)
    % c_2 = (V_i - c_i H_1) / (H_2  - H_1)
    
    % First we find the thickness category to which all the ice will go. 
    % It must be smaller than the currently too-small thickness 
    
    [a,b] = find(FSTD.Hmid < FSTD.H_max);
    ind = b(end); 
    
    % Here are our new thicknesses
    H_1 = FSTD.Hmid_i(ind);
    H_2 = FSTD.H_mid_max_i; 
    
    % Now we compute the new concentrations
    
    c_1 = (c_i * H_2 - V_i) / (H_2 - H_1); 
    c_2 = (V_i - c_i * H_1) / (H_2 - H_1); 
   
    % Now we reset the grid
    FSTD.H_max = FSTD.H_max_i; 
    FSTD.Hmid(end) = FSTD.H_mid_max_i; 
    FSTD.dA = FSTD.dA_i; 
    
    % Now we convert these to FSD values
    psi_1 = c_1 ./ FSTD.dA(:,ind); 
    psi_2 = c_2 ./ FSTD.dA(:,end); 
    
    FSTD.psi(:,end) = 0; 
    
    diffpsi = 0*FSTD.psi; 
    diffpsi(:,ind) = psi_1; 
    diffpsi(:,end) = psi_2 - FSTD.psi(:,end);
    
    FSTD.psi = FSTD.psi + diffpsi; 
    FSTD.V_max = sum(FSTD.psi(:,end) .* FSTD.dA(:,end)) * FSTD.H_max; 
    
    V_f = integrate_FSTD(FSTD.psi,FSTD.Hmid,FSTD.dA,0); 
    c_f = integrate_FSTD(FSTD.psi,1,FSTD.dA,0); 

    if abs(V_init - V_f) > 1e-4
        FSTD.eflag = 1;
        disp('Exchange of Volume at Thickest Size is a Problem')
    end
    
    if abs(c_init - c_f) > 1e-4
        FSTD.eflag = 1;
        disp('Exchange of Concentration at Thickest Size is a Problem')
    end
    
    
end



% now update the time!

FSTD.time_now = FSTD.time_now + OPTS.dt_temp;
