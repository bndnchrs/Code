%% Mech_Timestep
% This routine calculates the tendency at each floe size and thickness
% according to the FD parameterizations, and also updates the
% large-ice-thickness class appropriately.

% There are two parts to the tendency due to mechanical forcing
% First, there is the outgoing term. This is the total fraction at each
% size that goes to other floe sizs.

% Second, there is the incoming term. This is the determination of how the
% first term is distributed amongst larger floes

% Want to make a matrix which is equal to the total area lost at each
% timestep from the combination of r1,h1 and r2,h2.

% The matrices involved are:
% Size R refers to length(FSTD.Rmid)
% Size H refers to length(FSTD.Hmid)
% Size 1 refers to size 1

% 1: diagtwo - RRHH this handles double counting of the index (r,r,h,h)
% 2: pi * Rmid^2 - R - this is the floe size vector
% 3: MECH.Prob_Interact_raft - RR - the probability that floes of size r1
% and r2 will combine
% 4: MECH.gamma_raft - HH - The likelihood of rafting over ridging
% 5: FSTD.Numberdist - RH - The number distribution. This needs to be used
% twice as RH11 and 11RH

% All of these will be permuted to be the size RHRH where the third and
% fourth arguments refer to the second interaction pair.

% When floes move into the largest category, they have the thickness
% FSTD.H_max_i. However, the actual thickness in this category is
% FSTD.H_max. Therefore these floes will have to have the area they
% contribute reduced a bit to account for this. 

%% Rafting interaction

if MECH.rafting
    
    MECH.corr_hmax_raft = FSTD.H_max_i/FSTD.H_max * MECH.hmax_flag_raft ... 
    + 1 * (1 - MECH.hmax_flag_raft); 
    
    MECH.Out_raft = many_mult( ...
        permute(MECH.diagtwo,[1 3 2 4]), ...
        permute(pi*FSTD.Rmid.^2,[1 2 3 4])', ...
        permute(MECH.Prob_Interact_raft,[1 3 2 4]), ...
        permute(MECH.gamma_raft,[3 1 4 2]),...
        permute(FSTD.NumberDist,[1 2 3 4]), ...
        permute(FSTD.NumberDist,[3 4 1 2]) ...
        );
    
    MECH.In_raft_coeff = many_mult( ...
        permute(MECH.diagone,[1 3 2 4]), ...
        permute(MECH.gamma_raft,[3 1 4 2]), ...
        permute(MECH.Kfac_raft,[1 3 2 4]), ...
        permute(MECH.Prob_Interact_raft,[1 3 2 4]), ... 
        permute(FSTD.NumberDist,[1 2 3 4]), ...
        permute(FSTD.NumberDist,[3 4 1 2]), ...
        permute(MECH.corr_hmax_raft,[1 2 3 4]) ...
        );
    
    % Accumulate all entries in In_raft_coeff, and put them in position
    % according to the area weight they contribute
    MECH.In_raft = accumarray(MECH.S_out_raft(:),MECH.In_raft_coeff(:), size(FSTD.meshRmid(:)),[],0);
    MECH.In_raft = reshape(MECH.In_raft,size(FSTD.meshRmid));
    MECH.In_raft = bsxfun(@times,MECH.In_raft,pi*FSTD.Rmid'.^2);
    
end

%% Ridging interaction

if MECH.ridging
    
    MECH.corr_hmax_ridge = FSTD.H_max_i/FSTD.H_max * MECH.hmax_flag_ridge ... 
    + 1 * (1 - MECH.hmax_flag_ridge); 
    
    MECH.Out_ridge = many_mult( ...
        permute(MECH.diagtwo,[1 3 2 4]), ...
        permute(pi*FSTD.Rmid.^2,[1 2 3 4])', ...
        permute(MECH.Prob_Interact_ridge,[1 3 2 4]), ...
        permute(MECH.gamma_ridge,[3 1 4 2]),...
        permute(FSTD.NumberDist,[1 2 3 4]), ...
        permute(FSTD.NumberDist,[3 4 1 2]) ...
        );
    
    MECH.In_ridge_coeff = many_mult( ...
        permute(MECH.diagone,[1 3 2 4]), ...
        permute(MECH.gamma_ridge,[3 1 4 2]), ...
        permute(MECH.Kfac_ridge,[1 3 2 4]), ...
        permute(MECH.Prob_Interact_ridge,[1 3 2 4]), ... 
        permute(FSTD.NumberDist,[1 2 3 4]), ...
        permute(FSTD.NumberDist,[3 4 1 2]), ...
        permute(MECH.corr_hmax_ridge,[1 2 3 4]) ...
        );
    
    % Accumulate all entries in In_raft_coeff, and put them in position
    % according to the area weight they contribute
    MECH.In_ridge = accumarray(MECH.S_out_ridge(:),MECH.In_ridge_coeff(:), size(FSTD.meshRmid(:)),[],0);
    MECH.In_ridge = reshape(MECH.In_ridge,size(FSTD.meshRmid));
    MECH.In_ridge  = bsxfun(@times,MECH.In_ridge,pi*FSTD.Rmid'.^2);
    
    
end

%%
% Out is a nr by nh by nr by nh matrix, so we just need to sum it over
% the second interaction pairs to get the actual matrix. 
MECH.Out_raft = sum(sum(MECH.Out_raft,3),4);
MECH.Out_ridge = sum(sum(MECH.Out_ridge,3),4);

% In is all set. 

%% Now scale the interaction terms to the external forcing
MECH.In = MECH.In_raft + MECH.In_ridge;
MECH.Out = MECH.Out_raft + MECH.Out_ridge;

if sum(MECH.In(:)) > sum(MECH.Out(:))
    error('Creating Volume, In > Out')
end

MECH.diff = MECH.In - MECH.Out;
MECH.diff_raft = MECH.In_raft - MECH.Out_raft;
MECH.diff_ridge = MECH.In_ridge - MECH.Out_ridge;

% This is an ad-hoc way of doing the normalization, but fine here since
% the Kernel is normalized. If In = Out, use eps to make diff = 0.
sum(MECH.diff(:));
sum(abs(MECH.diff(:)));
diffeps = 0;

if sum(sum(abs(MECH.diff))) == 0
    diffeps = eps;
end


% diff_mech is the ridging mode and is normalized to -1
% It tells how much area must be redistributed.
normalizer = sum(sum(MECH.diff)) + diffeps;

MECH.diff = -MECH.diff / normalizer;

%% Now we consider how the opening due to divergence effects the FSTD
MECH.opening_coll = .5*(MECH.mag - MECH.eps_I);
MECH.opening_div = MECH.mag*MECH.alpha_0;

MECH.In = - MECH.mag*MECH.alpha_c*MECH.In / normalizer;
MECH.Out = - MECH.mag*MECH.alpha_c*MECH.Out / normalizer;

% Here is the "convergent mode" part of the total change in ice
% partial concentrations, which tells how much redistribution is
% done in a "volume conserving" way.

MECH.diff = MECH.mag*MECH.alpha_c*MECH.diff;

% At the moment, volume is not conserved: this is because some
% volume has left the "regular" floe sizes to reach the largest
% thickness category, and some of the area has left the largest
% thickness category, as well. We must update the ice thickness in
% this category to reflect these changes
% On the other hand, area has been correctly reported to all floe
% sizes.
MECH.V_max_in = -integrate_FSTD(MECH.diff(:,1:end-1),FSTD.Hmid(1:end-1),FSTD.dA(:,1:end-1),0);

% If there just isn't sufficient ice, we don't do anything
if sum(FSTD.psi(:)) <= 1e-8
    MECH.diff = 0*FSTD.psi;
    MECH.opening = 0;
    diffeps = eps;
    FSTD.psi = 0*psi;
    FSTD.openwater = 1;
end

%%
% This is the total amount of open water formed by the collisions. Positive
% always.
MECH.opening_coll = MECH.mag*MECH.alpha_c;

% This is the total amount of open water diverged or converged.

% opening_curr = mag*alpha_0;

% This is the net effect, the amount of open water that is formed.
% opening_mech = mag;

% The amount of divergence that isn't accounted for by the collision of
% floes. This leads directly to floes being exported.
MECH.divopening = MECH.eps_I;

%% Loss due to divergence of ice
% Loss in proportion to fractional area

% We need to distinguish between convergence and divergence.

% WHEN THERE IS AN ADVECTIVE PARAMETERIZATION THIS WILL NEED TO CHANGE
% DO NOT FORGET THIS!!!!

if MECH.eps_I < 0
    % When there is convergence, we have no new ice that is added so no
    % advective part
    MECH.convdiv = 0;
else
    % When there is divergence, we lose ice since there is advection out
    MECH.convdiv = 1;
end

if ~ADVECT.DO

MECH.diffadv = (FSTD.psi/(sum(sum(FSTD.psi))+diffeps))*MECH.divopening*MECH.convdiv;

else
    
    MECH.diffadv = 0*FSTD.psi; 

end

% Here is the amount of ice volume which is lost from the largest
% thickness category due to divergenceH
MECH.V_max_out = FSTD.H_max*sum(MECH.diffadv(:,end));

if MECH.V_max_in*OPTS.dt/FSTD.H_max < 1e-8
    MECH.V_max_in = 0;
end

if MECH.V_max_out < eps
    MECH.V_max_out = 0;
end

% Here, now, is the total change in partial concentrations across
% the board, accounting for mechanical combination and divergence
MECH.diff_noadv = MECH.diff;
MECH.diff = MECH.diff - MECH.diffadv;

MECH.opening = -sum(MECH.diff(:));