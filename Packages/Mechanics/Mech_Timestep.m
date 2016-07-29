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
    
    MECH.corr_hmax_raft = (FSTD.H_max_i / FSTD.H_max) * MECH.hmax_flag_raft ...
        + 1 * (1 - MECH.hmax_flag_raft);
    
    % Have now removed this correction.
    
    
    MECH.corr_hmax_raft = 0* MECH.corr_hmax_raft + 1;
    
    
    MECH.Out_raft = many_mult( ...
        permute(MECH.diagtwo,[1 3 2 4]), ...
        permute(pi*FSTD.Rmid.^2,[1 2 3 4])', ...
        permute(MECH.Prob_Interact_raft,[1 3 2 4]), ...
        permute(MECH.gamma_raft,[3 1 4 2]),...
        permute(FSTD.NumberDist,[1 2 3 4]), ...
        permute(FSTD.NumberDist,[3 4 1 2]), ...
        permute(FSTD.dA,[1 2 3 4]), ...
        permute(FSTD.dA,[3 4 1 2]) ...
        );
    
    MECH.In_raft_coeff = many_mult( ...
        permute(MECH.diagone,[1 3 2 4]), ...
        permute(MECH.gamma_raft,[3 1 4 2]), ...
        permute(MECH.Kfac_raft,[1 3 2 4]), ...
        permute(MECH.Prob_Interact_raft,[1 3 2 4]), ...
        permute(FSTD.NumberDist,[1 2 3 4]), ...
        permute(FSTD.NumberDist,[3 4 1 2]), ...
        permute(MECH.corr_hmax_raft,[1 2 3 4]), ...
        permute(FSTD.dA,[1 2 3 4]), ...
        permute(FSTD.dA,[3 4 1 2]) ...
        );
    
    % Accumulate all entries in In_raft_coeff, and put them in position
    % according to the area weight they contribute
    MECH.In_raft = accumarray(MECH.S_out_raft(:),MECH.In_raft_coeff(:), size(FSTD.meshRmid(:)),[],0);
    MECH.In_raft = reshape(MECH.In_raft,size(FSTD.meshRmid));
    MECH.In_raft = bsxfun(@times,MECH.In_raft,pi*FSTD.Rmid'.^2);
    
end

%% Ridging interaction

if MECH.ridging
    
    MECH.corr_hmax_ridge = (FSTD.H_max_i / FSTD.H_max) * MECH.hmax_flag_ridge ...
        + 1 * (1 - MECH.hmax_flag_ridge);
    
    % Have now removed this correction.
    MECH.corr_hmax_ridge = 0* MECH.corr_hmax_ridge + 1;
    
    
    MECH.Out_ridge = many_mult( ...
        permute(MECH.diagtwo,[1 3 2 4]), ...
        permute(pi*FSTD.Rmid.^2,[1 2 3 4])', ...
        permute(MECH.Prob_Interact_ridge,[1 3 2 4]), ...
        permute(MECH.gamma_ridge,[3 1 4 2]),...
        permute(FSTD.NumberDist,[1 2 3 4]), ...
        permute(FSTD.NumberDist,[3 4 1 2]), ...
        permute(FSTD.dA,[1 2 3 4]), ...
        permute(FSTD.dA,[3 4 1 2]) ...
        );
    
    MECH.In_ridge_coeff = many_mult( ...
        permute(MECH.diagone,[1 3 2 4]), ...
        permute(MECH.gamma_ridge,[3 1 4 2]), ...
        permute(MECH.Kfac_ridge,[1 3 2 4]), ...
        permute(MECH.Prob_Interact_ridge,[1 3 2 4]), ...
        permute(FSTD.NumberDist,[1 2 3 4]), ...
        permute(FSTD.NumberDist,[3 4 1 2]), ...
        permute(MECH.corr_hmax_ridge,[1 2 3 4]), ...
        permute(FSTD.dA,[1 2 3 4]), ...
        permute(FSTD.dA,[3 4 1 2]) ...
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

% Just remember that MECH.diff is df/dt * dr * dh so it has units of area
% and is not the MECH.diff that we will use to update the FSD later
MECH.diff = MECH.In - MECH.Out;

if sum(MECH.Out_ridge(:)) < sum(MECH.In_ridge(:))
    disp('Creating Volume, In > Out')
    FSTD.i
end

MECH.diff_raft = MECH.In_raft - MECH.Out_raft;
MECH.diff_ridge = MECH.In_ridge - MECH.Out_ridge;


%% Next we want to control how this redistribution affects the thickest
% class of floes.

% First we figure out what the thickness of the ice in that category should
% be given the interactions occuring only between floes in that thickest
% category

% This is the change in volume in that category because of changes in area
MECH.V_max_change = sum(MECH.Out(:,end).*FSTD.meshHmid(:,end));
% This is the area of floes that are deformed in that category
MECH.A_max_change_out = sum(MECH.Out(:,end));
% This is the area of floes that are now formed
MECH.A_max_change_in = sum(MECH.In(:,end));
% This is the "thickness" of those new floes.
MECH.H_max_change = FSTD.Hmid(end) * (MECH.A_max_change_out / MECH.A_max_change_in);

% This is the amount of volume that becomes these floes due to the
% interactions between thinner floes.
MECH.A_max_in = integrate_FSTD(MECH.diff(:,1:end-1),1,1,0);
MECH.V_max_coll = integrate_FSTD(MECH.diff(:,1:end-1),FSTD.meshHmid(:,1:end-1),1,0);

H_max_new = (MECH.V_max_change + MECH.V_max_coll)/(MECH.A_max_in + MECH.A_max_change_in);

%% MECH.diff is the time rate of change of f(r,h) dr dh. The sum over
% MECH.diff is now the total time rate of change of ice concentration.
% It has to be normalized to be -1 at this point, after which we will
% switch back to a tendency for f(r,h).

% This is an ad-hoc way of doing the normalization, but fine here since
% the Kernel is normalized. If In = Out, use eps to make diff = 0.
diffeps = 0;

if sum(sum(abs(MECH.diff))) == 0
    diffeps = eps;
end

% diff_mech is the ridging mode and much be normalized to -1
% It tells how much area must be redistributed.
% Remember MECH.diff right now is df/dt * dr * dh
normalizer = integrate_FSTD(MECH.diff,1,1,0)+ diffeps;

% Now we divide by dA, so that sum(MECH.diff.*FSTD.dA) = -1;
MECH.diff = -MECH.diff ./ (normalizer * FSTD.dA);
MECH.In = -MECH.mag*MECH.alpha_c*MECH.In ./ (normalizer*FSTD.dA);
MECH.Out = -MECH.mag*MECH.alpha_c*MECH.Out ./ (normalizer * FSTD.dA);

%% Now we consider how the opening due to divergence effects the FSTD
MECH.opening_coll = .5*(MECH.mag - MECH.eps_I);
MECH.opening_div = MECH.mag*MECH.alpha_0;

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

% This is the amount of ice volume moved to the thickest ice class by
% collisions between floes outside of that class

MECH.V_max_in = -integrate_FSTD(MECH.diff(:,1:end-1),FSTD.meshHmid(:,1:end-1),FSTD.dA(:,1:end-1),0);
MECH.A_max_in = sum(MECH.diff(:,end).*FSTD.dA(:,end));
MECH.H_max_in = MECH.V_max_in / MECH.A_max_in;

% This is ice with thickness FSTD.H_mid(end)

% Now we deal with what happens inside that category
% This is the total volume of ice that is fluxed out of this category of
% ice thickness (into the same category)



% if MECH.A_max_in > 0
%
%     MECH.H_max_in = MECH.V_max_in / MECH.A_max_in;
%
% end

% If there just isn't sufficient ice, we don't do anything
if sum(FSTD.psi(:).*FSTD.dA(:)) <= 1e-8
    MECH.diff = 0*FSTD.psi;
    MECH.opening = 0;
    diffeps = eps;
    FSTD.psi = 0*FSTD.psi;
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
    
    MECH.diffadv = (FSTD.psi/(sum(sum(FSTD.psi))+diffeps))*MECH.divopening*MECH.convdiv.*FSTD.dA;
    
else
    
    MECH.diffadv = 0*FSTD.psi;
    
end

% Here is the amount of ice volume which is lost from the largest
% thickness category due to divergenceH
MECH.V_max_out = FSTD.H_max*sum(MECH.diffadv(:,end).*FSTD.dA(:,end));

% if MECH.V_max_in*OPTS.dt/FSTD.H_max < 1e-8
%     MECH.V_max_in = 0;
%     MECH.A_max_in = 0;
%     MECH.diff(:,end) = 0;
% end

if MECH.V_max_out < eps
    MECH.V_max_out = 0;
end

% Here, now, is the total change in partial concentrations across
% the board, accounting for mechanical combination and divergence
MECH.diff_noadv = MECH.diff;
MECH.diff = MECH.diff - MECH.diffadv;

MECH.opening = integrate_FSTD(MECH.diff,1,FSTD.dA,0);