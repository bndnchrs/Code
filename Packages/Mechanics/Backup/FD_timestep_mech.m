%% Function FD_timestep_mech
% This routing calculates the tendency at each floe size and thickness
% according to the FD parameterizations, and also updates the
% large-ice-thickness class appropriately.

% There are two parts to the tendency due to mechanical forcing
% First, there is the outgoing term. This is the total fraction at each
% size that goes to other floe sizs. 

% Second, there is the incoming term. This is the determination of how the
% first term is distributed amongst larger floes. 

%% Calculate outgoing term

% K_raft(r1,r2,h1,h2) is the probability that the floe of size (r1,r2) will
% interact with a floe of size (h1,h2); 
if MECH.rafting
K_raft =  MECH.Prob_Interact_raft;
end

if MECH.ridging
MECH.Prob_Interact_ridge = MECH.Prob_Interact_ridge;
end

% For each combination of two floe sizes, we figure out which floe size
% will be formed when they combine. 
for r1 = 1:length(FSTD.Rmid)
    % First floe thickness
    for r2 = 1:length(FSTD.Rmid)
        % Second floe size
        
        % This is the size of a new floe formed by r1 and r2 in rafting
        if MECH.rafting
        raft_loc_r = MECH.S_R_raft(r1,r2);
        end
        % This is the size of a new floe formed by r1 and r2 in ridging        
     
        if MECH.ridging
        ridge_loc_r = MECH.S_R_ridge(r1,r2);
        end
        % Now we loop over each combination of thicknesses
        
        % Now for each thickness 
        for h1 = 1:length(FSTD.Hmid)
            % Second floe thickness
            for h2 = 1:length(FSTD.Hmid)
        
                
                % This is the thickness category of a new floe formed by
                % (r1,h1) and (r2,h2) when they ridge
                if MECH.ridging
                    ridge_loc_h = MECH.S_H_ridge(r1,r2,h1,h2);
                    % If the location of the new floe is in the thickest floe
                    % category, there is a correction that needs to be made to
                    % conserve volume, since this thickness can change.
                    if ridge_loc_h == length(FSTD.Hmid)
                    
                        if FSTD.H_max ~= 0
                            correct_Hmax_ridge = FSTD.H_max_i/FSTD.H_max;
                        else
                            correct_Hmax_ridge = 1;
                        end
                        
                    else
                        correct_Hmax_ridge = 1;
                    end
                    
                    
                end
                
                if MECH.rafting
                    raft_loc_h = MECH.S_H_raft(r1,r2,h1,h2);
                    
                    % Same thing for rafting
                    if raft_loc_h == length(FSTD.Hmid)
                        if FSTD.H_max ~= 0
                            correct_Hmax_raft = FSTD.H_max_i/FSTD.H_max;
                        else
                            correct_Hmax_raft = 1;
                        end
                    else
                        correct_Hmax_raft = 1;
                    end
                end
                
                
                
                %% Rafting Step
                if MECH.rafting % If we are doing rafting
                    
                    MECH.Out_raft(r1,h1) =  MECH.Out_raft(r1,h1) + ... % Subtracting from this category
                        MECH.diagtwo(r1,r2,h1,h2)*K_raft(r1,r2)*MECH.gamma_raft(h1,h2)*FSTD.NumberDist(r1,h1)*FSTD.NumberDist(r2,h2)*pi*FSTD.Rmid(r1)^2;
                    
                    % Term that gets added to higher sizes
                    MECH.In_raft(raft_loc_r,raft_loc_h) = MECH.In_raft(raft_loc_r,raft_loc_h) + ... % Adding in
                        MECH.gamma_raft(h1,h2)* ... % Likelihood of rafting
                        MECH.diagone(r1,r2,h1,h2)*correct_Hmax_raft* ... % Correction if going to thickest category and on diagonal
                        MECH.Kfac_raft(r1,r2,h1,h2)* ... % Correction for volume mis-match between incoming and outgoing
                        pi*FSTD.Rmid(raft_loc_r)^2 * ... % Area of the new floe
                        K_raft(r1,r2) * ... % Probability of pair interaction
                        FSTD.NumberDist(r1,h1)*FSTD.NumberDist(r2,h2); % Number distribution collision probability
                     
                    % Term that gets lost from this floe size
                    
                
                end
                
                %% Ridging Step
                if MECH.ridging
                    
                    % In from ridging combination of (r1,h1) and (r2,h2)
                    % We also have a probability of collision that depends
                    % on h_1 and h_2 now
                    MECH.In_ridge(ridge_loc_r,ridge_loc_h) = MECH.In_ridge(ridge_loc_r,ridge_loc_h) + ...
                        + MECH.gamma_ridge(h1,h2)* ...
                        MECH.diagone(r1,r2,h1,h2)*correct_Hmax_ridge* ...
                        MECH.Kfac_ridge(r1,r2,h1,h2)* ...
                        pi*FSTD.Rmid(ridge_loc_r)^2* ...
                        MECH.Prob_Interact_ridge(r1,r2)* ...
                        FSTD.NumberDist(r1,h1)*FSTD.NumberDist(r2,h2);
                    
                    % Out from ridging combination with (r2,h2)
                    MECH.Out_ridge(r1,h1) = MECH.Out_ridge(r1,h1) + ...
                        MECH.diagtwo(r1,r2,h1,h2)*MECH.Prob_Interact_ridge(r1,r2)*MECH.gamma_ridge(h1,h2)*FSTD.NumberDist(r1,h1)*FSTD.NumberDist(r2,h2)*pi*FSTD.Rmid(r1)^2;
                end
                
                % Test to see if it is broken
                
                
            end
            
        end
    end
end

%                if sum(Out_ridge(:)) - sum(In_ridge(:)) < 0
%                    error('Ridge Broken')
%                else
%                    if sum(Out_raft(:)) - sum(In_raft(:)) < 0
%                        error('Raft Broken')
%                    end
%                end

%% Here we handle what happens to the thickest floe class

% We now just treat the top row of floe thicknesses as its
% own FSD which adheres to the usual FSD equation. These
% thick floes will just stay in the thick category
for r1 = 1:length(FSTD.Rmid)
    for r2 = 1:length(FSTD.Rmid)
        
        if MECH.ridging
            
            ridge_loc_r = MECH.S_R_ridge(r1,r2);
            
        end
        
        if MECH.rafting
            
            raft_loc_r = MECH.S_R_raft(r1,r2);
            
        end
        diagone = .5;
        diagtwo = 1;
        
        
        
        if r1 == r2
            diagone = 1;
            diagtwo = 2;
        end
        
        %% Rafting Step
        if MECH.rafting
            
            % In from rafting combination of (r1,h_max) and (r2,h_max)
            MECH.In_raft(raft_loc_r,end) = MECH.In_raft(raft_loc_r,end) + ...
                + diagone*K_raft(r1,r2)*MECH.gamma_raft(end,end)*FSTD.NumberDist(r1,end)*FSTD.NumberDist(r2,end) * ...
                pi * FSTD.Rmid(raft_loc_r)^2 * MECH.Kfac_raft(r1,r2,end,end);
            
            % Out from rafting combination of (r1,h_max) and (r2,h_max)
            MECH.Out_raft(r1,end) = MECH.Out_raft(r1,end) + ...
                diagtwo * K_raft(r1,r2)*MECH.gamma_raft(end,end) * FSTD.NumberDist(r1,end) * FSTD.NumberDist(r2,end) * ...
                pi * FSTD.Rmid(r1)^2 ;
            
        end
        
        %% Ridging Step
        if MECH.ridging
            
            % In from ridging combination of (r1,h_max) and (r2,h_max)
            
            
            MECH.In_ridge(ridge_loc_r,end) = MECH.In_ridge(ridge_loc_r,end) + ...
                + diagone*MECH.Prob_Interact_ridge(r1,r2)*MECH.gamma_ridge(end,end)*FSTD.NumberDist(r1,end)*FSTD.NumberDist(r2,end) * ...
                pi * FSTD.Rmid(ridge_loc_r)^2 * MECH.Kfac_ridge(r1,r2,end,end);
            
            % Out from ridging combination of (r1,h_max) and (r2,h_max)
            MECH.Out_ridge(r1,end) = MECH.Out_ridge(r1,end) + ...
                diagtwo * MECH.Prob_Interact_ridge(r1,r2)*MECH.gamma_ridge(end,end) * FSTD.NumberDist(r1,end) * FSTD.NumberDist(r2,end) * ...
                pi * FSTD.Rmid(r1)^2 ;
        end
        
    end
end


%%
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
%%

% diff_mech is the ridging mode and is normalized to -1
% It tells how much area must be redistributed.
normalizer = sum(sum(MECH.diff)) + diffeps;

MECH.diff = -MECH.diff / normalizer;


%%
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