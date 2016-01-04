%% Mech_Init.m

% This code initializes the mechanical component of the model

%% Flags, parameters, etc.

if ~isfield(MECH,'do_thorndike')
    % Optional code to do the exact thorndike model for thickness
    % interactions
    MECH.do_thorndike = 0;
end



% Whether rafting is turned on or not
if ~isfield(MECH,'rafting')
    MECH.rafting = 1;
end

% Whether ridging is turned on or not
if ~isfield(MECH,'ridging')
    MECH.ridging = 1;
end

if ~isfield(EXFORC,'mech_on')
    % This tells us whether mechanics are on at each timestep
    EXFORC.mech_on = (MECH.rafting || MECH.ridging) * ones(OPTS.nt,1);
    
end

if ~isfield(MECH,'cut_down_interactions')
    MECH.cut_down_interactions = 0;
end

if ~isfield(MECH,'r_raft')
    % Size of a raft
    MECH.r_raft = 10;
end

if ~isfield(MECH,'r_ridge')
    % Size of a ridge
    MECH.r_ridge = 5;
end

if ~isfield(MECH,'k_ridge')
    % Thickness Increase Multiple
    MECH.k_ridge = 5;
    
end

if ~isfield(MECH,'H_raft')
    % Size at which rafting is suppressed
    MECH.H_raft = .3;
end




if ~isfield(EXFORC,'nu')
    % This is the time-series of strain rate invariants
    EXFORC.nu = zeros(OPTS.nt,2);
end

if ~isfield(MECH,'prescribenu')
    % This tells us how to get the ice strain rate.
    MECH.prescribenu = 0;
end



if ~MECH.prescribenu
    
    % Coefficient of drag between ocean and ice
    if ~isfield(OPTS,'ociccoeff')
        
        OPTS.ociccoeff = 1e-2;
        
    end
    % fall-off coefficient as concentration --> 1
    if ~isfield(OPTS,'ocicdelta')
        OPTS.ocicdelta = 10;
    end
    
    if ~isfield(OPTS,'ocicbeta')
        OPTS.ocicbeta = 1;
    end
    
    if ~OCEAN.DO
        error('need ocean to enable the simple ocean strain rate-ice strain rate code')
    end
    
end

if ~isfield(MECH,'dont_guarantee_bigger')
    MECH.dont_guarantee_bigger=0;
end

if ~isfield(MECH,'use_old_interactions')
    MECH.use_old_interactions = 0;
end

%% Calculating interaction matrices

MECH.gamma_raft = calc_gamma_raft_FD(FSTD.Hmid,FSTD.meshHmid,MECH.H_raft);
MECH.gamma_ridge = 1 - MECH.gamma_raft;


% Interaction Matrices for Rafting
if isfield(MECH,'D4_Ridging') && MECH.D4_Ridging == 1
    intstr = sprintf('Mechanics/Interaction_Matrices_4D/nr%drm%dnh%dhm%drra%drri%dkri%d',OPTS.nr,round(max(FSTD.Rmid)),OPTS.nh,round(max(FSTD.H)),MECH.r_raft,MECH.r_ridge,MECH.k_ridge);
    disp('4D_Ridging')
else
    intstr = sprintf('Mechanics/Interaction_Matrices/nr%drm%dnh%dhm%drra%drri%dkri%d',OPTS.nr,round(max(FSTD.Rmid)),OPTS.nh,round(max(FSTD.H)),MECH.r_raft,MECH.r_ridge,MECH.k_ridge);
end

if isfield(MECH,'use_old_interactions') && MECH.use_old_interactions == 1
    
    intstr = [intstr 'old'];
    
else
    
    intstr = [intstr 'new'];
    
end

if isfield(MECH,'dont_guarantee_bigger') && MECH.dont_guarantee_bigger == 1
    intstr = [intstr '_dgb'];
end

intstr = ([OPTS.path_of_code 'Packages/' intstr]);

%%

try
    
    load(intstr)
    fprintf('Was able to find an interaction scheme matching the initial conditions \n')
    
    if MECH.rafting
    MECH.S_R_raft = S_R_raft;
    MECH.S_H_raft = S_H_raft;
    MECH.Kfac_raft = Kfac_raft;
    MECH.Prob_Interact_raft = Prob_Interact_raft;
    end
    
    if MECH.ridging
    MECH.S_R_ridge = S_R_ridge;
    MECH.S_H_ridge = S_H_ridge;
    MECH.Kfac_ridge = Kfac_ridge;
    MECH.Prob_Interact_ridge = Prob_Interact_ridge;
    
    end
    
catch errloading
    
    % Make sure the file didn't exist and it isn't something else
    
    if strcmp(errloading.identifier ,'MATLAB:load:couldNotReadFile')
        
        if ~isfield(MECH,'rafting') || MECH.rafting == 1
            
            MECH.rafting = 1;
            
            disp('Calculating Interaction Matrices for Rafting ...')
            
            if ~(isfield(MECH,'use_old_interactions') && MECH.use_old_interactions == 1)
                disp('Using the new method!')
            end
            
            if isfield(MECH,'dont_guarantee_bigger') && MECH.dont_guarantee_bigger == 1
                disp('Not guaranteeing bigger floes!')
            end
            
            [MECH.S_R_raft,MECH.S_H_raft,MECH.Kfac_raft,MECH.Prob_Interact_raft] = ...
                calc_sizes_raft_FD(FSTD.Rmid,MECH.r_raft,OPTS.domain_width^2,FSTD.Hmid,MECH.dont_guarantee_bigger,MECH.use_old_interactions);
            
            S_R_raft = MECH.S_R_raft;
            S_H_raft = MECH.S_H_raft;
            Kfac_raft = MECH.Kfac_raft;
            Prob_Interact_raft = MECH.Prob_Interact_raft;
            
            save(intstr,'S_R_*','S_H_*','Kfac*','Prob_Interact_*');
        end
        
        
        % Interaction Matrices for Ridging
        if ~isfield(MECH,'ridging') || MECH.ridging == 1
            
            MECH.ridging = 1;
            
            disp('Calculating Interaction Matrices for Ridging ...')
            
            if ~(exist('use_old_interactions','var') && use_old_interactions == 1)
                disp('Using the new method!')
            end
            
            if exist('dont_guarantee_bigger','var') && dont_guarantee_bigger == 1
                disp('Not guaranteeing bigger floes!')
            end
            
            if exist('D4_Ridging','var') && D4_Ridging == 1
                [MECH.S_R_ridge,MECH.S_H_ridge,MECH.Kfac_ridge,MECH.Prob_Interact_ridge] = ...
                    calc_sizes_ridge_FD_2(FSTD.Rmid,OPTS.domain_width^2,FSTD.Hmid,MECH.dont_guarantee_bigger,MECH.use_old_interactions);
                
                S_R_ridge = MECH.S_R_ridge;
                S_H_ridge = MECH.S_H_ridge;
                Kfac_ridge = MECH.Kfac_ridge;
                Prob_Interact_ridge = MECH.Prob_Interact_ridge;
                
                save(intstr,'S_R_*','S_H_*','Kfac*','Prob_Interact_*');
                
            else
                
                [MECH.S_R_ridge,MECH.S_H_ridge,MECH.Kfac_ridge,MECH.Prob_Interact_ridge] = ...
                    calc_sizes_ridge_FD(FSTD.Rmid,MECH.r_ridge,OPTS.domain_width^2,MECH.k_ridge,FSTD.Hmid,MECH.dont_guarantee_bigger,MECH.use_old_interactions);
                
                S_R_ridge = MECH.S_R_ridge;
                S_H_ridge = MECH.S_H_ridge;
                Kfac_ridge = MECH.Kfac_ridge;
                Prob_Interact_ridge = MECH.Prob_Interact_ridge;
                
                save(intstr,'S_R_*','S_H_*','Kfac*','Prob_Interact_*');
                
            end
            
        end
        
    else
        
        % Not sure what this error may be
        throw(errloading)
    end
    
end

% Create diagonal matrices for correctly counting sums
MECH.diagtwo = ones([length(FSTD.Rmid) length(FSTD.Rmid) length(FSTD.Hmid) length(FSTD.Hmid)]);
MECH.diagone = .5*MECH.diagtwo;

for i = 1:length(FSTD.Rmid)
    for j = 1:length(FSTD.H)
        MECH.diagone(i,i,j,j) = 1;
        MECH.diagtwo(i,i,j,j) = 2;
    end
end

% This matrix is nr by nr by nh by nh, and defines an ambiguous index into
% which each of the floes formed by the collision of (r1,h1) and (r2,h2)
% will go. 

% i = MECH.S_out(r1,r2,h1,h2) is a number between 1 and nr +
% nh*(length(FSTD.H)-1)
% The usage is when we are accumulating all of the out data. All of the
% floes formed with size r3 and h3, with index i3, will be accumulated
% together
if MECH.rafting
MECH.S_out_raft = bsxfun(@plus,length(FSTD.Rmid) * (MECH.S_H_raft -1),MECH.S_R_raft);
MECH.S_out_raft = permute(MECH.S_out_raft,[1 3 2 4]); 
MECH.hmax_flag_raft = (MECH.S_H_raft == length(FSTD.Hmid)); 
MECH.hmax_flag_raft = permute(MECH.hmax_flag_raft,[1 3 2 4]); 
MECH.hmax_flag_raft(:,end,:,end) = 0; 
end

if MECH.ridging
MECH.S_out_ridge = bsxfun(@plus,length(FSTD.Rmid) * (MECH.S_H_ridge -1),MECH.S_R_ridge);
MECH.S_out_ridge = permute(MECH.S_out_ridge,[1 3 2 4]); 
MECH.hmax_flag_ridge = (MECH.S_H_ridge == length(FSTD.Hmid)); 
MECH.hmax_flag_ridge = permute(MECH.hmax_flag_ridge,[1 3 2 4]); 
MECH.hmax_flag_ridge(:,end,:,end) = 0; 
end

%% Matrices involved in updating psi. These should be initialized

MECH.In = zeros(length(FSTD.Rmid),length(FSTD.Hmid));
MECH.In_raft = MECH.In;
MECH.In_ridge = MECH.In;
MECH.Out = MECH.In;
MECH.Out_raft = MECH.In;
MECH.Out_ridge = MECH.In;
