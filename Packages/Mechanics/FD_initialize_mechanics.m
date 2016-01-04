

if FSTD.DO
    
    %% Strain Rate
    
    if ~isfield(MECH,'eps_I')
        MECH.eps_I = -1e-7;
    end
    
    if ~isfield(MECH,'eps_I_0')
        MECH.eps_I_0 = MECH.eps_I;
    end
    
    if ~isfield(MECH,'eps_II');
        MECH.eps_II = 0;
    end
    
    if ~isfield(MECH,'eps_II_0')
        MECH.eps_II_0 = MECH.eps_II;
    end
    
    if ~isfield(MECH,'do_thorndike')
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
    
    if ~isfield(MECH,'nu')
        % Size at which rafting is suppressed
        MECH.nu = zeros(OPTS.nt,2);
        MECH.nu(:,1) = MECH.eps_I;
        MECH.nu(:,2) = MECH.eps_II;
        
    end
    
    if ~isfield(MECH,'simple_oc_sr')
        % Infer the strain rate tensor simply from the ocean
        MECH.simple_oc_sr = 0;
        
    end
    
    if MECH.simple_oc_sr
        
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
    
    %%
    
    % Ridging and rafting matrices
    MECH.gamma_raft = calc_gamma_raft_FD([FSTD.H FSTD.H_max],FSTD.meshH,MECH.H_raft);
    MECH.gamma_ridge = 1 - MECH.gamma_raft;
    
    
    % Interaction Matrices for Rafting
    if isfield(MECH,'D4_Ridging') && MECH.D4_Ridging == 1
        intstr = sprintf('./Mechanics/Interaction_Matrices_4D/nr%drm%dnh%dhm%drra%drri%dkri%d',OPTS.nr,round(max(FSTD.R)),OPTS.nh,round(max(FSTD.H)),MECH.r_raft,MECH.r_ridge,MECH.k_ridge);
        disp('4D_Ridging')
    else
        intstr = sprintf('./Mechanics/Interaction_Matrices/nr%drm%dnh%dhm%drra%drri%dkri%d',OPTS.nr,round(max(FSTD.R)),OPTS.nh,round(max(FSTD.H)),MECH.r_raft,MECH.r_ridge,MECH.k_ridge);
    end
    
    if isfield(MECH,'use_old_interactions') && MECH.use_old_interactions == 1
        
        intstr = [intstr 'old'];
        
    else
        
        intstr = [intstr 'new'];
        
    end
    
    if isfield(MECH,'dont_guarantee_bigger') && MECH.dont_guarantee_bigger == 1
        intstr = [intstr '_dgb'];
    end
    
    %%
    
    try
        
        load(intstr)
        fprintf('Was able to find an interaction scheme matching the initial conditions \n')
        
        MECH.S_R_raft = S_R_raft;
        MECH.S_H_raft = S_H_raft;
        MECH.Kfac_raft = Kfac_raft;
        MECH.Prob_Interact_raft = Prob_Interact_raft;
        
        MECH.S_R_ridge = S_R_ridge;
        MECH.S_H_ridge = S_H_ridge;
        MECH.Kfac_ridge = Kfac_ridge;
        MECH.Prob_Interact_ridge = Prob_Interact_ridge;
        
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
                    calc_sizes_raft_FD(FSTD.R,MECH.r_raft,OPTS.A_tot,[FSTD.H FSTD.H_max],MECH.dont_guarantee_bigger,MECH.use_old_interactions);
                
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
                        calc_sizes_ridge_FD_2(FSTD.R,OPTS.A_tot,[FSTD.H FSTD.H_max],MECH.dont_guarantee_bigger,MECH.use_old_interactions);
                    
                    S_R_ridge = MECH.S_R_ridge;
                    S_H_ridge = MECH.S_H_ridge;
                    Kfac_ridge = MECH.Kfac_ridge;
                    Prob_Interact_ridge = MECH.Prob_Interact_ridge;
                    
                    save(intstr,'S_R_*','S_H_*','Kfac*','Prob_Interact_*');
                    
                else
                    
                    [MECH.S_R_ridge,MECH.S_H_ridge,MECH.Kfac_ridge,MECH.Prob_Interact_ridge] = ...
                        calc_sizes_ridge_FD(FSTD.R,MECH.r_ridge,OPTS.A_tot,MECH.k_ridge,[FSTD.H FSTD.H_max],MECH.dont_guarantee_bigger,MECH.use_old_interactions);
                    
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
    
    
    % Initial Pressure
    MECH.P_0 = OPTS.H_0*exp(-20*(1-FSTD.conc));
    
    %% Matrices involved in updating psi
    MECH.In = zeros(length(FSTD.R),length(FSTD.H)+1);
    MECH.In_raft = MECH.In;
    MECH.In_ridge = MECH.In;
    MECH.Out = MECH.In;
    MECH.Out_raft = MECH.In;
    MECH.Out_ridge = MECH.In;
    
    
else
    
    % We didn't initialize the model!
    fprintf('NO MAIN MODEL ENABLED... WONT PROCEED FOR MECHANICS \n ')
    
end