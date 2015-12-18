function dt_temp = calc_max_timestep(psi,diff_FD,dt,flag,debug)
%% dt_temp = cut_timestep(psi,diff,dt_sub)
% This function counts the maximum timestep possible (< dt_sub) which
% allows the solutions to remain so that the total ice concentration is
% bounded above by 1 and each individual ice concentration is bounded below
% by 0. 

if nargin < 4
    % If we don't specify flag, set it to zero
    flag = 0; 
    debug = 0; 
end

% Flag allows us to ignore the error that comes when psi is greater than
% one. Useful for re-using this code to deal with thickness. 

% This is the potential value for psi
psitemp = psi+dt*diff_FD;

% The temporary value
dt_temp = dt; 

% If the minimum value of all is less than zero, we need to cut the
% timestep

%%

if min(psitemp(:)) < 0
    
    % if psi is < 0 , we take the maximum timestep possible in order to
    % keep all values >= 0
    
    % All locations that become smaller than zero after the long timestep
    less_than_zeros = psitemp(psitemp < 0);
    % The total difference at each of those points
    diff_ltz = diff_FD(psitemp < 0);
    
    % This is the timestep that would take each value exactly to zero
    dt_temp_ltzs = -(less_than_zeros./diff_ltz) + dt;
    % We have to take the minimum timestep, so that all values after
    % dt_temp will become larger than one
    dt_temp = min(dt_temp_ltzs);
    
    if debug
    
    disp('cutting timestep, too small')
    
    end
    %% If still, we are less than zero, we have to quit
    
    
    % Now do the same thing in case values get larger than 1
else if max(psitemp(:)) > 1 && flag == 0
        % flag is if we are doing thickness calculation: flag will be turned on then as
        % there is no upper bound on ice thickness
        
        % If psi is greater than one, we take the maximum timestep possible
        % in order to keep all values <= 1
        if debug
            disp('cutting timestep, too high')
        end
        
        greater_than_ones = psitemp(psitemp > 1) - 1;
        diff_gto = diff_FD(psitemp > 1);
        dt_temp_gtos = -(greater_than_ones./diff_gto) + dt;
        dt_temp = min(dt_temp_gtos);
        

        
    end
end

psitemp = psi+dt_temp*diff_FD;

% This can become a problem sometimes. We allow machine precision errors in
% the ice advection to take the ice concentration over 1 at times. This
% almost never occurs unless for some special initialization

if sum(psitemp(:)) > 1 + 1e-8 && flag == 0
    
    concex = sum(psitemp(:)) - 1;
    diff_gto = sum(diff_FD(:));
    diff_temp = -(concex./diff_gto) + dt; 
    dt_temp = min(diff_temp); 
   
    if debug
    disp('too much conc, cutting timestep')
    end
    
end


% Now we have calculated a minimum timestep that ensures that all values
% will be greater than zero. If all values would have been greater than
% zero, we ensure that the timestep is small enough so all values are less
% than one.
% 
% if dt_temp <= 0
%     dt_temp = NaN;
% end

if dt_temp <= 0
    sprintf('Cut Timestep is negative, equal to %d',dt_temp)
    dt_temp = NaN;
    FSTD.eflag = 1;
end

end