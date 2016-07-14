function dt_temp = calc_max_timestep(FSTD,OPTS,OCEAN)

% Calculate the maximal timestep possible so that the FSTD is never
% smaller than zero anywhere.
dt_temp = compute_max_timestep(FSTD.psi,FSTD.diff,OPTS.dt_sub,0,FSTD.dA,OPTS.debug);
% Now calculate the maximal timestep so that the volume in the highest
% thickness category is >= 0.
dt_temp = compute_max_timestep(FSTD.V_max,FSTD.dV_max,dt_temp,1,FSTD.dA,OPTS.debug);

if OCEAN.DO
    dt_temp = compute_max_timestep(OCEAN.T-OCEAN.Tfrz,OCEAN.dTdt,dt_temp,0,0,OPTS.debug);
    dt_temp = compute_max_timestep(OCEAN.H_ml,OCEAN.w,dt_temp,0,0,OPTS.debug);
end

% If there is an error, dt_temp comes back as a string. We then error
% ourselves out.
if strcmp(OPTS.dt_temp,'dt')
    FSTD.eflag = 1;
    fprintf('Cut timestep is negative at timestep %d',FSTD.i);
end

end

function dt_temp = compute_max_timestep(psi,diff_FD,dt,flag,dA,debug)
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

alpha = min(psitemp(:));

% Now we check to see if any values for psi are less than zero.
if alpha < -eps
    
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
    
end

% This becomes our temporary psi
psitemp = psi+dt_temp*diff_FD;

%%

% However, we can also have an ice concentration that is in excess of 1.
% This needs to be corrected as well!
concnew = integrate_FSTD(psitemp,1,dA,0);
conc = integrate_FSTD(psi,1,dA,0);

% When flag == 1 we are only stopping ice volume from becoming less than 0
% So flag == 0 means we are considering the case when conc > 1
if concnew > 1 + 1e-10 && flag == 0
    
    % This is the excess
    concex = concnew - 1;
    % This is the total difference in the FSTD
    diff_gto = (concnew - conc)/dt_temp;
    %
    %    diff_temp = -(concex./diff_gto) + dt;
    dt_temp_2 = (1 - conc) / diff_gto;
    
    dt_temp = min(dt_temp,dt_temp_2);
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

dt_temp = min(dt_temp,dt);

if dt_temp <= 0
    sprintf('Cut Timestep is negative, equal to %d',dt_temp)
    dt_temp = NaN;
    FSTD.eflag = 1;
end

end