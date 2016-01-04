%% Initializing the advection scheme, used to transport ice in from other locations
% We need:
% To know whether we will be advecting things in
% Also need the strain invariants that will prognose the convergence


if ~isfield(OCEAN,'eps_II')
    % Ocean Velocities
    OCEAN.eps_II = zeros(OPTS.nt,1); 
    
end

if ~isfield(OCEAN,'UVEL')
    OCEAN.UVEL =  zeros(2,OPTS.nt); 
end

if ~isfield(ADVECT,'FSTD_in')
    % This is the FSTD of the ice on the western boundary
    ADVECT.FSTD_in = 0 * FSTD.meshR + 1/numel(FSTD.meshR);
end

if ~isfield(OPTS,'ociccoeff')
    % This is the coefficient which translates ocean straining to ice
    % straining
    OPTS.ociccoeff = 1e-2;
    
end

if ~isfield(OPTS,'ocicdelta')
    % Falloff coefficient of strain rate dependance on concentration
    if MECH.DO && MECH.simple_oc_sr
        % We wont need it
        % OPTS.ocicdelta = % MECH.ocicdelta;
    else
        OPTS.ocicdelta = 10;
    end
end

% This is a cutoff of concentration for the advection code. 
if ~isfield(OPTS,'ocicbeta')
     if MECH.DO && MECH.simple_oc_sr
        % We wont need it
     else
        OPTS.ocicbeta = 1;
     end
end

if ~isfield(OPTS,'Domainwidth')
    OPTS.Domainwidth = 1e4; 
end

