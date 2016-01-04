%% get_ice_strain_rate

% This routine pulls or creates the strain rate state at each large-scale
% timestep, and assorted necessasry quantities.

% There are several options:
% 1: MECH.prescribenu == 1
% Here we prescribe the strain rate as a nt by 2 matrix EXFORC.nu(i,j)

% 2: MECH.prescribenu == 0
% Now we decide to calculate the strain rate. There are two sub-options

% 2a: ADVECT.DO && ~MECH.prescribeOCnu
% Now we calculate the strain rate by using the ocean velocity field. This
% means we nead ADVECT on. 

% 2b: MECH.prescribeOCnu

% Here we specify the ocean strain rate invariants, and then using the
% scaling from before turn this into an ice strain rate invariant. This is
% independent of the ocean or advection code

% This is usually the case. We prescribe these fields through the tensor nu
if isfield(MECH,'prescribenu') && MECH.prescribenu
    
    MECH.eps_I = EXFORC.nu(FSTD.i,1);
    MECH.eps_II = EXFORC.nu(FSTD.i,2);
    
else
    
    % On the other hand, if we do not want to do this, we calculate the strain
    % via the simple ocean strain rate below.
    
    % Strain rate is inferred from the ocean strain rate tensor
    
    % oc_to_ic is the coefficient that relates the ocean straining to the
    % ice straining
    MECH.oc_to_ic = OPTS.ociccoeff * (FSTD.conc * OPTS.H_0/FSTD.Hmean) * ( 1 - exp(-(1 + OPTS.ocicbeta-FSTD.conc)/OPTS.ocicdelta));
    
    % If we include advection, then we calculate eps_I using this code
    if ADVECT.DO
        
        % The straining of the ice pack is related to the ocean straining
        % by the oc_to_ic terms. 
        
        MECH.eps_I = OCEAN.StrainInvar(1) * MECH.oc_to_ic;
        MECH.eps_II = OCEAN.StrainInvar(2) * MECH.oc_to_ic;
        
    end
    
    
end


if isfield(MECH,'rescale_eps') && MECH.rescale_eps
    % A seperate rescaling in the ice thickness and concentration
    MECH.eps_I = MECH.eps_I*(FSTD.conc)*(OPTS.H_0/FSTD.Hmean);
    MECH.eps_II = MECH.eps_II*(FSTD.conc)*(OPTS.H_0/FSTD.Hmean);
    
end

%% This is for cases in which we want to turn off and on mech
MECH.eps_I = MECH.eps_I * EXFORC.mech_on(FSTD.i);
MECH.eps_II = MECH.eps_II * EXFORC.mech_on(FSTD.i);

%% Just some simple calculations from the strain rate tensor

% Magnitude of Strain Rate Tensor
MECH.mag = sqrt(MECH.eps_I^2 + MECH.eps_II^2);
% Ratio of Divergence to Strain
costheta = MECH.eps_I/MECH.mag;

% Opening and closing coefficients
MECH.alpha_0 = .5*(1 + costheta);
MECH.alpha_c = .5*(1 - costheta);
