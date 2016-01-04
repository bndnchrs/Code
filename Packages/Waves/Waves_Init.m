% Swell_Init
% This routine initializes the swell fracture code

OPTS.g = 9.81; 

if ~isfield(OPTS,'Domainwidth')
    % Width of Domain
    OPTS.Domainwidth = 1e4; 
end

if ~isfield(WAVES,'Per')

    % Unless pre-defined, we choose the exact period bins so that each
    % wavelength will fracture into one and only one radius bin
    WAVES.Per = (4*pi*FSTD.Rmid/OPTS.g).^(1/2);
     
end

WAVES.Lambda = WAVES.Per.^(2)*OPTS.g/(2*pi);

if ~isfield(WAVES,'prescribe_spec')
    % Whether the wave spectrum is imposed or not. If  it is not imposed we
    % set it to be equal to the bretschneider spectrum. 
    WAVES.prescribe_spec = 0; 
end

if WAVES.prescribe_spec == 0
    
    if ~isfield(WAVES,'H_s')
        WAVES.H_s = 2; % Significant Wave Height
    end
    
    if ~isfield(WAVES,'P_z')
        WAVES.P_z = 6; % Zero-crossing period
    end 
    
end

if ~isfield(EXFORC,'stormy')
    % Flag for whether we are allowing fracture to occur
    EXFORC.stormy = ones(1,OPTS.nt); 
end

% If we want to interpolate the attenuation coefficients, use this. 
if isfield(WAVES,'Do_interp_atten') && WAVES.Do_interp_atten 
    
    load('Swell/interp_coeff')
    
end

if ~isfield(WAVES,'epscrit')

    % The critical strain rate
WAVES.epscrit = 3e-5;

end

if ~isfield(WAVES,'v_group_coeff')
    % Fraction of group velocity that waves move at
    WAVES.v_group_coeff = 1; 
end

if ~isfield(WAVES,'bandwidth')
    % The minimum distance between two successive peaks
    WAVES.bandwidth = 10; 
end

if ~isfield(WAVES,'maxcounts')
   % The maximum number of surface wave field realizations
    WAVES.maxcounts = 5; 
end

if ~isfield(WAVES,'smoothing')
    WAVES.smoothing = 1; 
end

%% Matrices for Handling In and Out
WAVES.In = 0*FSTD.psi;
WAVES.Out = 0*FSTD.psi;
