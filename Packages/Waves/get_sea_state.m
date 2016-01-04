% get_sea_state
% This function retrieves the wave spectrum at a given time. If it is
% unspecified in EXFORC it is simple the Bretschneider spectrum

if WAVES.prescribe_spec
    % If we want to prescribe the spectrum
    WAVES.spec = EXFORC.wavespec(FSTD.i,:);
else
    % Bretschneider Spectrum
    WAVES.spec = (1/(4*OPTS.g))*(WAVES.H_s^2/WAVES.P_z^2).*(WAVES.Per/WAVES.P_z).^2 .* ... 
    exp((-1/pi) * (WAVES.Per/WAVES.P_z).^4); 

end

