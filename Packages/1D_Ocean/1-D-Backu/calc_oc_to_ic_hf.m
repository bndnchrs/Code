function Q_oi = calc_oc_to_ic_hf(OCEAN)
% If the ocean is on, transmit heat from the ocean to the ice. Q_oi is the
% heat flux exchanged between the ocean and ice, per square meter of ice

% Calculate the density of the sea water
dens = OCEAN.EOS(OCEAN.T,OCEAN.S);
% This is just a prefactor, cp * rho * ustar
prefac = OCEAN.cp_w * dens * OCEAN.ustar_oceice;

% Heat flux from ocean to ice is the exchange calculated at the ice
% base, where the temperature is equal to the freezing temperature.
Q_oi =  prefac * (OCEAN.T - OCEAN.Tfrz);

% We have the option to say that this term is equal to zero always
if OCEAN.no_oi_hf
    Q_oi = 0;
end


end