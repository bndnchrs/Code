if FSTD.DO
    
    if ~isfield(THERMO,'deltamerge')
        % Floes of size r will merge and become floes of size between r and
        % r x deltamerge
        THERMO.deltamerge = 5; 
    end
      
    
end