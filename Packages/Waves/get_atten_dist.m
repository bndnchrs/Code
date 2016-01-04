% get_atten_dist

% This function calculates the wave attenuation distance W(\lambda) as well
% as the swell fracture timescales \tau(\lambda) using an approximation to
% the results obtained by Kohout and Meylan (2008).


% The expression is ln(alpha) = hbar - (1/6)sqrt(2 \pi \lambda/g) - 3 +
% ln(c/2 x r bar)

WAVES.Lambda = WAVES.Per.^(2)*OPTS.g/(2*pi);

% This is the integral over N x R / integral over N, so the mean floe size
% by number
rbar = integrate_FSTD(FSTD.NumberDist,FSTD.Rmid',FSTD.dA,1)
rbar = integrate_FSTD(FSTD.psi,FSTD.Rmid',FSTD.dA,1)

% This is the integral over psi x H / integral over psi, so the mean ice
% thickness. This is also the same as the volume (psi x H) divided by the
% concentration (psi)
hbar = integrate_FSTD(FSTD.psi,FSTD.Hmid,FSTD.dA,1);

%%

% This will be called if we want to do full interpolation
load('int_interp_coeff')

%logalpha = interp2(interp_H,interp_P,atten,hbar,Per,'linear',-.5)';

if isfield(WAVES,'polyfitter')
    
    if WAVES.polyfitter == 1
        WAVES.logalpha = polyval2(pv1,hbar,WAVES.Per);
    else if WAVES.polyfitter == 2
            WAVES.logalpha = polyval2(pv2,hbar,WAVES.Per);
        else if WAVES.polyfitter == 3
                WAVES.logalpha = polyval2(pv3,hbar,WAVES.Per);
            else if WAVES.polyfitter == 4
                    %% Exact Value from KM
                    WAVES.logalpha = -2.38 + 0*polyval2(pv3,hbar,WAVES.Per);
                end
            end
        end
    end
else
    WAVES.logalpha = polyval2(pv2,hbar,WAVES.Per);
end

if isfield(WAVES,'dobennetts') && WAVES.dobennetts == 1
   
    load bennetts_coeffs
    % Use the wave attenuation model of Bennetts et al (2012). 
    
    F = griddedInterpolant(T',H',logalpha');
    % F(per,thick) now returns ln(alpha)
   
    for kk = 1:length(WAVES.Per)
        
        WAVES.logalpha(kk) = F(WAVES.Per(kk),hbar);
        
    end
    
    clear kk
    
end
 
%% 
% Now we have the attenuation coefficients, which are per floe

alpha_notilde = exp(1).^(WAVES.logalpha);

% The attenuation coefficient as a function of lambda.
WAVES.alpha_atten = alpha_notilde * FSTD.conc/(2*rbar);