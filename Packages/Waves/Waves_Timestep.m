%This code executes one timestep of the swell fracture code

% First we need to generate the spectrum of sea surface heights

% Generate the sea surface height spectrum
% This is the zero-crossing wavelength
WAVES.Lambdaz = OPTS.g* WAVES.P_z^2 / (2*pi);

% Tau is the domain fraction divided by the group velocity of the
% wavelength that is the same as the zero-crossing period.
WAVES.tau = OPTS.Domainwidth / (WAVES.v_group_coeff * (1/2)*(OPTS.g*WAVES.Lambdaz/(2*pi)).^(1/2));

% cts contains the total number of new floes formed in each size class
WAVES.bin_cts = zeros(length(FSTD.Hmid),length(FSTD.Rmid));

WAVES.Omega = 0*FSTD.psi;

% dX is the vector of all of the fracture lengths (i.e. the distance
% between subsequent fracture lengths. It will be added to on each timestep
dX = cell(length(FSTD.H),1);

%% First we create a synthetic domain with vector X.

% The discretization is equal to the mean floe size spacing divided by two
dx = mean(FSTD.dR)/2;
% This creates the vector, which has length equal to Domainwidth
X = 0:dx:OPTS.Domainwidth; % Domain discretization

% Now we get the Fourier representation of the wave field
dlambda = [WAVES.Lambda(1) diff(WAVES.Lambda)]; % dlambda
dlambda(1) = dlambda(2);

% These are the Fourier amplitudes inferred from the wave spectrum

% The coefficients of each fourier term are equal to the attenuation at
% that specific wavelength multiplied by the spectral amplitude.
WAVES.coeffs = bsxfun(@times,sqrt(2*WAVES.spec .* dlambda),exp(1).^(-bsxfun(@times,WAVES.alpha_atten,X))');

% The counter just counts how many realizations of the wave field we have
% done
counter = 0;


% This is the vector of floe sizes we will use in the code. We will then
% interpolate the FSD formed over this field into the actual floe size bins

Rhist = linspace(FSTD.Rint(1),FSTD.Rint(end),length(FSTD.Rmid));
% This is the spacing of the discretization
dRHist = [Rhist(1) diff(Rhist)];
dRHist(1) = dRHist(2);

%% Now we generate the surface wave field and calculate the breaking

% We do this up to a specified value. It could be 1.
while counter < WAVES.maxcounts
    
    counter = counter + 1;
    
    % We create a random phase for each fourier component
    randphase = 2*pi*rand(1,length(WAVES.spec)); % Random phases
    
    % Construct Sea Surface Height Record
    WAVES.fourier = cos(bsxfun(@plus,bsxfun(@times,2*pi./WAVES.Lambda,X'),randphase));
    
    % Sea surface height is obtained by summing the Fourier components
    WAVES.eta = sum(WAVES.coeffs.*WAVES.fourier,2);
    
    % We use the peakfinder code to find peaks of the spectrum
    [maxloc, minloc] = peakfinder2(X,WAVES.eta,WAVES.bandwidth);
    
    % This now is the array where pairs of values are maxes and mins
    extremelocs = sort([minloc maxloc]);
    
    % The values in space of each extreme point, for each H
    extremex = repmat(X(extremelocs),[length(FSTD.H) 1]) ;
    % The SSH values of each extreme point, for each H
    extremeeta = repmat(WAVES.eta(extremelocs)',[length(FSTD.H) 1]) ;
    
    
    %%
    % First increment
    dx1 = (extremex(:,2:end-1) - extremex(:,1:end-2));
    % Second increment
    dx2 = (extremex(:,3:end) - extremex(:,2:end-1));
    
    % Numerator of strain
    D2Y = 2*(extremeeta(:,1:end-2).*dx2 - extremeeta(:,2:end-1).*(dx2 + dx1) + extremeeta(:,3:end).*dx1);
    
    % Denominator of strain
    D2X = dx1 .*dx2 .* (dx1 + dx2);
    
    if isempty(D2X)
        D2X = Inf;
        D2Y = 0;
    end
    
    strain = .5*(bsxfun(@times,FSTD.Hmid',abs(D2Y./D2X)));  % Calculate strain magnitude
    % strain is now a NFRAC by H size array
    
    % Those values exceeding strain threshold
    abovelocs = strain > WAVES.epscrit;
    
    % For each H, we pick out all of the locations that will lead to
    % fracture, and put that in Xlocs. Then we take the distance between
    % them and call that the fracture lengths, and store those in a new
    % array
    if sum(abovelocs(:)) ~= 0
        for jind = 1:size(abovelocs,1)
            
            % The extreme X vales for ice of this thickness
            Xlocs = extremex(jind,abovelocs(jind,:));
            % This is now the vector or fracture lengths implied by the sea
            % surface height field.
            dX{jind} = [dX{jind} diff(Xlocs)];
            
        end
    end
    
end

% To find out the relative likelihood of a floe fracturing, we need to bin
% these counts back onto the original discretization. To do this, we first
% take the histogram of counts onto the even discretization Rhist, which
% becomes the cell array A.

% Then we interpolate onto the actual discretization Rmid by using a kernel
% smoother
%%

for jind = 1:size(abovelocs,1)
    
    
    if sum(dX{jind}) > 1
        
        [A{jind},edges{jind}]  = histcounts(dX{jind},[FSTD.Rint Inf]);
        WAVES.bin_cts(jind,:) = A{jind};
        cens = dX{jind} > max(FSTD.Rmid); 
        smoothcts(jind,:) = ksdensity(dX{jind},FSTD.Rmid);
        if sum(smoothcts(jind,:)) == 0
            smoothcts(jind,:) = A{jind}; 
        end
        
    else
        
        smoothcts(jind,:) = 0*WAVES.bin_cts(jind,:);
        
    end
end

%%
smoothcts = bsxfun(@rdivide,smoothcts,sum(bsxfun(@times,smoothcts+eps,FSTD.dR),2));

%%
if WAVES.smoothing
    
    WAVES.bin_cts = smoothcts;
    
end

%%

Frac = bsxfun(@rdivide,bsxfun(@times,smoothcts,FSTD.Rmid),sum(bsxfun(@times,smoothcts+eps,FSTD.Rmid.*FSTD.dR),2));

WAVES.In = 0*FSTD.psi;
WAVES.Out = 0*FSTD.psi;

if sum(abovelocs(:)) > 1
    
    % Now we determine two things:
    % 1: What is the FSD formed when a floe of size (i,j) fractures? This
    % is called WAVES.F
    for i = 1:length(FSTD.Rmid)
        for j = 1:length(FSTD.H)
            % The FSD is the relative proportion of counts for all sizes
            % smaller than size i, and it is normalized to one.
            WAVES.F(i,j,:) = [Frac(j,1:i) zeros(1,length(FSTD.Rmid)-i)]';
            WAVES.F(i,j,:) = squeeze(WAVES.F(i,j,:))'.*FSTD.dR / sum(eps + squeeze(WAVES.F(i,j,:))'.*FSTD.dR);
            
            % If there aren't any of these, so that we get a NaN, just set
            % it to zero.
            if isnan(sum(squeeze(WAVES.F(i,j,:)).*FSTD.dR'))
                WAVES.F(i,j,:) = zeros(1,length(FSTD.dR));
            end
            
            % Now the fraction of domain covered in floes of size i that
            % will fracture is the total length of fracture lengths smaller
            % than i, divided by the total length of all fractures.
            WAVES.Omega(i,j) = sum(FSTD.Rmid(1:i).*WAVES.bin_cts(jind,1:i)) / sum(eps + FSTD.Rmid.*WAVES.bin_cts(jind,:));
            
        end
    end
    
    %% Now we need to obtain the actual in and out matrix
    
    for r1 = 1:length(FSTD.Rmid)
        % Integrate over floe sizes
        for h1 = 1:length(FSTD.H)
            % The total rate of loss of ice concentration due to fracture
            % is (1/tau) * OMEGA(i,j) * Psi(i,j) * dA(i,j)
            
            % The total rate of gain of concentration into floe size
            % category (k,j) is
            % (1/tau) * OMEGA(i,j) * Psi(i,j) * dA(i,j) * F(i,j,k)
            
            % The time rate of change of psi(k,j), therefore, is this
            % divided by the bin area,
            % (1/tau) * OMEGA(i,j) * Psi(i,j) * dA(i,j) * F(i,j,k) / dA(k,j)
            dpsi = (1/WAVES.tau) * WAVES.Omega(r1,h1) * FSTD.psi(r1,h1)  ;
            
            
            WAVES.In(:,h1) = WAVES.In(:,h1) + dpsi * FSTD.dA(r1,h1) * squeeze(WAVES.F(r1,h1,:)) ./ FSTD.dA(:,h1);
            % And the outgoing is simply Omega/tau * psi
            WAVES.Out(r1,h1) = dpsi;
            
        end
    end
    
    
end

WAVES.diff = WAVES.In - WAVES.Out;

if abs(sum(WAVES.diff.*FSTD.dA)) > eps
    disp('asd')
end

if isnan(WAVES.diff)
    disp('uhoh')
end

% Keeping track of the volume in the largest floe thickness
WAVES.V_max_in = FSTD.H_max*sum(WAVES.In(:,end).*FSTD.dA(:,end));
WAVES.V_max_out = FSTD.H_max*sum(WAVES.Out(:,end).*FSTD.dA(:,end));
