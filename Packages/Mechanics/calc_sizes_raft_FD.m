function [S_R,S_H, Kfac, Prob_Interact] = calc_sizes_raft_FD(R,maxr,A_tot,H,dont_guarantee_bigger,use_old_interactions)
%%
% Calculates the location of the new floe formed by rafting of floes
% of sizes Rt and the likelihood of interaction

Hstar = H;
H_max = H(end);
H = H(1:end-1);

S_R = zeros(length(R),length(R));
Prob_Interact = S_R;
S_H = zeros(length(R),length(R),length(H)+1,length(H)+1);
Kfac = 1 + S_H;

%%

for i = 1:length(R)
    for j = 1:length(R)
        %% Interaction Probability
        

        
        d = min(R(i),R(j));
        s = max(R(i),R(j));
        
        if d > maxr
            % If bigger than raft size, delta is the difference
            delta_d = maxr;
        else
            % Otherwise width is the whole thing
            delta_d = d;
        end
        
        % Core is the non-delta part
        acore_d = (d - delta_d)^2;
        % Contact zone is the rest
        acz_d = d^2 - acore_d;
        
        % Bigger floe has same radius contact zone
        acore_s = (s - delta_d)^2;
        acz_s = s^2 - acore_s;
        
        
        num = acz_d*acz_s;
        denom = (A_tot - acore_s - acore_d)^2;
        Prob_Interact(i,j) = A_tot^2*num/denom;
        

        
        
        
        
        %% Area Loss
        
        % Radius of the newly formed floe
        if d^2 < maxr^2/2
            rnew = sqrt(s^2 + d^2 -    d^2/2);
        else
            rnew = sqrt(s^2 + d^2 - maxr^2/2);
        end
        %%
        if exist('dont_guarantee_bigger','var') && dont_guarantee_bigger
            % Picking the new floe location
            % The first floe category thicker than that
            temp = rnew - R;
            % if > 0, then R is smaller. If < 0, send to infinity
            temp(temp < 0) = Inf;
            % The largest R such that R < rnew
            [~, loc] = min(temp);
            % No larger than the length... actually impossible. 
            loc = min(loc,length(R));
            S_R(i,j) = loc;
            
        else

           % The first floe category thicker than that
            temp = rnew - R;
            % if > 0, then R is smaller. If < 0, send to infinity
            temp(temp < 0) = Inf;
            % The largest R such that R < rnew
            [~, loc] = min(temp);
            
            % Floe size must increase. 
            loc  = max([i+1,j+1,loc]);
            loc = min(loc,length(R));
            S_R(i,j) = loc;
            
        end
        %% Now onto determining outgoing thickness/area loss
        
        % Note: this following sections are not dependent on the type of
        % interaction being considered: it requires only the conservation
        % of volume.
        
        % Outgoing floe size and floe category are determined using the
        % floe-floe interaction statistics. The outgoing thickness will be
        % determined via conservation of volume.
        
        for k = 1:length(H)
            for l = 1:length(H)
                %%
                % The important thing to do is to conserve volume. We
                % accomplish this in spite of fixed thickness and area
                % categories by putting "fractional" amounts of ice into
                % each category. This is the role of K_fac.
                
                
                % Initial Volume of the colliding floes
                V_i = H(k)*pi*R(i)^2 + H(l)*pi*R(j)^2;
                
                % The predicted thickness given the fixed area category
                H_an = V_i / (pi*R(loc)^2);
                
               % Now Calculate Thickness Category
                temp = H_an - Hstar; 
                temp(temp < 0) = Inf; 
                [~, loc2] = min(temp); 
                
                % Thickness may increase up to our maximum thickness
                % category (k+1,l+1 ) < loc2 < length(H_max)
                
                if exist('use_old_interactions','var') && use_old_interactions == 1
                    % If we want everything to go up in size and thickness
                    loc2 = max([k+1 l+1 loc2]);
                    
                else
                    % If we don't
                    % The mean thickness of the incoming floes
                    hbarin = (V_i / (pi*R(i)^2 + pi*R(j)^2));                
                    % The first thickness category thicker than that
                    temp = hbarin - Hstar; 
                    temp(temp > -eps) = -Inf; 
                    [~, volloc] = max(temp);
                    
                    % We send ice to at least the size it comes from, and
                    % the area must not increase at the new size/thickness
                    loc2  = max(min(k,l),max(loc2,volloc));
                end
                % This is where the thickness will go
                
                S_H(i,j,k,l) = loc2;
                
                % Here is the thickness which is generated in that category
                V_out = Hstar(loc2)*pi*R(loc)^2;
                
                % If we are in normal areas, simply take the ratio of the
                % incoming to outgoing volumes to conserve volume
                Kfac(i,j,k,l) = V_i / V_out;
                
                % In all of these cases, the thickness is increased, so
                % that we lose area automagically. However it is possible
                % that the thickness is not able to be increased, which
                % requires some tact.
                
            end
            
            % We must consider what happens when one of the two interaction
            % partners is a floe belonging to the thickest category. The
            % area will be determined by the probability of interaction and
            % the thickness will be deposited in the largest thickness
            % category
            
            S_H(i,j,k,end) = length(H) + 1;
            S_H(i,j,end,k) = length(H) + 1;
            
            % First floe is regular
            V_i = H(k)*pi*R(i)^2 + H_max*pi*R(j)^2;
            
            % Outgoing floe has a reduced area rnew
            V_out = H_max*rnew^2;
            
            Kfac(i,j,k,end) = V_i / V_out;
            
            V_i = H_max*pi*R(i)^2 + H(k)*pi*R(j)^2;
            V_out = H_max*R(loc)^2;
            
            Kfac(i,j,end,k) = V_i / V_out;
            
            
            
        end
        
        % Last case is if both floes are the largest thickness, trivially
        % they put ice in the thickest category
        S_H(i,j,end,end) = length(H)+1;
        
        % The ratio of volumes is equal to just the ratio of incident areas
        % since the thickness doesn't change and so we are in the FSD
        % setting, where we will later  bump up the thickness to account
        %for the decreased area but not decreased volume
        
        if loc == length(R)
            Kfac(i,j,end,end) = (rnew^2/R(loc)^2);
        end
        
        Kstar = Kfac(i,j,end,end)*R(loc)^2 / (R(i)^2 + R(j)^2);
        if Kstar > 1
            Kfac(i,j,end,end) = (rnew^2/R(loc)^2);
        end
        
    end
end

check_ridgemat_out;

end
