function [S_R,S_H, Kfac, Prob_Interact] = calc_sizes_ridge_FD(R,maxr,A_tot,k_ridge,H,dont_guarantee_bigger,use_old_interactions)
%%
% Calculates the location of the new floe formed by ridging of floes
% of sizes Rt and the likelihood of interaction

H_max = H(end);

% This is the matrix of locations into which a floe size will go
S_R = zeros(length(R),length(R));
Prob_Interact = S_R;
% This is the matrix of locations into which floe thickness will go
S_H = zeros(length(R),length(R),length(H),length(H));
Kfac = 1 + S_H;

%%

% For each floe size, we can calculate the probability that two such floes
% will collide with one another. This is Prob_Interact. 

for r1 = 1:length(R)
    for r2 = 1:length(R)
        %% Interaction Probability
        % This is simply something like A_1 * A_2, but they have been
        % decremented according to the size of a contact zone. 
        
        d = min(R(r1),R(r2));
        s = max(R(r1),R(r2));
        
        if d > maxr
            % If bigger than ridge size, delta is the difference
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
        
        if s > maxr
            acore_s = s^2 - maxr^2;
            acz_s = maxr^2;
        else
            acore_s = 0;
            acz_s = s^2;
        end
        
        % This is now the probability, whre A_tot is used to make these
        % numbers large enough to work with numerically. 
        Prob_Interact(r1,r2) = A_tot^2*(acz_d*acz_s)/(A_tot - acore_s - acore_d)^2;
        
        %% Area Loss
        % This is the new area: it is less than the area of the two
        % colliding floes. 
        
        % Ridge is fixed multiple of initial thickness of small floe.
        Anew = s^2 + d^2 - (1 - 1/k_ridge)*acz_d;
        rnew = sqrt(Anew);
        
        % To figure out where the floe will go, we pick the first new
        % category larger than rnew
        temp = rnew - R;
        % if > 0, then R is smaller. If < 0, send to infinity
        temp(temp < 0) = Inf;
        % The largest R such that R < rnew
        [~, loc] = min(temp);
        % We also guarantee that floe sizes will increase
        % Floe size must increase.
        loc  = max([r1+1,r2+1,loc]);
        loc = min(loc,length(R));
        S_R(r1,r2) = loc;
        
        % We can specify that floes don't need to get bigger if we want
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
            S_R(r1,r2) = loc;
        end
        
        %% Now onto determining outgoing thickness/area loss
        
        % Note: this following sections are not dependent on the type of
        % interaction being considered: it requires only the conservation
        % of volume. 
        
        % Outgoing floe size and floe category are determined using the
        % floe-floe interaction statistics. The outgoing thickness will be
        % determined via conservation of volume.
        
        % These loops concern all interactions between floes that are
        % thinner than the thickest floe class. 
        
        for h1 = 1:length(H)
            for h2 = 1:length(H)
                %% 
                % The important thing to do is to conserve volume. We
                % accomplish this in spite of fixed thickness and area
                % categories by putting "fractional" amounts of ice into
                % each category. This is the role of K_fac.
                
                
                % Initial Volume of the colliding floes
                V_i = H(h1)*pi*R(r1)^2 + H(h2)*pi*R(r2)^2;
                
                % We can estimate what the thickness of the new floe will
                % be by dividing the incident volume by the expected area
                H_an = V_i / (pi*R(loc)^2);
                
                % Now Calculate Thickness Category by doing twho things. 
               
                % First, we find the first thickness category
                % that is larger than the thickness required by moving the
                % floe to its desired size category. 
                temp = H_an - H; 
                temp(temp < 0) = Inf; 
                [~, loc2] = min(temp); 
                
                % Then we find the ean thickness of the incoming floes. It
                % is required that the new thickness be thicker than this
                % mean. Otherwise, the total area would increase, which is
                % a no-no
                
                % The mean thickness of the incoming floes
                hbarin = (V_i / (pi*R(r1)^2 + pi*R(r2)^2));

                % The first thickness category larger than this is volloc
                temp = hbarin - H;
                temp(temp > -eps) = -Inf;
                [~, volloc] = max(temp);
                
                % We send ice to at least the thickness it comes from,
                % requiring that the thickness must increase to be at least
                % as thick as the mean thickness before
                loc2  = max(min(h1,h2),max(loc2,volloc));
                    
                    
                % This is where the thickness will go
                S_H(r1,r2,h1,h2) = loc2;
                
                % Now we can determine the volume of a floe that is formed
                V_out = H(loc2)*pi*R(loc)^2;
                
                % If we are in normal areas, simply take the ratio of the
                % incoming to outgoing volumes to conserve volume
                Kfac(r1,r2,h1,h2) = V_i / V_out;
                
                % In all of these cases, the thickness is increased, so
                % that we lose area automagically. However it is possible
                % that the thickness is not able to be increased, which
                % requires some tact.
                
                % Since there will be a variable thickness category, it can
                % be the case that less area is added than is 
                
            end
        end
        
        %% DEALING WITH THE FINAL THICKNESS CATEGORY
        % We must consider what happens when one of the two interaction
        % partners is a floe belonging to the thickest category. The
        % area will be determined by the probability of interaction and
        % the thickness will be deposited in the largest thickness
        % category
        
%         for h1 = 1:length(H)
%             % If one of the interaction pairs has the thickest thickness
%             % then we must put all of the potential interactions into this
%             % class
%             
%             S_H(r1,r2,h1,end) = length(H);
%             S_H(r1,r2,end,h1) = length(H);
%             
%             % The incident volume is assuming the second floe is at thickness
%             % H_max.
%             V_i = H(h1)*pi*R(r1)^2 + H_max*pi*R(r2)^2;
%             
%             % The area of the formed floe is equal to rnew^2 multiplied by the
%             % thickness of the ice there.
%             V_out = H_max*rnew^2;
%             
%             % Then we need to conserve volume at the thickest category
%             Kfac(r1,r2,h1,end) = V_i / V_out;
%             
%             % We must also do this as if the first floe is in the thickest
%             % category.
%             
%             V_i = H_max*pi*R(r1)^2 + H(h1)*pi*R(r2)^2;
%             V_out = H_max*rnew^2;
%             
%             Kfac(r1,r2,end,h1) = V_i / V_out;
%             
%         end
        
 
        % This factor is the ratio of the new area to the incident area
        Kstar = Kfac(r1,r2,end,end)*R(loc)^2 / (R(r1)^2 + R(r2)^2);
        
        % If this factor is greater than one, we have to reduce it. This is
        % done by assuming that the actual area transfer is smaller, and is
        % equal to the predicted area rnew^2. 
        
        if Kstar > 1
            Kfac(r1,r2,end,end) = (rnew^2/R(loc)^2);
        end
        
        
    end
end

check_ridgemat_out;

end
