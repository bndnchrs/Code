function [du_out,lossl,lossu] = advect2_upwind(psi,dA,rind,hind,dx,dy,vr,vh,flagrcorr,flaghcorr)
%% advect2_upwind(psi,dA,R,[H H_max],v_r,v_h,dt_sub)
% Here the x,y correspond to the edges of each successive grid cell

% We assume that the incoming function psi is defined so that psi(i,j) is
% the total ice area between rind(i) and rind(i+1) and hind(j) and hind(j+1).
% The imperatives of this routine are the following:

% I: CONSERVE AREA in advection
% II: CONSERVE VOLUME "

% To accomplish I, this simply means that the total area advected out of
% one category must be replaced into the category it enters.

% We would accomplish this, assuming ice is being advected from category
% (i,j) to (i+1,j) by requiring

% dpsi(i,j)*dA(i,j) = dpsi(i+1,j)*dA(i+1,j)

% Accomplishing II is more tricky. The volume in a given category is equal
% to the midpoint thickness multiplied by the concentration. Suppose that a
% piece of ice of area 1 and thickness H(i+1)-eps increases by a thickness 2 * epsilon.
% From the perspective of the model, the mean thickness is now H(i+1)_1/2
% instead of H(i)_1/2, even though the total ice thickness has only changed by a value of 2 epsilon
% This requires some care, and means that the ice area transfer is diminshed
% A bit than might be expected.

% We calculate three terms:

% 1: The area occupied per size: psi(i,j)*dA(i,j)
% 2: The length of the interval across which ice will be exchanged:
% DX = .5 * (dr(i) + dr(j))
% 3: The speed of the advection: vr(i).

% This conserves both volume and area.

% If there is a negative tendency vr(1,:) or vh(:,1), then we have a choice
% if flagrcorr = 1, then this ice that becomes smaller than R(1) will be
% erased from the face of the earth.
% if flaghcorr = 1, then the ice that becomes smaller than H(1) will
% similarly be deleted forever, never to freeze again

if nargin < 7
    flagrcorr = 0;
    flaghcorr = 0;
end

nx = length(rind);
ny = length(hind);

%% We assume that C(i,j) is the total concentration between x_i and x_i+1, and y_j and y_j+1

% Now create a matrix PSI which has 4 extra categories in size and
% thickness
PSI = single(zeros(length(rind) + 4, length(hind) + 4));
DA = PSI; 
VX = PSI;
VY = PSI;

% Do the same for DX
DX = single(zeros(length(rind)+4,1) + dx(end));
% And DY
DY = single(zeros(length(hind)+4,1) + dy(end));

% The first two indices will be the same as the third
DX(3:end-2) = single(dx);
DX(1:2) = DX(3);
DY(3:end-2) = single(dy);
DY(1:2) = DY(3);

% Same here, psi(1:2,1:2) = 0, psi(end-1:end,end-1:end) = 0. 
PSI(3:end-2,3:end-2) = single(psi).*single(dA);

VX(3:end-2,3:end-2) = vr;
VY(3:end-2,3:end-2) = vh;

vxpos = sign(VX);
vypos = sign(VY);

du = zeros(size(PSI));
du_loss = du; 
du_gain = du; 
du_gain_x = du;
du_gain_y = du; 
%% Calculate tendencies
% We need to conserve volume and area

denDX = .5*(DX(1:end-1) + DX(2:end));
denDX = 1./denDX;
denDX = single(denDX); 


denDY = .5*(DY(1:end-1) + DY(2:end));
denDY = 1./denDY;
denDY = single(denDY); 

posupx = vxpos(2:nx+3,2:ny+3); 
posupx = single(posupx == 1); 

posupy = vypos(2:nx+3,2:ny+3); 
posupy = single(posupy == 1); 

%% This ugly and frustrating routing is to ensure that sum(du(:)) = 0 whenever we aren't advecting things out


du_loss_x = abs(VX(2:nx+3,2:ny+3) .* PSI(2:nx+3,2:ny+3) .* ...
    (bsxfun(@times,posupx,denDX(2:nx+3)) + bsxfun(@times,1-posupx,denDX(1:nx+2)))); 

du_loss_y = abs(VY(2:nx+3,2:ny+3) .* PSI(2:nx+3,2:ny+3) .* ...
    (bsxfun(@times,posupy,denDY(2:ny+3)') + bsxfun(@times,1-posupy,denDY(1:ny+2)'))); 

du_max = max(du_loss_x(:) + du_loss_y(:)); 

du_loss(2:nx+3,2:ny+3) = du_loss_x + du_loss_y; 
du_gain_x(3:nx+4,2:ny+3) = du_loss_x .* posupx;
du_gain_x(1:nx+2,2:ny+3) = du_gain_x(1:nx+2,2:ny+3) + du_loss_x .* (1-posupx); 

du_gain_y(2:nx+3,3:ny+4) = du_loss_y .* posupy; 
du_gain_y(2:nx+3,1:ny+2) = du_gain_y(2:nx+3,1:ny+2) + du_loss_y .* (1 - posupy); 

du_gain = du_gain_x + du_gain_y; 

du = du_gain - du_loss; 

if sum(du(:)) ~= 0
 
a = size(du);

[~,b] = max(du(:)); 
resid = sum(du(:)); 
du(b) = du(b) - resid; 

end



[du_out,lossl,lossu] = correctedge(du,flagrcorr,flaghcorr);

du_out_1 = du_out;
lossl = lossl;
lossu = lossu;

% there is the potential that the division by dA is incorrect. 
du_out = du_out ./ dA; 
lossl = lossl ./ dA(1,:); 
lossu = lossu ./ dA(:,1); 
end