function [edgnorm_u,edgnorm_z] = goafem_stochcol_L2edge(xy, sols_u_loc, sols_z_loc, coords, evt)
% GOAFEM_STOCHCOL_L2EDGE computes local norms for both primal and
% dual problems
%
% input: 
%              xy     mesh coordinates
%      sols_u_loc     pirmal solution at collocation points
%      sols_z_loc     dual solution at colocation points
%          coords     collocation points
%             evt     elements to vertex mapping
%
% output:
%       edgnorm_u     array of edge norm for primal solution 
%       edgnorm_z     array of edge norm for dual solution 
%
% TR; 28 September 2022
% Initialisation
ncoords = size(coords,1);
nel = size(evt,1);

edgnorm_u = zeros(nel,size(coords,1));
edgnorm_z = zeros(nel,size(coords,1));

x = xy(:,1);
y = xy(:,2);

% Recover local coordinates
xl_v = zeros(nel,3);
yl_v = zeros(nel,3);
for i = 1:3
    xl_v(:,i) = x(evt(:,i));
    yl_v(:,i) = y(evt(:,i));
end

% Compute Edge Lengths
hx_v = xl_v(:,3) - xl_v(:,2);
hy_v = yl_v(:,3) - yl_v(:,2);
he_v(:,1) = sqrt(hx_v.*hx_v + hy_v.*hy_v);

hx_v = xl_v(:,1) - xl_v(:,3);
hy_v = yl_v(:,1) - yl_v(:,3);
he_v(:,2) = sqrt(hx_v.*hx_v + hy_v.*hy_v);

hx_v = xl_v(:,2) - xl_v(:,1);
hy_v = yl_v(:,2) - yl_v(:,1);
he_v(:,3) = sqrt(hx_v.*hx_v + hy_v.*hy_v);

% Construct the integration rule
nngpt = 7;
[oneg,onew] = gausspoints_oned(nngpt);

for k = 1:ncoords % loop over collocation points
    u_gal_loc = sols_u_loc{k};
    z_gal_loc = sols_z_loc{k};
    u_loc_norm = zeros(nel,1);
    z_loc_norm = zeros(nel,1);
    % Loop over gauss points
    for igpt = 1:nngpt
        sigpt = oneg(igpt);
        wght = onew(igpt);
        s = (1.0 + sigpt)/2.0;
        smat = [s, 1-s; 0, s; s, 0];
        urt = zeros(nel,3);
        zrt = zeros(nel,3);
        for i = 1:3
            [jac_v,~,phi_v,~,~] = tderiv(smat(i,1), smat(i,2), xl_v, yl_v);
            for ivtx = 1:3
                urt = urt + u_gal_loc(:,ivtx).*phi_v;
                zrt = zrt + z_gal_loc(:,ivtx).*phi_v;
            end
        end
        for i = 1:3
        	u_loc_norm = u_loc_norm + wght*urt(:,i).^2.*jac_v.*he_v(:,i);
            z_loc_norm = z_loc_norm + wght*zrt(:,i).^2.*jac_v.*he_v(:,i);
        end
    end
    % Divide by 2 due to quadrature rule.
     edgnorm_u(:,k) = u_loc_norm(:)/2;
     edgnorm_z(:,k) = z_loc_norm(:)/2;
end

end  % end function