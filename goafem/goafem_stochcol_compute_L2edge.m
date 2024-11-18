function [edgnorm_u,edgnorm_z] = goafem_stochcol_compute_L2edge(sols_u_diff, sols_z_diff, evt)
% GOAFEM_STOCHCOL_COMPUTE_L2EDGE computes local norms for both primal and
% dual problems
%
% input:
%     sols_u_diff     pirmal solution at new colocation points
%     sols_z_diff     dual solution at new colocation points
%             evt     elements to vertex mapping
%
% output:
%       edgnorm_u     array of edge norm for primal solution 
%       edgnorm_z     array of edge norm for dual solution 
%
% TR; 28 September 2022

% Initialisation
ncoords = zeros((size(sols_u_diff,2)));
nel = size(evt,1);

edgnorm_u = zeros(nel,size(coords_diff,1));
edgnorm_z = zeros(nel,size(coords_diff,1));

% Construct the integration rule
nngpt = 7;
[oneg,onew] = gausspoints_oned(nngpt);

for k = 1:ncoords % loop over collocation points
    u_gal = sols_u_diff(:,k);
    z_gal = sols_z_diff(:,k);
    u_gal_loc = zeros(nel,3);
    z_gal_loc = zeros(nel,3);
    for ivtx = 1:3 % loop over local nodes
        u_gal_loc(:,ivtx) = u_gal(evt(:,ivtx));
        z_gal_loc(:,ivtx) = z_gal(evt(:,ivtx));
    end
    u_loc_norm = zeros(nel,3);
    z_loc_norm = zeros(nel,3);
    for igpt = 1:nngpt % loop over gauss points
        sigpt = oneg(igpt);
        wght = onew(igpt);
        s = (1.0 + sigpt)/2.0;
        smat = [s, 1-s; 0, s; s, 0];
        urt = 0;
        zrt = 0;
        for ivtx = 1:3
            [jac_v,~,phi_v,~,~] = tderiv(smat(ivtx,1), smat(ivtx,2), xl_v, yl_v);
             urt = urt + u_gal_loc(:,ivtx).*phi_v(:,ivtx);
             zrt = zrt + z_gal_loc(:,ivtx).*phi_v(:,ivtx);
        end
        u_loc_norm = u_loc_norm + wght*urt.^2.*jac_v;
        z_loc_norm = z_loc_norm + wght*zrt.^2.*jac_v;
    end
    edgnorm_u(:,k) = u_loc_norm(:);
    edgnorm_z(:,k) = z_loc_norm(:);
end

end  % end function