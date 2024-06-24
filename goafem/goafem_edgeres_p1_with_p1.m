function [edgeres_u,edgeres_z] = goafem_edgeres_p1_with_p1(xy, evt, eboundt, u_gal, z_gal, eex, ...
    tve, els, coeff_fun)
% COMMENTS NEED FINISHING
%
%   input:
%               xy    vertex coordinate vector
%              evt    element mapping matrix
%          eboundt    element boundary mapping matrix
%            u_gal    vertex primal solution vector
%            z_gal    vertex dual solution vector
%              eex    element connectivity array
%              tve    edge location array
%              els    elementwise edge lengths
%        coeff_fun    function handle for the diffusion coefficient
%
%   output:
%          edgeres_u   edge residuals
%          edgeres_z   edge residuals
%
% Function(s) called: gausspoints_oned
%                     tderiv
%                     p1fluxjmps

x = xy(:,1);
y = xy(:,2);
nel = size(evt,1);

% Construct the 1D integration rule
ngpt = 7;
[oneg,onew] = gausspoints_oned(ngpt);

% Initialisation
edgeres_u = zeros(nel,3);
edgeres_z = zeros(nel,3);

xl_v = zeros(nel,3);
yl_v = zeros(nel,3);
% Recover local coordinates
for ivtx = 1:3
    xl_v(:,ivtx) = x(evt(:,ivtx));
    yl_v(:,ivtx) = y(evt(:,ivtx));
end

%fprintf('computing P1 flux jumps... ')

% Loop over Gaussian points
for igpt = 1:ngpt
    
    sigpt = oneg(igpt);
    sigpt_ref_edge = (1.0 + sigpt)/2.0; % [-1,1] -> [0,1]    reference-edge map
    sigpt_l = (1.0 + sigpt)/4.0;        % [-1,1] -> [0,1/2]   LEFT sub-edge map
    sigpt_r = (3.0 + sigpt)/4.0;        % [-1,1] -> [1/2,1]  RIGHT sub-edge map
    wigpt = onew(igpt);
    
    % First edge
    [~,~,phi_v_1,~,~] = tderiv(sigpt_ref_edge,1-sigpt_ref_edge,xl_v,yl_v);
    
    % Second edge
    [~,~,phi_v_2,~,~] = tderiv(0,sigpt_ref_edge,xl_v,yl_v);
    
    % Third edge
    [~,~,phi_v_3,~,~] = tderiv(sigpt_ref_edge,0,xl_v,yl_v);
    
    % Jump of the finite element solution over the left and right sub-
    % element edges for the primal and dual problems
    [njmp_l_u] = p1fluxjmps(u_gal, eex, xy, evt, eboundt, tve, sigpt_l, ...
        coeff_fun);
    [njmp_r_u] = p1fluxjmps(u_gal, eex, xy, evt, eboundt, tve, sigpt_r, ...
        coeff_fun);
    [njmp_l_z] = p1fluxjmps(z_gal, eex, xy, evt, eboundt, tve, sigpt_l, ...
        coeff_fun);
    [njmp_r_z] = p1fluxjmps(z_gal, eex, xy, evt, eboundt, tve, sigpt_r, ...
        coeff_fun);
    
    % Contribution of the first edge
    edgeres_u(:,1) = edgeres_u(:,1) + (1/2) * wigpt * njmp_l_u(:,1) .* phi_v_1(:,2) .* (els(:,1)/4);
    edgeres_u(:,1) = edgeres_u(:,1) + (1/2) * wigpt * njmp_r_u(:,1) .* phi_v_1(:,3) .* (els(:,1)/4);
    edgeres_z(:,1) = edgeres_z(:,1) + (1/2) * wigpt * njmp_l_z(:,1) .* phi_v_1(:,2) .* (els(:,1)/4);
    edgeres_z(:,1) = edgeres_z(:,1) + (1/2) * wigpt * njmp_r_z(:,1) .* phi_v_1(:,3) .* (els(:,1)/4);
    
    % Contribution of the second edge
    edgeres_u(:,2) = edgeres_u(:,2) + (1/2) * wigpt * njmp_l_u(:,2) .* phi_v_2(:,3) .* (els(:,2)/4);
    edgeres_u(:,2) = edgeres_u(:,2) + (1/2) * wigpt * njmp_r_u(:,2) .* phi_v_2(:,1) .* (els(:,2)/4);
    edgeres_z(:,2) = edgeres_z(:,2) + (1/2) * wigpt * njmp_l_z(:,2) .* phi_v_2(:,3) .* (els(:,2)/4);
    edgeres_z(:,2) = edgeres_z(:,2) + (1/2) * wigpt * njmp_r_z(:,2) .* phi_v_2(:,1) .* (els(:,2)/4);
    
    % Contribution of the third edge
    edgeres_u(:,3) = edgeres_u(:,3) + (1/2) * wigpt * njmp_l_u(:,3) .* phi_v_3(:,2) .* (els(:,3)/4);
    edgeres_u(:,3) = edgeres_u(:,3) + (1/2) * wigpt * njmp_r_u(:,3) .* phi_v_3(:,1) .* (els(:,3)/4);
    edgeres_z(:,3) = edgeres_z(:,3) + (1/2) * wigpt * njmp_l_z(:,3) .* phi_v_3(:,2) .* (els(:,3)/4);
    edgeres_z(:,3) = edgeres_z(:,3) + (1/2) * wigpt * njmp_r_z(:,3) .* phi_v_3(:,1) .* (els(:,3)/4);
    
end
% end loop over Gaussian points

%fprintf('done\n')

end  % end function
