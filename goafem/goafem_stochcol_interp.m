function [sols_u_interp, sols_z_interp] = goafem_stochcol_interp(solsu, solsz, ...
    paras_sg, paras_fem, sols_ref_u, sols_ref_z, paras_sg_ref, ...
    paras_fem_ref, polys)
% GOAFEM_STOCHCOL_INTERP interpolates a solution to a reference solution
%
%   input:  
%             solsu   primal solution
%             solsz   dual solution
%          paras_sg   sparse grid parameters
%         paras_fem   spatial mesh parameter
%        sols_ref_u   reference primal solution
%        sols_ref_z   reference dual solution
%      paras_sg_ref   reference sparse grid parameters
%     paras_fem_ref   reference spatial mesh parameter
%             polys   list of Lagrange polynomials
%            listy2   modified cell array associated with (quadratic) expansions of diffusion coeffients
%
%   output:
%     sols_u_interp   interpolated primal solution
%     sols_u_interp   interpolated dual solution
%
%
% TR; 28 September 2022

gridd          = paras_sg{4};
clincombiset   = paras_sg{5};
indlincombiset = paras_sg{6};
coords_ref = paras_sg_ref{end};
[N_ref, ~] = size(coords_ref);
xy        = paras_fem{1};
xy_ref    = paras_fem_ref{1};
% interpolate the approximation in reference collocation points and
% reference finite element grid
sols_u_interp = zeros(size(sols_ref_u));
sols_z_interp = zeros(size(sols_ref_z));
for k = 1:N_ref
    LL = stochcol_getinterpolant_2(gridd, ...
        clincombiset, indlincombiset, coords_ref(k,:), polys);
    uyy = solsu*LL;
    zyy = solsz*LL;
    uyy_interp = griddata(xy(:,1), xy(:,2), uyy, xy_ref(:,1), xy_ref(:,2));
    zyy_interp = griddata(xy(:,1), xy(:,2), zyy, xy_ref(:,1), xy_ref(:,2));
    sols_u_interp(:,k) = uyy_interp;
    sols_z_interp(:,k) = zyy_interp;
end

end
