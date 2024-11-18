function [Fform] = goafem_compute_F(paras_fem, paras_sg, ...
    sols_z, rhs_fun, qoi_fun, X, fun_p, rule_id, aa, polys)
% GOAFEM_COMPUTE_F computes the linear functional associated with the
% correction term for the error in the goal functional.
%
% [Fform] = goafem_compute_F(paras_fem, paras_sg, ...
%    sols_z, rhs_fun, qoi_fun, X, fun_p, rule_id, aa, polys)
%
% input:
%        paras_fem    spatial mesh parameter cell array
%         paras_sg    sparse grid parameter cell array
%           sols_z    per-collocation point dual solution vectors
%           rhs_fun   function handle cell array for the right-hand-side
%           qoi_fun   function handle cell array for the quantity of interest
%           X         index set
%           rule_id   rule for generation of sampling points
%           aa        diffusion coefficient handle
%           polys     list of lagrange polynomials
%
% output:
%             Fform   linear functional F(sols_z)
%

%
% SEE ALSO: goafem_singlelevelSC
%
% AS; 28 June 2024

M = size(paras_sg{9},2);


        X_diff = stochcol_rmargin_set(X, stochcol_margin_set(X), 0);
        X_reft = stochcol_getgrid([X; X_diff]);
        paras_sg_reft = stochcol_sg(X_reft, rule_id);
        coords_reft = paras_sg_reft{9};     
        
        gridd = paras_sg{4};
        clincombiset = paras_sg{5};
        indlincombiset = paras_sg{6};       
        coords = paras_sg{9};
        
        K = size(coords,1);
        K2 = size(coords_reft,1);      
        IntA = zeros(K2,1);
        IntY = zeros(K2,1);
        [~, KA, ~] = intersect(coords_reft, coords, 'rows');
        [~,~,~,~,IntY(:)] = stochcol_multilag(paras_sg_reft{4}, paras_sg_reft{1}, paras_sg_reft{2}, paras_sg_reft{3}, rule_id, fun_p);
        for k = 1:K2
            rhs_fun_tmp = rhs_fun;
            rhs_fun_tmp{1} = @(x1, x2) rhs_fun{1}(x1, x2, coords_reft(k,:));
            [~, zgalreft, ~, ~, ~, fk] = goafem_stochcol_fem_solver(coords_reft(k,:), ...
                            paras_fem, aa, rhs_fun_tmp, qoi_fun);
            if sum(k == KA)
                IntA(k) = fk' * zgalreft;
            else
                Intk = zeros(K, 1);
                for j = 1:K
                    z_gal = sols_z(:,j);
                    Intk(j) = fk' * z_gal;
                end
                Lk = stochcol_getinterpolant_2(gridd, clincombiset, indlincombiset, coords_reft(k,:), polys);
                %Lk = Lk*Lk';
                IntA(k) = Lk'*Intk; %sum(Lk.*Intk,[1,2]);
            end
        end
        Fform = IntA'*IntY;
end