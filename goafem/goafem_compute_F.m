function [Bform] = goafem_compute_F(KL_DATA, paras_fem, paras_sg, ...
    list, listy, listy2, sols_u, sols_z, rhs_fun, qoi_fun, rf_type, X, fun_p, rule_id, aa, polys)
% GOAFEM_DOBILINEARFORM computes the bilinear form associated with the
% correction term for the error in the goal functional.
%
% Bform = goafem_doBilinearForm(KL_DATA, paras_sg, paras_fem, ...
%    list, listy, sols_u, sols_z, rhs_fun, qoi_fun)
%
% input:
%          KL_DATA    function cell for the expanded diffusion coefficient
%        paras_fem    spatial mesh parameter cell array
%         paras_sg    sparse grid parameter cell array
%             list    integral cell array of products of Lagrange functions
%            listy    modified cell array associated with (linear) expansions of diffusion coeffients
%           listy2    modified cell array associated with (quadratic) expansions of diffusion coeffients
%           sols_u    per-collocation point primal solution vectors
%           sols_z    per-collocation point dual solution vectors
%           rhs_fun   function handle cell array for the right-hand-side
%           qoi_fun   function handle cell array for the quantity of interest
%           rf_type   expansion type
%
% output:
%             Bform   bilinear form B(sols_u, sols_z)
%
% Function(s) called: goafem_stochcol_Mmatrices
%                     goafem_stochcol_fem_setup
%
% SEE ALSO: goafem_singlelevelSC
%
% TR; 14 October 2023

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
        Bform = IntA'*IntY;
end