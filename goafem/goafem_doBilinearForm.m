 function [Bform] = goafem_doBilinearForm(KL_DATA, paras_fem, paras_sg, ...
    list, listy, listy2, sols_u, sols_z, rhs_fun, qoi_fun, rf_type, X, fun_p, rule_id, aa, polys, varargin)
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

if nargin(rhs_fun{1}) == 3
    rhs_fun{1} = @(x1, x2) rhs_fun{1}(x1, x2, [0, 0]);
end
switch rf_type
    
    case 3 % linear coefficient
        M = size(paras_sg{9},2);
        Am = cell(M+1,1);
        IntY = cell(M+1,1);
        IntA = cell(M+1,1);
        Bform = 0;

        for m = 1:M+1
            IntY{m} = goafem_stochcol_Mmatrices(paras_sg{4}, paras_sg{5}, ...
                            paras_sg{6}, list, m-1, listy);     
            [Am{m}, ~, ~] = goafem_stochcol_fem_setup(paras_fem{1}, paras_fem{2}, ...
                            KL_DATA.coeff{m}, rhs_fun);
            IntA{m} = sols_u' * Am{m} * sols_z;
            Bform = Bform + sum(dot(IntY{m}, sols_u' * Am{m} * sols_z));
        end
        
    case 2 % quadratic coefficient
        Am = cell(M+1,1);
        IntY = cell(M+1,1);
        am = cell(M+1,1);
        IntA = cell(M+1,1);
        Bform1 = 0;
        Bform2 = 0;

        for m = 1:M+1 % a0^2 and a0*am contributions
            IntY{m} = goafem_stochcol_Mmatrices(paras_sg{4}, paras_sg{5}, ...
                            paras_sg{6}, list, m-1, listy);
            am{m} = @(x,y) KL_DATA.coeff{m}(x,y).*KL_DATA.coeff{1}(x,y);
            [Am{m}, ~, ~] = goafem_stochcol_fem_setup(paras_fem{1}, paras_fem{2}, ...
                            am{m}, rhs_fun);
             if m > 1
                 Am{m} = 2*Am{m};
             end
            IntA{m} = sols_u' * Am{m} * sols_z;
%             IntA{m} = IntA{m}.*(abs(IntA{m}) > 1e-16);
            Bform1 = Bform1 + sum(dot(IntY{m}, sols_u' * Am{m} * sols_z));
        end
        
        am2 = cell(M,M);
        Am2 = cell(M,M);
        IntY2 = cell(M,M);
        IntA2 = cell(M,M);
        for m = 1:M % am*am' and am^2 contribution
            for m2 = m:M
                IntY2{m,m2} = goafem_stochcol_M2matrices(paras_sg{4}, paras_sg{5}, ...
                                paras_sg{6}, list, m, m2, listy, listy2);
                am2{m,m2} = @(x,y) KL_DATA.coeff{m+1}(x,y).*KL_DATA.coeff{m2+1}(x,y); 
                [Am2{m,m2}, ~, ~] = goafem_stochcol_fem_setup(paras_fem{1}, paras_fem{2}, ...
                            am2{m,m2}, rhs_fun);
                if m2 < m
                    Am2{m,m2} = 2*Am2{m,m2};
                end
                IntA2{m,m2} = sols_u' * Am2{m,m2} * sols_z;
%                 IntA2{m,m2} = IntA2{m,m2}.*(abs(IntA2{m,m2}) > 1e-16);
                Bform2 = Bform2 + sum(dot(IntY2{m,m2}, sols_u' * Am2{m,m2} * sols_z),[1,2]);
            end
        end
        
        % [Bform1, Bform2]
        Bform = Bform1 + Bform2;
        
    otherwise % incl. exponential coefficient (rf_type = 1)
        if numel(varargin) > 0
            Q = varargin{1};
        end
        X_diff = stochcol_rmargin_set(X, stochcol_margin_set(X), Q);
        X_reft = stochcol_getgrid([X, ones(size(X, 1), Q); X_diff]);
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
        parfor k = 1:K2
            [ugalreft, zgalreft, Ak] = goafem_stochcol_fem_solver(coords_reft(k,:), ...
                            paras_fem, aa, rhs_fun, qoi_fun);
            if sum(k == KA)
                IntA(k) = ugalreft' * Ak * zgalreft;
            else
                Intk = zeros(K,K);
                for i = 1:K
                    for j = 1:K
                        u_gal = sols_u(:,i);
                        z_gal = sols_z(:,j);
                        Intk(i,j) = u_gal' * Ak * z_gal;
                    end
                end
                Lk = stochcol_getinterpolant_2(gridd, clincombiset, indlincombiset, coords_reft(k,:), polys);
                Lk = Lk*Lk';
                IntA(k) = sum(Lk.*Intk,[1,2]);
            end
        end
        Bform = IntA'*IntY;
end

return
