function [errest_u, errest_z, elerr_u, elerr_z, ederr_u, ederr_z, ...
            Msetu, Msetz] = ...
          goafem_stochcol_fem_estimator(pmethod, ...
            paras_fem_errest, paras_fem, paras_detail, u_gal, z_gal, ...
            coord, aa, aax1, aax2, rhs_fun, qoi_fun, tout)
%GOAFEM_STOCHCOL_FEM_ESTIMATOR computes spatial errors for given
%collocation point coord
%       
% input:
%            pmethod    defines P1/P2 approximation
%   paras_fem_errest, 
%          paras_fem,
%       paras_detail    parameters of FEM and cooresponding mesh
%              u_gal    primal solution at collocation point coord
%              z_gal    dual solution at collocation point coord
%              coord    collocation point
%     aa, aax1, aax2    handles for dissusion coefficient and componets of its gradient
%            rhs_fun    function handle for the right-hand side
%            qoi_fun    function handle for the quantity of interest
%               tout    switcher to posteriori error estimation
% output:
%            elerr_u    vector of 2-level element primal indicators
%            elerr_z    vector of 2-level element dual indicators
%            ederr_u    vector of 2-level edge primal indicators
%            ederr_z    vector of 2-level edge dual indicators
%           errest_u    global 2-level primal error estimate
%           errest_z    global 2-level dual error estimate
%              Msetu    set of marked elements/edges for primal problem
%              Msetz    set of marked elements/edges for dual problem
%
% SEE ALSO: diffpost_p1_with_p1_2level
%
% TR; 13 July 2022

if nargin < 13, tout = 0; end
xy      = paras_fem{1};
evt     = paras_fem{2};
eboundt = paras_fem{4};
eex     = paras_detail{1};
tve     = paras_detail{2};
els     = paras_detail{3};
evtY    = paras_detail{4};
xyY     = paras_detail{5};
boundY  = paras_detail{6};
Ybasis  = paras_detail{7};
coeff_fun = @(x1, x2) aa(x1, x2, coord);
gradcoeffx1_fun = @(x1, x2) aax1(x1, x2, coord);
gradcoeffx2_fun = @(x1, x2) aax2(x1, x2, coord);
subdivPar = paras_fem_errest(1);

% -------------------------------------------------------------------
% ESTIMATE
% -------------------------------------------------------------------
if tout
fprintf('\nA posteriori error estimation\n'); end
%
% Compute a posteriori error estimation
errorTime = tic;
%
if pmethod == 1
    % P1-error estimation
    pestim = paras_fem_errest(2);
    if pestim == 1
        estimtype = paras_fem_errest(3);
        markedgelem = paras_fem_errest(4);
        markstrat = paras_fem_errest(5);
        threshold = paras_fem_errest(6);
        % Using linear midpoint Hat functions
        if tout
        fprintf('Error estimation using 3 edge midpoint linear functions\n'); end
        if estimtype == 1
            % Hierarchical eY estimator: elementwise residual problem
            if tout
            fprintf('Hierarchical eY estimator: solving elementwise residual problems\n'); end
            [elerr_u, elerr_z, fe, ge, ae] = ...
                goafem_diffpost_p1_with_p1(xy, evt, eex, tve, ...
                els, eboundt, u_gal, z_gal, subdivPar, coeff_fun, ...
                gradcoeffx1_fun, gradcoeffx2_fun, rhs_fun, qoi_fun);
            [~, elerr_u] = diffpost_p1_bc(ae, fe, elerr_u, xy, evt, eboundt);
            [~, elerr_z] = diffpost_p1_bc(ae, ge, elerr_z, xy, evt, eboundt);
            %
            % Global error estimate
            errest_u = norm(elerr_u,2);
            errest_z = norm(elerr_z,2);
        elseif estimtype == 2
            % Hierarchical eY estimator: solving the assembled linear system
            if tout
            fprintf('Hierarchical eY estimator: solving assembled linear system\n'); end
            [elerr_u, ederr_u, elerr_z, ederr_z, errest_u, errest_z] = ...
                goafem_diffpost_p1_with_p1_linsys(xy, evt, eboundt, u_gal, z_gal, ...
                evtY, xyY, boundY, Ybasis, subdivPar, coeff_fun, rhs_fun, qoi_fun);
        else% estimtype == 3
            % 2-level error estimators
            if tout
            fprintf('Two-level estimator\n'); end
            [elerr_u, ederr_u, elerr_z, ederr_z, errest_u, errest_z] = ...
                goafem_diffpost_p1_with_p1_2level(xy, evt, eboundt, ...
                u_gal, z_gal, evtY, xyY, boundY, Ybasis, subdivPar, ...
                coeff_fun, rhs_fun, qoi_fun);
        end
    elseif pestim == 2
        markedgelem = paras_fem_errest(3);
        markstrat = paras_fem_errest(4);
        threshold = paras_fem_errest(5);
        % Using quadratic Midpoint Bubble functions
        if tout
        fprintf('Error estimation using 4 quadratic bubble functions\n'); end
        [elerr_u, elerr_z, fe, ge, ae] = ...
            goafem_diffpost_p1_with_p2(xy, evt, eex, tve, ...
            els, eboundt, u_gal, z_gal, coeff_fun, gradcoeffx1_fun, ...
            gradcoeffx2_fun, rhs_fun, qoi_fun);
        [~, elerr_u] = diffpost_p1_bc(ae, fe, elerr_u, xy, evt, eboundt);
        [~, elerr_z] = diffpost_p1_bc(ae, ge, elerr_z, xy, evt, eboundt);
        %
        % Global error estimate
        errest_u = norm(elerr_u,2);
        errest_z = norm(elerr_z,2);
    end
else% pmethod == 2
    markedgelem = paras_fem_errest(2);
    markstrat = paras_fem_errest(3);
    threshold = paras_fem_errest(4);
    % P2-error estimation
    if tout
    fprintf('Error estimation using quartic bubble functions\n'); end
    %
    % check performance with that of legacy code
    [elerr_u, elerr_z, ~, ~, ~] = ...
    goafem_diffpost_p2_with_p4(xy, evt, eex, tve, els, eboundt, ...
    u_gal, z_gal, coeff_fun, rhs_fun, qoi_fun);
    %
    % Global error estimate
    errest_u = norm(elerr_u,2);
    errest_z = norm(elerr_z,2);
end % end if
%
if tout
fprintf('Estimated energy error (primal): %10.4e\n',errest_u);
fprintf('Estimated energy error (dual): %10.4e\n',errest_z);
fprintf('Estimation took %.5f sec\n',toc(errorTime)); 
end

% -------------------------------------------------------------------
% MARK (need to incorporate dual properly)
% -------------------------------------------------------------------
if markedgelem == 1
    % Marking elements
    [Msetu] = marking_strategy_fa(elerr_u, markstrat, threshold);
    [Msetz] = marking_strategy_fa(elerr_z, markstrat, threshold);
    
else%markedgelem == 2
    % Marking edges
    [Msetu] = marking_strategy_fa(ederr_u, markstrat, threshold);
    [Msetz] = marking_strategy_fa(ederr_z, markstrat, threshold);
end

end
