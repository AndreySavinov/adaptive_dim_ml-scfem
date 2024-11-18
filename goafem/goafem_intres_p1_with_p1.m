function [intresu,intresz] = goafem_intres_p1_with_p1(xy, xl_s, yl_s, ...
                             evt, u_gal, z_gal, subdivPar, gradcoeffx_fun, ...
                             gradcoeffy_fun, rhs_fun, qoi_fun)
%   [intresu,intresz] = goafem_intres_p1_with_p1(xy, xl_s, yl_s, ...
%                       evt, u_gal, z_gal, subdivPar, gradcoeffx_fun, ...
%                       gradcoeffy_fun, rhs_fun, qoi_fun)
%
%   input:
%               xy    vertex coordinate vector
%             xl_s    x-coordinates physical sub-elements
%             yl_s    y-coordinates physical sub-elements
%              evt    element mapping matrix
%            u_gal    vertex primal solution vector
%            z_gal    vertex dual solution vector
%        subdivPar    red/bisec3 uniform sub-division switch
%   gradcoeffx_fun    function handle for the first component of the gradient of the coefficient
%   gradcoeffy_fun    function handle for the second component of the gradient of the coefficient
%          rhs_fun    function handle for the right-hand side
%          qoi_fun    function handle for the quantity of interest
%   output:
%          intresu    interior residual (primal)
%          intresz    interior residual (dual)
%
% Function(s) called: triangular_gausspoints
%                     subelt_transf
%                     tderiv
%                     tgauss
%
% TR; 13 July 2022

x = xy(:,1);
y = xy(:,2);
nel = size(evt,1);

% Construct the integration rule
nngpt = 7;
[s,t,wt] = triangular_gausspoints(nngpt);

% Recover local coordinates and local solution
xl_v = zeros(nel,3);
yl_v = zeros(nel,3);
sl_u = zeros(nel,3);
sl_z = zeros(nel,3);
for ivtx = 1:3
    xl_v(:,ivtx) = x(evt(:,ivtx));
    yl_v(:,ivtx) = y(evt(:,ivtx));
    sl_u(:,ivtx) = u_gal(evt(:,ivtx));
    sl_z(:,ivtx) = z_gal(evt(:,ivtx));
end

% Preallocate matrices
intresu = zeros(nel,3);
intresz = zeros(nel,3);
%
bdem = zeros(nel,4,3,3);
bde = zeros(nel,3,3);
gdem = zeros(nel,4,3,3);
gde = zeros(nel,3,3);
%
fdem = zeros(nel,4,3);
fde = zeros(nel,3);
%
xl_m = zeros(nel,3);
yl_m = zeros(nel,3);

% Loop over sub-elements
for subelt = 1:4
    % Recover local coordinates over sub-element
    for ivtx = 1:3
        xl_m(:,ivtx) = xl_s(:,subelt,ivtx);
        yl_m(:,ivtx) = yl_s(:,subelt,ivtx);
    end
    % Loop over Gauss points
    for igpt = 1:nngpt
        sigpt = s(igpt);
        tigpt = t(igpt);
        wght = wt(igpt);
        [sigptloc,tigptloc] = subelt_transf(sigpt,tigpt,subelt,subdivPar);
        
        % Evaluate derivatives, gradient of the coefficients, and rhs
        [~,invjac_v,~,dphidx_v,dphidy_v] = tderiv(sigptloc,tigptloc,xl_v,yl_v);
        [jac_m, ~, phi_m, ~, ~] = tderiv(sigpt, tigpt, xl_m, yl_m);
        [diffx] = tgauss(sigpt, tigpt, xl_m, yl_m, gradcoeffx_fun);
        [diffy] = tgauss(sigpt, tigpt, xl_m, yl_m, gradcoeffy_fun);
        [rhs_m, qoi_m] = goafem_stochcol_tgauss(sigpt, tigpt, xl_m, yl_m, rhs_fun, qoi_fun);
        
        % Loop over mid-edge hat functions
        for j = 1:3
            
            % Compute rhs-contribution from the source ( L2 <=(:,1) and DIV-H1 <=(:,4) )
            fdem(:,subelt,j) = fdem(:,subelt,j) + wght * (rhs_m(:,1) + rhs_m(:,4)) .* phi_m(:,j) .* jac_m(:);
            gdem(:,subelt,j) = gdem(:,subelt,j) + wght * (qoi_m(:,1) + qoi_m(:,4)) .* phi_m(:,j) .* jac_m(:);
            
            % Compute div(a*grad)-contribution = grad(a)*grad(u_h): loop
            % over vertices old hat functions
            for i = 1:3
                bdem(:,subelt,j,i) = bdem(:,subelt,j,i) + wght * diffx(:) .* dphidx_v(:,i) .* phi_m(:,j) .* invjac_v(:) .* jac_m(:);
                bdem(:,subelt,j,i) = bdem(:,subelt,j,i) + wght * diffy(:) .* dphidy_v(:,i) .* phi_m(:,j) .* invjac_v(:) .* jac_m(:);
            end
            % end vertices old hat functions loop
        end
        % end mid-edge hat functions loop
    end
    % end Gauss points loop
end
% end sub-elements loop

% -----------------------------------------------------------------------------
% Manual assembly of subelement contributions
% -----------------------------------------------------------------------------
if subdivPar == 1
    %
    % Red sub-division: assembling
    %
    % First edge
    bde(:,1,1) = bdem(:,2,3,1) + bdem(:,3,2,1) + bdem(:,4,1,1);
    bde(:,1,2) = bdem(:,2,3,2) + bdem(:,3,2,2) + bdem(:,4,1,2);
    bde(:,1,3) = bdem(:,2,3,3) + bdem(:,3,2,3) + bdem(:,4,1,3);
    fde(:,1)   = fdem(:,2,3)   + fdem(:,3,2)   + fdem(:,4,1);
    gde(:,1)   = gdem(:,2,3)   + gdem(:,3,2)   + gdem(:,4,1);
    % Second edge
    bde(:,2,1) = bdem(:,1,3,1) + bdem(:,3,1,1) + bdem(:,4,2,1);
    bde(:,2,2) = bdem(:,1,3,2) + bdem(:,3,1,2) + bdem(:,4,2,2);
    bde(:,2,3) = bdem(:,1,3,3) + bdem(:,3,1,3) + bdem(:,4,2,3);
    fde(:,2)   = fdem(:,1,3)   + fdem(:,3,1)   + fdem(:,4,2);
    gde(:,2)   = gdem(:,1,3)   + gdem(:,3,1)   + gdem(:,4,2);
    % Third edge
    bde(:,3,1) = bdem(:,1,2,1) + bdem(:,2,1,1) + bdem(:,4,3,1);
    bde(:,3,2) = bdem(:,1,2,2) + bdem(:,2,1,2) + bdem(:,4,3,2);
    bde(:,3,3) = bdem(:,1,2,3) + bdem(:,2,1,3) + bdem(:,4,3,3);
    fde(:,3)   = fdem(:,1,2)   + fdem(:,2,1)   + fdem(:,4,3);
    gde(:,3)   = gdem(:,1,2)   + gdem(:,2,1)   + gdem(:,4,3);
    
else
    %
    % Bisec3 sub-division: assembling
    %
    % First edge
    bde(:,1,1) = bdem(:,3,2,1) + bdem(:,4,2,1);
    bde(:,1,2) = bdem(:,3,2,2) + bdem(:,4,2,2);
    bde(:,1,3) = bdem(:,3,2,3) + bdem(:,4,2,3);
    fde(:,1)   = fdem(:,3,2)   + fdem(:,4,2);
    gde(:,1)   = gdem(:,3,2)   + gdem(:,4,2);
    % Second edge
    bde(:,2,1) = bdem(:,1,3,1) + bdem(:,2,1,1) + bdem(:,3,3,1) + bdem(:,4,1,1);
    bde(:,2,2) = bdem(:,1,3,2) + bdem(:,2,1,2) + bdem(:,3,3,2) + bdem(:,4,1,2);
    bde(:,2,3) = bdem(:,1,3,3) + bdem(:,2,1,3) + bdem(:,3,3,3) + bdem(:,4,1,3);
    fde(:,2)   = fdem(:,1,3)   + fdem(:,2,1)   + fdem(:,3,3)   + fdem(:,4,1);
    gde(:,2)   = gdem(:,1,3)   + gdem(:,2,1)   + gdem(:,3,3)   + gdem(:,4,1);
    % Third edge
    bde(:,3,1) = bdem(:,1,2,1) + bdem(:,2,2,1);
    bde(:,3,2) = bdem(:,1,2,2) + bdem(:,2,2,2);
    bde(:,3,3) = bdem(:,1,2,3) + bdem(:,2,2,3);
    fde(:,3)   = fdem(:,1,2)   + fdem(:,2,2);
    gde(:,3)   = gdem(:,1,2)   + gdem(:,2,2);
end

% Assemble interior residuals from rhs source contribution
for i = 1:3
    intresu(:,i) = intresu(:,i) + fde(:,i);
    intresz(:,i) = intresz(:,i) + gde(:,i);
end

% Assemble by adding the div(a*grad)-contribution times Galerkin solution
for j = 1:3
    for k = 1:3
        intresu(:,j) = intresu(:,j) + bde(:,j,k).*sl_u(:,k);
        intresz(:,j) = intresz(:,j) + bde(:,j,k).*sl_z(:,k);
    end
end

end  % end function