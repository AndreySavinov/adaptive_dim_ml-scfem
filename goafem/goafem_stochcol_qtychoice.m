% GOAFEM_STOCHCOL_QTYCHOICE
% Function handles w.r.t. spatial variables are defined according to an
% L2/H1 splitting, as in the examples given by Mommer, Stevenson.
%
% TR; 13 September 2023

% Function handles for the RHS {L2, H1(1), H1(2), DIV-H1}
fprintf('\n Choose right-hand side function:\n')
fprintf('1. Unit RHS function \n')
fprintf('2. Mommer-Stevenson based Examplse \n')
% 2.1. Modified Mommer-Stevenson type problem with added unit L2 part
% 2.2. Variable Mommer-Stevenson (define cut-off variable)
fprintf('3. Variable RHS \n')
fprintf('4. One-Peak problem \n')
fprintf('Or have a look into goafem_stochcol_qtychoice.m for other choices \n')
% 3.1. Quadratic variable RHS
% 3.2. Exponential variable RHS
% 9.9. Unsmiley (easteregg)
RHS_type = default('(default is 1)',1);


switch RHS_type
    case 1 % Unit RHS function
        F_rhs = {{@(x1,x2) ones(size(x1)), @(x1,x2) zeros(size(x1)), @(x1,x2) zeros(size(x1)), @(x1,x2) zeros(size(x1))}, ...
                 {@(x1,x2) ones(size(x1)), @(x1,x2) zeros(size(x1)), @(x1,x2) zeros(size(x1)), @(x1,x2) zeros(size(x1))}, ...
                 {@(x1,x2) ones(size(x1)), @(x1,x2) zeros(size(x1)), @(x1,x2) zeros(size(x1)), @(x1,x2) zeros(size(x1))}, ...
                 {@(x1,x2) ones(size(x1)), @(x1,x2) zeros(size(x1)), @(x1,x2) zeros(size(x1)), @(x1,x2) zeros(size(x1))}};
    case 2 % Mommer-Stevenson based example
        F_rhs = {{@(x1,x2) zeros(size(x1)), @(x1,x2) 1.*(x2 <= 0.5 - x1),  @(x1,x2) zeros(size(x1)), @(x1,x2) zeros(size(x1))}, ...
                 {@(x1,x2) zeros(size(x1)), @(x1,x2) 1.*(x2 <= -1.0 - x1), @(x1,x2) zeros(size(x1)), @(x1,x2) zeros(size(x1))}, ...
                 {@(x1,x2) zeros(size(x1)), @(x1,x2) 1.*(x2 >= 1.5 + x1),  @(x1,x2) zeros(size(x1)), @(x1,x2) zeros(size(x1))}, ...
                 {@(x1,x2) zeros(size(x1)), @(x1,x2) 1.*(x2 >= 1.5 + x1),  @(x1,x2) zeros(size(x1)), @(x1,x2) zeros(size(x1))}};
    case 2.1 % Mommer-Stevenson with unit L2
        F_rhs = {{@(x1,x2) ones(size(x1)), @(x1,x2) 1.*(x2 <= 0.5 - x1),  @(x1,x2) zeros(size(x1)), @(x1,x2) zeros(size(x1))}, ...
                 {@(x1,x2) ones(size(x1)), @(x1,x2) 1.*(x2 <= -1.0 - x1), @(x1,x2) zeros(size(x1)), @(x1,x2) zeros(size(x1))}, ...
                 {@(x1,x2) ones(size(x1)), @(x1,x2) 1.*(x2 >= 1.5 + x1),  @(x1,x2) zeros(size(x1)), @(x1,x2) zeros(size(x1))}, ...
                 {@(x1,x2) ones(size(x1)), @(x1,x2) 1.*(x2 >= 1.5 + x1),  @(x1,x2) zeros(size(x1)), @(x1,x2) zeros(size(x1))}};
    case 2.2 % Variable Mommer-Stevenson
        vms = default('Variable indicator sub-domain size (default is 0.5)', 0.5);
        mag = default('Magnitude of indicator (default is 1)', 1.0);
        F_rhs = {{@(x1,x2) zeros(size(x1)), @(x1,x2) mag.*(x2 <= vms - x1),      @(x1,x2) zeros(size(x1)), @(x1,x2) zeros(size(x1))}, ...
                 {@(x1,x2) zeros(size(x1)), @(x1,x2) mag.*(x2 <= -2 + vms - x1), @(x1,x2) zeros(size(x1)), @(x1,x2) zeros(size(x1))}, ...
                 {@(x1,x2) zeros(size(x1)), @(x1,x2) mag.*(x2 >= 2 - vms + x1),  @(x1,x2) zeros(size(x1)), @(x1,x2) zeros(size(x1))}, ...
                 {@(x1,x2) zeros(size(x1)), @(x1,x2) mag.*(x2 >= 2 - vms + x1),  @(x1,x2) zeros(size(x1)), @(x1,x2) zeros(size(x1))}};
    case 3 % Linear variable RHS
        F_rhs = {{@(x1,x2) x1 - x2, @(x1,x2) zeros(size(x1)), @(x1,x2) zeros(size(x1)), @(x1,x2) zeros(size(x1))}, ...
                 {@(x1,x2) x1 - x2, @(x1,x2) zeros(size(x1)), @(x1,x2) zeros(size(x1)), @(x1,x2) zeros(size(x1))}, ...
                 {@(x1,x2) x1 - x2, @(x1,x2) zeros(size(x1)), @(x1,x2) zeros(size(x1)), @(x1,x2) zeros(size(x1))}, ...
                 {@(x1,x2) x1 - x2, @(x1,x2) zeros(size(x1)), @(x1,x2) zeros(size(x1)), @(x1,x2) zeros(size(x1))}};
    case 3.1 % Quadratic variable RHS
        F_rhs = {{@(x1,x2) 2.*ones(size(x1)) - x1.^2 - x2.^2, @(x1,x2) zeros(size(x1)), @(x1,x2) zeros(size(x1)), @(x1,x2) zeros(size(x1))}, ...
                 {@(x1,x2) 2.*ones(size(x1)) - x1.^2 - x2.^2, @(x1,x2) zeros(size(x1)), @(x1,x2) zeros(size(x1)), @(x1,x2) zeros(size(x1))}, ...
                 {@(x1,x2) 2.*ones(size(x1)) - x1.^2 - x2.^2, @(x1,x2) zeros(size(x1)), @(x1,x2) zeros(size(x1)), @(x1,x2) zeros(size(x1))}, ...
                 {@(x1,x2) 2.*ones(size(x1)) - x1.^2 - x2.^2, @(x1,x2) zeros(size(x1)), @(x1,x2) zeros(size(x1)), @(x1,x2) zeros(size(x1))}};
    case 3.2 % Exponential variable RHS
        F_rhs = {{@(x1,x2) exp(-((x1 + 0.5).^2 + (x2 - 0.5).^2)), @(x1,x2) zeros(size(x1)), @(x1,x2) zeros(size(x1)), @(x1,x2) zeros(size(x1))}, ...
                 {@(x1,x2) exp(-((x1 + 0.5).^2 + (x2 - 0.5).^2)), @(x1,x2) zeros(size(x1)), @(x1,x2) zeros(size(x1)), @(x1,x2) zeros(size(x1))}, ...
                 {@(x1,x2) exp(-((x1 + 0.5).^2 + (x2 - 0.5).^2)), @(x1,x2) zeros(size(x1)), @(x1,x2) zeros(size(x1)), @(x1,x2) zeros(size(x1))}, ...
                 {@(x1,x2) exp(-((x1 + 0.5).^2 + (x2 - 0.5).^2)), @(x1,x2) zeros(size(x1)), @(x1,x2) zeros(size(x1)), @(x1,x2) zeros(size(x1))}};
    case 4 % One peak RHS
        F_rhs = {{@(x1,x2,y) one_peak_rhs(x1, x2, y), @(x1,x2) zeros(size(x1)), @(x1,x2) zeros(size(x1)), @(x1,x2) zeros(size(x1))}, ...
                 {@(x1,x2,y) one_peak_rhs(x1, x2, y), @(x1,x2) zeros(size(x1)), @(x1,x2) zeros(size(x1)), @(x1,x2) zeros(size(x1))}, ...
                 {@(x1,x2,y) one_peak_rhs(x1, x2, y), @(x1,x2) zeros(size(x1)), @(x1,x2) zeros(size(x1)), @(x1,x2) zeros(size(x1))}, ...
                 {@(x1,x2,y) one_peak_rhs(x1, x2, y), @(x1,x2) zeros(size(x1)), @(x1,x2) zeros(size(x1)), @(x1,x2) zeros(size(x1))}};
    case 9.9 % Unsmiley (easteregg)
        FACE_f = {@(x1,x2) 1.*(0.008 < (x1-0.7).^2+(x2-0.7).^2 & (x1-0.7).^2+(x2-0.7).^2 < 0.012) + ...
                           1.*(0.008 < (x1-0.3).^2+(x2-0.7).^2 & (x1-0.3).^2+(x2-0.7).^2 < 0.012), ...
                  @(x1,x2) 1.*(0.02 < (x1-0.4).^2+(x2-0.4).^2 & (x1-0.4).^2+(x2-0.4).^2 < 0.04) + ...
                           1.*(0.02 < (x1+0.4).^2+(x2-0.4).^2 & (x1+0.4).^2+(x2-0.4).^2 < 0.04), ...
                  @(x1,x2) 1.*((x1-0.1).^2+(x2-0.6).^2 < 0.03), ...
                  @(x1,x2) 1.*(0.02 < (x1-0.4).^2+(x2-0.4).^2 & (x1-0.4).^2+(x2-0.4).^2 < 0.04) + ...
                           1.*(0.02 < (x1+0.4).^2+(x2-0.4).^2 & (x1+0.4).^2+(x2-0.4).^2 < 0.04)};
        if dom_type == 1
            FACE_f = FACE_f{1};
        elseif dom_type == 5
            FACE_f = FACE_f{2};
        elseif dom_type == 2
            FACE_f = FACE_f{3}; 
        elseif dom_type == 4
            FACE_f = FACE_f{4};
        else
            error('For this domain the procedure is unavailable');
        end
        F_rhs = {{@(x1,x2) FACE_f(x1,x2),  @(x1,x2) zeros(size(x1)), @(x1,x2) zeros(size(x1)), @(x1,x2) zeros(size(x1))}, ...
                 {@(x1,x2) FACE_f(x1,x2), @(x1,x2) zeros(size(x1)), @(x1,x2) zeros(size(x1)), @(x1,x2) zeros(size(x1))}, ...
                 {@(x1,x2) FACE_f(x1,x2), @(x1,x2) zeros(size(x1)), @(x1,x2) zeros(size(x1)), @(x1,x2) zeros(size(x1))}, ...
                 {@(x1,x2) FACE_f(x1,x2), @(x1,x2) zeros(size(x1)), @(x1,x2) zeros(size(x1)), @(x1,x2) zeros(size(x1))}};
end

if dom_type == 1
    F_rhs = F_rhs{1};
elseif dom_type == 5
    F_rhs = F_rhs{2};
elseif dom_type == 2
    F_rhs = F_rhs{3}; 
elseif dom_type == 4
    F_rhs = F_rhs{4};
elseif dom_type == 3
    F_rhs = F_rhs{1};
else
    error('For this domain the procedure is unavailable');
end
    

% Function handles for the QOI {L2, H1(1), H1(2), DIV-H1}
fprintf('\n Choose quantity of interest:\n')
fprintf('1. Square Sub-Domain Integrals \n')
% 1. Rectangular Sub-Domain Integral
% 1.1. Unit QOI for Full Domain Integral
% 1.2. Circular Sub-Domain Integral
% 1.3. Elliptic Sub-Domain Integral
fprintf('2. Mollifiers for point-wise estimation (also redefines RHS fn.)\n')
fprintf('3. ( int w(x)u )^2 for Subdomain  \n')
fprintf('4. Convection term for Subdomain  \n')
fprintf('5. u^2 for Sub Domain Integral \n')
fprintf('6. Variance with Mollifier weight  \n')
fprintf('7. int u^2 for the enormous domain  \n')
fprintf('Or have a look into goafem_stochcol_qtychoice.m for other choices \n')
%fprintf('7.x. Variable RHS \n')
% 7.1. Quadratic variable RHS
% 7.2. Exponential variable RHS

%fprintf('8.x. Mommer-Stevenson based Examples \n')
% 8.1. Modified Mommer-Stevenson type problem with added unit L2 part
% 8.2. Variable Mommer-Stevenson (define cut-off variable)
% 9.9. Unsmiley (easteregg)
%fprintf('0.0. More info \n')

QOI_type = default('(default is 1)',1);

switch QOI_type
     case 1 % Rectangular Sub-Domain Integral
        x1_0 = default('Domain x1-coord lower limit (default is 0.25)', 0.25);
        x1_1 = default('Domain x1-coord upper limit (default is 0.75)', 0.75);
        x2_0 = default('Domain x2-coord lower limit (default is 0.25)', 0.25);
        x2_1 = default('Domain x2-coord upper limit (default is 0.75)', 0.75);
        supp = @(x1,x2) 1.*(x1_0 < x1 & x1 < x1_1 & x2_0 < x2 & x2 < x2_1);
        area = (x1_1 - x1_0)*(x2_1 - x2_0);
        G_rhs = {{@(x1,x2) (1/area).*supp(x1,x2), @(x1,x2) zeros(size(x1)), @(x1,x2) zeros(size(x1)), @(x1,x2) zeros(size(x1))}, ...
                 {@(x1,x2) (1/area).*supp(x1,x2), @(x1,x2) zeros(size(x1)), @(x1,x2) zeros(size(x1)), @(x1,x2) zeros(size(x1))}, ...
                 {@(x1,x2) (1/area).*supp(x1,x2), @(x1,x2) zeros(size(x1)), @(x1,x2) zeros(size(x1)), @(x1,x2) zeros(size(x1))}, ...
                 {@(x1,x2) (1/area).*supp(x1,x2), @(x1,x2) zeros(size(x1)), @(x1,x2) zeros(size(x1)), @(x1,x2) zeros(size(x1))}};
    case 1.1 % Unit QOI function
        G_rhs = {{@(x1,x2) ones(size(x1)), @(x1,x2) zeros(size(x1)), @(x1,x2) zeros(size(x1)), @(x1,x2) zeros(size(x1))}, ...
                 {@(x1,x2) ones(size(x1)), @(x1,x2) zeros(size(x1)), @(x1,x2) zeros(size(x1)), @(x1,x2) zeros(size(x1))}, ...
                 {@(x1,x2) ones(size(x1)), @(x1,x2) zeros(size(x1)), @(x1,x2) zeros(size(x1)), @(x1,x2) zeros(size(x1))}, ...
                 {@(x1,x2) ones(size(x1)), @(x1,x2) zeros(size(x1)), @(x1,x2) zeros(size(x1)), @(x1,x2) zeros(size(x1))}};
     
    case 1.2 % Circular Sub-Domain Integral
        x1_0 = default('Integral region centre x1-coord (default is 0.45)', 0.45);
        x2_0 = default('Integral region centre x2-coord (default is 0.65)', 0.65);
        r = default('Integral region radius (default is 0.25)',0.25);
        supp = @(x1,x2) (x1_0 - x1).^2 + (x2_0 - x2).^2 < r.^2;
        area = pi*r^2;
        G_rhs = {{@(x1,x2) (1/area).*supp(x1,x2), @(x1,x2) zeros(size(x1)), @(x1,x2) zeros(size(x1)), @(x1,x2) zeros(size(x1))}, ...
                 {@(x1,x2) (1/area).*supp(x1,x2), @(x1,x2) zeros(size(x1)), @(x1,x2) zeros(size(x1)), @(x1,x2) zeros(size(x1))}, ...
                 {@(x1,x2) (1/area).*supp(x1,x2), @(x1,x2) zeros(size(x1)), @(x1,x2) zeros(size(x1)), @(x1,x2) zeros(size(x1))}, ...
                 {@(x1,x2) (1/area).*supp(x1,x2), @(x1,x2) zeros(size(x1)), @(x1,x2) zeros(size(x1)), @(x1,x2) zeros(size(x1))}};
   
    case 1.3 % Ellipse Sub-Domain Integral
        x1_0 = default('Integral region centre x1-coord (default is 0.4)', 0.4);
        x2_0 = default('Integral region centre x2-coord (default is 0.7)', 0.7);
        x1_1 = default('Horizontal semi-axis (default is 0.1)', 0.1);
        x2_1 = default('Vertical semi-axis (default is 0.3)', 0.3);
        supp = @(x1,x2) ((x1_0 - x1)/x1_1).^2 + ((x2_0 - x2)/x2_1).^2 < 1.^2;
        area = pi*x1_1*x2_1;
        G_rhs = {{@(x1,x2) (1/area).*supp(x1,x2), @(x1,x2) zeros(size(x1)), @(x1,x2) zeros(size(x1)), @(x1,x2) zeros(size(x1))}, ...
                 {@(x1,x2) (1/area).*supp(x1,x2), @(x1,x2) zeros(size(x1)), @(x1,x2) zeros(size(x1)), @(x1,x2) zeros(size(x1))}, ...
                 {@(x1,x2) (1/area).*supp(x1,x2), @(x1,x2) zeros(size(x1)), @(x1,x2) zeros(size(x1)), @(x1,x2) zeros(size(x1))}, ...
                 {@(x1,x2) (1/area).*supp(x1,x2), @(x1,x2) zeros(size(x1)), @(x1,x2) zeros(size(x1)), @(x1,x2) zeros(size(x1))}};
    case 2 % Mollifiers
        F_rhs = {@(x1,x2) ones(size(x1)), @(x1,x2) zeros(size(x1)), @(x1,x2) zeros(size(x1)), @(x1,x2) zeros(size(x1))};
        x1_0 = default('Pointwise estimation: x1-coord (default is 0.4)', 0.4);
        x2_0 = default('Pointwise estimation: x2-coord (default is -0.5)', -0.5);
        r = default('Mollifier Radius (default is 0.15)',0.15);
        Mollconst = 2.1436*r^-2;
        supp = @(x1,x2) (x1_0 - x1).^2 + (x2_0 - x2).^2 < r.^2;
        moll0 = @(x1,x2) exp(-r^2./...
                            (r^2 - ((x1_0 - x1).^2 + (x2_0 - x2).^2)).*((x1_0 - x1).^2 + (x2_0 - x2).^2 < r.^2) + ...
                            1.*((x1_0 - x1).^2 + (x2_0 - x2).^2 >= r.^2)) ;
        moll = @(x1,x2) Mollconst *supp(x1,x2) .* moll0(x1,x2);
        G_rhs = {{@(x1,x2) moll(x1,x2), @(x1,x2) zeros(size(x1)), @(x1,x2) zeros(size(x1)), @(x1,x2) zeros(size(x1))}, ...
                 {@(x1,x2) moll(x1,x2), @(x1,x2) zeros(size(x1)), @(x1,x2) zeros(size(x1)), @(x1,x2) zeros(size(x1))}, ...
                 {@(x1,x2) moll(x1,x2), @(x1,x2) zeros(size(x1)), @(x1,x2) zeros(size(x1)), @(x1,x2) zeros(size(x1))}, ...
                 {@(x1,x2) moll(x1,x2), @(x1,x2) zeros(size(x1)), @(x1,x2) zeros(size(x1)), @(x1,x2) zeros(size(x1))}};
    case 8 % Mommer-Stevenson based example
        G_rhs = {{@(x1,x2) zeros(size(x1)), @(x1,x2) 1.*(x2 >= 1.5 - x1), @(x1,x2) zeros(size(x1)), @(x1,x2) zeros(size(x1))}, ...
                 {@(x1,x2) zeros(size(x1)), @(x1,x2) 1.*(x2 >= 1.0 - x1), @(x1,x2) zeros(size(x1)), @(x1,x2) zeros(size(x1))}, ...
                 {@(x1,x2) zeros(size(x1)), @(x1,x2) 1.*(x2 <= x1 - 1.5), @(x1,x2) zeros(size(x1)), @(x1,x2) zeros(size(x1))}, ...
                 {@(x1,x2) zeros(size(x1)), @(x1,x2) 1.*(x2 <= x1 - 1.5), @(x1,x2) zeros(size(x1)), @(x1,x2) zeros(size(x1))}};
    case 8.1 % Mommer-Stevenson with unit L2
        G_rhs = {{@(x1,x2) ones(size(x1)), @(x1,x2) 1.*(x2 >= 1.5 - x1), @(x1,x2) zeros(size(x1)), @(x1,x2) zeros(size(x1))}, ...
                 {@(x1,x2) ones(size(x1)), @(x1,x2) 1.*(x2 >= 1.0 - x1), @(x1,x2) zeros(size(x1)), @(x1,x2) zeros(size(x1))}, ...
                 {@(x1,x2) ones(size(x1)), @(x1,x2) 1.*(x2 <= x1 - 1.5), @(x1,x2) zeros(size(x1)), @(x1,x2) zeros(size(x1))}, ...
                 {@(x1,x2) ones(size(x1)), @(x1,x2) 1.*(x2 <= x1 - 1.5), @(x1,x2) zeros(size(x1)), @(x1,x2) zeros(size(x1))}};
    case 8.2 % Variable Mommer-Stevenson
        if ~exist('vms','var')
            vms = default('Variable indicator sub-domain size (default is 0.5)', 0.5);
            mag = default('Magnitude of indicator (default is 1)', 1.0);
        end
        G_rhs = {{@(x1,x2) zeros(size(x1)), @(x1,x2) mag.*(x2 >= 2 - vms - x1), @(x1,x2) zeros(size(x1)), @(x1,x2) zeros(size(x1))}, ...
                 {@(x1,x2) zeros(size(x1)), @(x1,x2) mag.*(x2 >= 2 - vms - x1), @(x1,x2) zeros(size(x1)), @(x1,x2) zeros(size(x1))}, ...
                 {@(x1,x2) zeros(size(x1)), @(x1,x2) mag.*(x2 <= x1 - 2 + vms), @(x1,x2) zeros(size(x1)), @(x1,x2) zeros(size(x1))}, ...
                 {@(x1,x2) zeros(size(x1)), @(x1,x2) mag.*(x2 <= x1 - 2 + vms), @(x1,x2) zeros(size(x1)), @(x1,x2) zeros(size(x1))}};
    case 7 
        G_rhs = {{@(x1,x2) ones(size(x1)), @(x1,x2) zeros(size(x1)), @(x1,x2) zeros(size(x1)), @(x1,x2) zeros(size(x1))}, ...
                 {@(x1,x2) ones(size(x1)), @(x1,x2) zeros(size(x1)), @(x1,x2) zeros(size(x1)), @(x1,x2) zeros(size(x1))}, ...
                 {@(x1,x2) ones(size(x1)), @(x1,x2) zeros(size(x1)), @(x1,x2) zeros(size(x1)), @(x1,x2) zeros(size(x1))}, ...
                 {@(x1,x2) ones(size(x1)), @(x1,x2) zeros(size(x1)), @(x1,x2) zeros(size(x1)), @(x1,x2) zeros(size(x1))}};
    
%     case 7 % Linear variable RHS
%         G_rhs = {{@(x1,x2) x2 - x1, @(x1,x2) zeros(size(x1)), @(x1,x2) zeros(size(x1)), @(x1,x2) zeros(size(x1))}, ...
%                  {@(x1,x2) x2 - x1, @(x1,x2) zeros(size(x1)), @(x1,x2) zeros(size(x1)), @(x1,x2) zeros(size(x1))}, ...
%                  {@(x1,x2) x2 - x1, @(x1,x2) zeros(size(x1)), @(x1,x2) zeros(size(x1)), @(x1,x2) zeros(size(x1))}, ...
%                  {@(x1,x2) x2 - x1, @(x1,x2) zeros(size(x1)), @(x1,x2) zeros(size(x1)), @(x1,x2) zeros(size(x1))}};
%     case 7.1 % Quadratic variable RHS
%         G_rhs = {{@(x1,x2) 2.*ones(size(x1)) - x1.^2 - x2.^2, @(x1,x2) zeros(size(x1)), @(x1,x2) zeros(size(x1)), @(x1,x2) zeros(size(x1))}, ...
%                  {@(x1,x2) 2.*ones(size(x1)) - x1.^2 - x2.^2, @(x1,x2) zeros(size(x1)), @(x1,x2) zeros(size(x1)), @(x1,x2) zeros(size(x1))}, ...
%                  {@(x1,x2) 2.*ones(size(x1)) - x1.^2 - x2.^2, @(x1,x2) zeros(size(x1)), @(x1,x2) zeros(size(x1)), @(x1,x2) zeros(size(x1))}, ...
%                  {@(x1,x2) 2.*ones(size(x1)) - x1.^2 - x2.^2, @(x1,x2) zeros(size(x1)), @(x1,x2) zeros(size(x1)), @(x1,x2) zeros(size(x1))}};
%     case 7.2 % Exponential variable RHS
%         G_rhs = {{@(x1,x2) exp(-((x1 + 0.5).^2 + (x2 - 0.5).^2)), @(x1,x2) zeros(size(x1)), @(x1,x2) zeros(size(x1)), @(x1,x2) zeros(size(x1))}, ...
%                  {@(x1,x2) exp(-((x1 + 0.5).^2 + (x2 - 0.5).^2)), @(x1,x2) zeros(size(x1)), @(x1,x2) zeros(size(x1)), @(x1,x2) zeros(size(x1))}, ...
%                  {@(x1,x2) exp(-((x1 + 0.5).^2 + (x2 - 0.5).^2)), @(x1,x2) zeros(size(x1)), @(x1,x2) zeros(size(x1)), @(x1,x2) zeros(size(x1))}, ...
%                  {@(x1,x2) exp(-((x1 + 0.5).^2 + (x2 - 0.5).^2)), @(x1,x2) zeros(size(x1)), @(x1,x2) zeros(size(x1)), @(x1,x2) zeros(size(x1))}};             
    case 5 % Non-linear (second moment of true solution in sub-domain)
        weight_goal = @(x1, x2) 4.*(0.25 < x1 & x1 < 0.75 & 0.25 < x2 & x2 < 0.75);
        G_rhs = {{weight_goal, @(x1,x2) zeros(size(x1)), @(x1,x2) zeros(size(x1)), @(x1,x2) zeros(size(x1))}, ...
                 {weight_goal, @(x1,x2) zeros(size(x1)), @(x1,x2) zeros(size(x1)), @(x1,x2) zeros(size(x1))}, ...
                 {weight_goal, @(x1,x2) zeros(size(x1)), @(x1,x2) zeros(size(x1)), @(x1,x2) zeros(size(x1))}, ...
                 {weight_goal, @(x1,x2) zeros(size(x1)), @(x1,x2) zeros(size(x1)), @(x1,x2) zeros(size(x1))}};
        %
        
    case 4 % Non-linear (square of convection term in sub-domain) %%% TODO
        %weight_goal = @(x1,x2) 2.*(x2 >= 1 - x1);
        G_rhs = {{@(x1,x2) 2.*(x2 >= 1 - x1), @(x1,x2) zeros(size(x1)), @(x1,x2) zeros(size(x1)), @(x1,x2) zeros(size(x1))}, ...
                 {@(x1,x2) 2.*(x2 >= 1 - x1), @(x1,x2) zeros(size(x1)), @(x1,x2) zeros(size(x1)), @(x1,x2) zeros(size(x1))}, ...
                 {@(x1,x2) 2.*(x2 >= 1 - x1), @(x1,x2) zeros(size(x1)), @(x1,x2) zeros(size(x1)), @(x1,x2) zeros(size(x1))}, ...
                 {@(x1,x2) 2.*(x2 >= 1 - x1), @(x1,x2) zeros(size(x1)), @(x1,x2) zeros(size(x1)), @(x1,x2) zeros(size(x1))}};
             
    case 3 % Non-linear (square of integral in sub-domain)
        %weight_goal = @(x1,x2) 8.*(x2 >= 1.5 - x1);
        G_rhs = {{@(x1,x2) 8.*(x2 >= 1.5 - x1), @(x1,x2) zeros(size(x1)), @(x1,x2) zeros(size(x1)), @(x1,x2) zeros(size(x1))}, ...
                 {@(x1,x2) 8.*(x2 >= 1.5 - x1), @(x1,x2) zeros(size(x1)), @(x1,x2) zeros(size(x1)), @(x1,x2) zeros(size(x1))}, ...
                 {@(x1,x2) 8.*(x2 >= 1.5 - x1), @(x1,x2) zeros(size(x1)), @(x1,x2) zeros(size(x1)), @(x1,x2) zeros(size(x1))}, ...
                 {@(x1,x2) 8.*(x2 >= 1.5 - x1), @(x1,x2) zeros(size(x1)), @(x1,x2) zeros(size(x1)), @(x1,x2) zeros(size(x1))}};

    case 6 %Variance on linear functional + Molifier weight
        x1_0 = default('Pointwise estimation: x1-coord (default is 0.4)', 0.4);
        x2_0 = default('Pointwise estimation: x2-coord (default is -0.5)', -0.5);
        r = default('Mollifier Radius (default is 0.15)',0.15);
        E = 0.1511;%default('Mean value of int_D w(x) u(x,y). Can be found by running the same problem with linear functional int_\Gamma int_D w(x) u(x,y) (default is 1)', 1);
        Mollconst = 2.1436*r^-2;
        supp = @(x1,x2) (x1_0 - x1).^2 + (x2_0 - x2).^2 < r.^2;
        moll0 = @(x1,x2) exp(-r^2./...
                            (r^2 - ((x1_0 - x1).^2 + (x2_0 - x2).^2)).*((x1_0 - x1).^2 + (x2_0 - x2).^2 < r.^2) + ...
                            1.*((x1_0 - x1).^2 + (x2_0 - x2).^2 >= r.^2)) ;
        moll = @(x1,x2) Mollconst *supp(x1,x2) .* moll0(x1,x2);
        G_rhs = {{@(x1,x2) moll(x1,x2), @(x1,x2) zeros(size(x1)), @(x1,x2) zeros(size(x1)), @(x1,x2) zeros(size(x1))}, ...
                 {@(x1,x2) moll(x1,x2), @(x1,x2) zeros(size(x1)), @(x1,x2) zeros(size(x1)), @(x1,x2) zeros(size(x1))}, ...
                 {@(x1,x2) moll(x1,x2), @(x1,x2) zeros(size(x1)), @(x1,x2) zeros(size(x1)), @(x1,x2) zeros(size(x1))}, ...
                 {@(x1,x2) moll(x1,x2), @(x1,x2) zeros(size(x1)), @(x1,x2) zeros(size(x1)), @(x1,x2) zeros(size(x1))}};
             
    case 9.9 % Unsmiley (easteregg)
        FACE_q = {@(x1,x2) 1.*(0.32-2.15.*(0.5-x1).^2 < x2 & 0.4-3.*(0.5-x1).^2 > x2), ...
                  @(x1,x2) 1.*(0.32-2.15.*(0.5*x1).^2 < 0.5*x2 + 0.5 & 0.4-3.*(0.5*x1).^2 > 0.5*x2 + 0.5), ...
                  @(x1,x2) 1.*((x1-0.6).^2+(x2-0.1).^2 < 0.03), ...
                  @(x1,x2) 1.*(0.32-2.15.*(0.5*x1).^2 < 0.5*x2 + 0.5 & 0.4-3.*(0.5*x1).^2 > 0.5*x2 + 0.5)};
        if dom_type == 1
            FACE_q = FACE_q{1};
        elseif dom_type == 5
            FACE_q = FACE_q{2};
        elseif dom_type == 2
            FACE_q = FACE_q{3}; 
        elseif dom_type == 4
            FACE_q = FACE_q{4};
        else
            error('For this domain the procedure is unavailable');
        end
        G_rhs = {{@(x1,x2) FACE_q(x1,x2), @(x1,x2) zeros(size(x1)), @(x1,x2) zeros(size(x1)), @(x1,x2) zeros(size(x1))}, ...
                 {@(x1,x2) FACE_q(x1,x2), @(x1,x2) zeros(size(x1)), @(x1,x2) zeros(size(x1)), @(x1,x2) zeros(size(x1))}, ...
                 {@(x1,x2) FACE_q(x1,x2), @(x1,x2) zeros(size(x1)), @(x1,x2) zeros(size(x1)), @(x1,x2) zeros(size(x1))}, ...
                 {@(x1,x2) FACE_q(x1,x2), @(x1,x2) zeros(size(x1)), @(x1,x2) zeros(size(x1)), @(x1,x2) zeros(size(x1))}};
end

if dom_type == 1
    G_rhs = G_rhs{1};
elseif dom_type == 5
    G_rhs = G_rhs{2};
elseif dom_type == 2
    G_rhs = G_rhs{3}; 
elseif dom_type == 4
    G_rhs = G_rhs{4};
elseif dom_type == 3
    G_rhs = G_rhs{1};
else
    error('For this domain the procedure is unavailable');
end
G_rhs{5} = -1;

if QOI_type == 5 || QOI_type == 7
    G_rhs{5} = 0;
elseif QOI_type == 4
    G_rhs{5} = 1;
elseif QOI_type == 3
    G_rhs{5} = 2;
elseif QOI_type == 6
    G_rhs{5} = 3;
    G_rhs{6} = E;
end
