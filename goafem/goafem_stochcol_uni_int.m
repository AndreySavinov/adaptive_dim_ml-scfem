function [list, listy, listy2, listy3, listy4] = ...
    goafem_stochcol_uni_int(fun_p, polys, L)
% GOAFEM_STOCHCOL_UNI_INT calculates a cell containing integrals of
% products of Lagrange basis functions.
%
%   input:
%             fun_p   function handle for the probability density function
%             polys   one-dimensional Lagrange polynomials
%                 L   integral limit (parameter space is of form [-L, L]^M)
%
%   output:
%              list   integral cell array of products of Lagrange functions
%             listy   modified cell array associated with linear expansions of diffusion coeffients
%            listy2   modified cell array associated with quadratic expansions of diffusion coeffients
%            listy3   modified cell array associated with expansions of diffusion coeffients
%            listy4   modified cell array associated with expansions of diffusion coeffients
%
% Typically, L = 1.
% 
% In each of the outputs, a given cell in the array will contain each
% integral for every combination of integrals in that collocation point
% level-pair.
%
% Function(s) called: goafem_stochcol_fem_setup
%                     goafem_stochcol_fem_solver
%
% SEE ALSO: stochcol_onedlagpolys
%           stochcol_uni_int
%
% TR; 28 September 2022

startpreTime = tic;

level = length(polys);
list = cell(level, level);
listy = cell(level, level);
listy2 = cell(level, level);
listy3 = cell(level, level);
listy4 = cell(level, level);

progress = waitbar(0, 'Generating precomputation data...');
counter = 0;
tot = 0;

for i = 1:level
    for j = 1:i
        M = length(polys{i});
        N = length(polys{j});
        tot = tot + M*N;
    end
end
        
for i = 1:level
    for j = 1:i
        M = length(polys{i});
        N = length(polys{j});
        list{i,j} = zeros(M, N);
        listy{i,j} = zeros(M, N);
        listy2{i,j} = zeros(M, N);
        for m = 1:M
            for n = 1:N
                integrand = @(x) fun_p(x).* polys{i}{m}(x).*polys{j}{n}(x);
                integrand1 = @(x) x.* integrand(x);
                integrand2 = @(x) x.* x.* integrand(x);            
                list{i,j}(m,n) = integral(integrand,-L,L,'ArrayValued',true);
                listy{i,j}(m,n) = integral(integrand1,-L,L,'ArrayValued',true);
                listy2{i,j}(m,n) = integral(integrand2,-L,L,'ArrayValued',true);                
                counter = counter + 1;
                try 
                    waitbar(counter/(2*tot),progress,'Generating polynomial data...') 
                catch
                end
            end
        end
        listy{i,j}(abs(listy{i,j}) <= 1e-16) = 0.0;
        listy2{i,j}(abs(listy2{i,j}) <= 1e-16) = 0.0;
        list{j,i} = transpose(list{i,j});
        listy{j,i} = transpose(listy{i,j});
        listy2{j,i} = transpose(listy2{i,j});
        listy3{i,j} = zeros(M, N);
        listy4{i,j} = zeros(M, N);
        for m = 1:M
            for n = 1:N
                integrand3 = @(x) x.* x.* x.* integrand(x);
                integrand4 = @(x) x.* x.* x.* x.* integrand(x);
                listy3{i,j}(m,n) = integral(integrand3,-L,L,'ArrayValued',true);
                listy4{i,j}(m,n) = integral(integrand4,-L,L,'ArrayValued',true);
                counter = counter + 1;
                try 
                    waitbar(counter/(2*tot),progress,'Generating polynomial data...')
                catch
                end
            end
        end
        listy3{i,j}(abs(listy3{i,j}) <= 1e-16) = 0.0;
        listy4{i,j}(abs(listy4{i,j}) <= 1e-16) = 0.0;
        listy3{j,i} = transpose(listy3{i,j});
        listy4{j,i} = transpose(listy4{i,j});
    end
end

close(progress)
endpreTime = toc(startpreTime);
fprintf('\n<strong>Precomputation time:</strong> %.2f sec\n',endpreTime);

end