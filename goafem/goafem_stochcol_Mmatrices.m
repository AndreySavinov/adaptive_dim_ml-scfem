function G = goafem_stochcol_Mmatrices(grid,clcset,indlcset,list,m,listy)
% GOAFEM_STOCHCOL_MMATRICES calculates modified integrals of Lagrange basis
% functions associated with linear expansions of diffusion coefficients
%
% G = goafem_stochcol_Mmatrices(grid,clcset,indlcset,list,m,listy)
%
%   input:
%              grid   sparse grid parameter
%            clcset   sparse grid parameter
%          indlcset   sparse grid parameter
%              list   integral cell array of products of Lagrange functions
%                 m   Parameter direction for which the integral is modified
%             listy   modified cell array associated with (linear) expansions of diffusion coeffients
%
%   output:
%                 G   modified integrals corresponding to linear expansions
%

%
% SEE ALSO: stochcol_gmatrices
%           goafem_stochcol_uni_int
%           goafem_doBilinearForm
%
% TR; 12 August 2022
K = size(grid, 1);
N = size(grid, 2);
G = zeros(K, K);
for k = 1:K
    for l = 1:k
        temp = 0;
        cset_k = clcset{k};
        cset_l = clcset{l};
        indset_k = indlcset{k};
        indset_l = indlcset{l};
        for i = 1:length(cset_k)
            for j = 1:length(cset_l)
                temp_int = 1;
                for n = 1:N
                    if n == m
                        temp_int = temp_int*listy{indset_k(i,n), ...
                        indset_l(j,n)}(grid(k,n),grid(l,n));
                    else % n =/= m
                        temp_int = temp_int*list{indset_k(i,n), ...
                        indset_l(j,n)}(grid(k,n),grid(l,n));
                    end
                end
                temp = temp + temp_int*cset_k(i)*cset_l(j);
            end
        end
        G(k,l) = temp;
        G(l,k) = temp;
    end
end