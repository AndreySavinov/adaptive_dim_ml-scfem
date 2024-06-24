function G = goafem_stochcol_matrix1(grid,clcset,indlcset,list,m,listy)
% GOAFEM_STOCHCOL_MATRIX calculates modified integrals of Lagrange basis
% functions
%
%  G = goafem_stochcol_matrix1(grid,clcset,indlcset,list)
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
%                 G   modified integrals values
%
%
% SEE ALSO: stochcol_gmatrices
%           goafem_stochcol_uni_int
%           goafem_doBilinearForm
%
% TR; 24 November 2022

K = size(grid, 1);
N = size(grid, 2);
G = zeros(K, 1);
for k = 1:K
    temp = 0;
    cset_k = clcset{k};
    indset_k = indlcset{k};
    for i = 1:length(cset_k)
        temp_int = 1;
        for n = 1:N
            if n == m
                temp_int = temp_int*listy{indset_k(i,n),1}(grid(k,n),1);
            else % n =/= m
                temp_int = temp_int*list{indset_k(i,n),1}(grid(k,n),1);
            end
        end
        temp = temp + temp_int*cset_k(i);
    end
    G(k) = temp;
end
    
end