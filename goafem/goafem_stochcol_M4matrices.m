function G = goafem_stochcol_M4matrices(grid,clcset,indlcset,list,m,m2,m3,m4,listy,listy2,listy3,listy4)
% GOAFEM_STOCHCOL_M2MATRICES calculates modified integrals of Lagrange basis
% functions associated with quartic expansions of diffusion coefficients
%
%   input:
%              grid   sparse grid parameter
%            clcset   sparse grid parameter
%          indlcset   sparse grid parameter
%              list   integral cell array of products of Lagrange functions
%                 m   First parameter direction for which the integral is modified
%                m2   Second parameter direction for which the integral is modified
%                m3   Third parameter direction for which the integral is modified
%                m4   Fourth parameter direction for which the integral is modified
%             listy   modified cell array associated with (linear) expansions of diffusion coeffients
%            listy2   modified cell array associated with (quadratic) expansions of diffusion coeffients
%            listy3   modified cell array associated with (cubic) expansions of diffusion coeffients
%            listy4   modified cell array associated with (quartic) expansions of diffusion coeffients
%
%   output:
%                 G   modified integrals corresponding to quadratic expansions
%
% G = goafem_stochcol_Mmatrices(grid,clcset,indlcset,list,m,n,listy,listy2)
%
% SEE ALSO: stochcol_gmatrices
%           goafem_stochcol_Mmatrices
%           goafem_stochcol_uni_int
%           goafem_doBilinearForm
%
% TR; 18 November 2022
%
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
                    if n == m || n == m2 || n == m3 || n == m4
                        if isequal(m,m2,m3,m4) % equal to all four
                            temp_int = temp_int*listy4{indset_k(i,n), ...
                            indset_l(j,n)}(grid(k,n),grid(l,n));
                        elseif isequal(n,m,m2,m3) || isequal(n,m2,m3,m4) || isequal (n,m3,m4,m) || isequal(n,m4,m,m2) % equal to three
                            temp_int = temp_int*listy3{indset_k(i,n), ...
                            indset_l(j,n)}(grid(k,n),grid(l,n));
                        elseif isequal(n,m,m2) || isequal(n,m2,m3) || isequal (n,m3,m4) || isequal(n,m4,m) || isequal(n,m,m3) || isequal(n,m2,m4) % equal to two
                            temp_int = temp_int*listy2{indset_k(i,n), ...
                            indset_l(j,n)}(grid(k,n),grid(l,n));
                        else % equal to one
                            temp_int = temp_int*listy{indset_k(i,n), ...
                            indset_l(j,n)}(grid(k,n),grid(l,n));
                        end
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