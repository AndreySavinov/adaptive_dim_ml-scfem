function G = goafem_stochcol_matrix2(grid,clcset,indlcset,list,m,m2,listy,listy2)
%GOAFEM_STOCHCOL_MATRIX single matrix integral calculator
%
%  G = goafem_stochcol_matrix(grid,clcset,indlcset,list)
%
%
%

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
            if n == m || n == m2
                if m == m2
                    temp_int = temp_int*listy2{indset_k(i,n),1}(grid(k,n),1);
                else % m =/= m2
                    temp_int = temp_int*listy{indset_k(i,n),1}(grid(k,n),1);
                end
            else % n =/= m
                temp_int = temp_int*list{indset_k(i,n),1}(grid(k,n),1);
            end
        end
        temp = temp + temp_int*cset_k(i);
    end
    G(k) = temp;
end
    
end