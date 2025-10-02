function list = stochcol_uni_int(fun_p,polys,L)
% STOCHCOL_UNI_INT (pre)computes integrals of products of basis functions

level = length(polys);
list = cell(level, level);
for i = 1:level
    for j = 1:i
        M = length(polys{i});
        N = length(polys{j});
        list{i,j} = zeros(M, N);
        for m = 1:M
            for n = 1:N
                integrand = @(x) fun_p(x).*polys{i}{m}(x).*polys{j}{n}(x);
                for L_current = -L:L-1
                    list{i,j}(m,n) = list{i,j}(m,n) + integral(integrand,L_current,L_current+1,'ArrayValued',true, 'RelTol', 1e-20);
                end
                if abs(list{i,j}(m,n)) < 1e-10
                    list{i,j}(m,n) = 0;
                end
            end
        end
        list{j,i} = transpose(list{i,j});
    end
end
