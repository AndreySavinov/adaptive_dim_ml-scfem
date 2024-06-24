function plotstats(M,paras_sg,coords,fun_p,paras_fem,sols, dom_type)
% postprocess to generate statistics

multiweights = paras_sg{8};
MEAN = 0;
VAR  = 0;
for k = 1:size(coords,1)
    % total pdf
    pyy = 1;
    for m = 1:M
        pyy = pyy * fun_p(coords(k,m));
    end
    MEAN = MEAN + pyy * multiweights(k) * sols(:,k); % mean
end

for k = 1:size(coords,1)
    % total pdf
    pyy = 1;
    for m = 1:M
        pyy = pyy * fun_p(coords(k,m));
    end
    VAR = VAR + pyy * multiweights(k) * (sols(:,k)-MEAN).^2; % mean
end
STD = VAR.^(1/2);
sc_plotsol_p1(dom_type, MEAN, VAR,paras_fem{2},paras_fem{1})
end