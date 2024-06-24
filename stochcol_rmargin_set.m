function Y = stochcol_rmargin_set(X, XMargin, Q)
%STOCHCOL_RMARGIN_SET compute the reduced margin set of the downward closed set
%  Y = stochcol_rmargin_set(X, XMargin)
%   inputs: X         index set
%           XMargin   margin set of X
%
%   TIFISS function: FX 10 October 2019
% Copyright (c) 2019 F. Xu

% prepare X
X = sortrows(X);
% test whether X is downward closed;
Xdc = stochcol_getgrid(X);
if ~isequal(X, Xdc)
    error('The input is not a downward closed set!')
end

[r, c] = size(XMargin);
UnitSet = eye(c);
Y = [];
for i = 1:r % test each index in margin set
    index = XMargin(i,:);
    diffset = ones(c, 1)*index - UnitSet;
    diffset_new = [];
    for j = 1:c
        if min(diffset(j,:))~=0
            diffset_new = [diffset_new;diffset(j,:)];
        end
    end
    if ismember(diffset_new, X, 'rows')
        Y = unique([Y; index], 'rows');
    end
end
% candidates = [ones(Q+1, size(Y,2)), ones(Q+1, Q+1)+flip(eye(Q+1))];
% count_cand = 0;
% for i = 1:Q+1
%     if ismember(diffset_new, X, 'rows')
Y = [[ones(Q, size(Y,2)), ones(Q, Q)+flip(eye(Q))];...
    Y, ones(size(Y,1), Q);
    ];
    
    
    
    