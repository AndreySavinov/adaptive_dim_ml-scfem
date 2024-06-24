function spind = goafem_indicator(su, sz, pu, pz, rule)
% GOAFEM_INDICATOR computes an indicator to determine whether the next
% iteration in the adaptive loop should use spatial or parametric
% refinement.
%
% spind = goafem_indicator(sx, sz, px, pz, rule)
%
% input:
%         su    primal spatial estimate
%         sz    dual spatial estimate
%         pu    primal parametric estimate
%         pz    dual parametric estimate
%       rule    specific rule for selection
%       
% output:
%      spind    indicator to determine spatial(0)/parametric(1) refinement
%
% SEE ALSO: goafem_singlelevelSC
%
% TR; 17 july 2022

if rule == 1 % determine by combined product
    spind = su * sz < pu * pz;
elseif rule == 2 % determine by largest quantity
    spind = prod([su, sz, pu] <= pz) + prod([su, sz, pz] <= pu);
    % Eliminates any ties
    spind = any(spind);
elseif rule == 3
    spind = su^2 + sz^2 < pu^2 + pz^2;
end
end

