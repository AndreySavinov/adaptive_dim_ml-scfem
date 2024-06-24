function [Mset1, Mset2] = goafem_mark_combine(Mu, Mz, mel, rule, ...
                          eu, ez, eu2, ez2, paras_detail)
% GOAFEM_MARK_COMBINE combines primal and dual marked sets for both spatial
% and parametric cases. Different rules for merging marked sets are
% provided.
%
% [Mset1, Mset2] = goafem_mark_combine(Mu, Mz, mel, rule, ...
%                  eu, ez, eu2, ez2, paras_detail)
% 
% input:
%               Mu    primal marked set
%               Mz    dual marked set
%              mel    switch for marking elements(1)/edges(2)
%             rule    rule for merging marked sets
%               eu    primal (local) error indicators
%               ez    dual (local) error indicators
%              eu2    second primal error indicators (spatial only)
%              ez2    second dual error indicators (spatial only)
%     paras_detail    detail space parameters (spatial only)
%             
% output:
%            Mset1    a marked set for refinement
%            Mset2    a marked set for refinement (spatial only)
%
% SEE ALSO: goafem_singlelevelSC
%
% TR; 13 July 2022

% Spatial refinement
if nargin == 9 
    evtY    = paras_detail{4};
    Ybasis  = paras_detail{7};
% Rule 1 - simple union of primal/dual sets
    if rule == 1
        Mset = union(Mu, Mz, 'rows');
% Rule 2 - union of smaller set and subset of larger set
    elseif rule == 2 % 
        % Switch to edge indicators if edge-based marking is used
        if mel == 2
            eu = eu2;
            ez = ez2;
        end
        % Define larger/smaller sets
        if size(Mu,1) <= size(Mz,1)
            M1 = Mu;
            M2 = Mz;
            e = ez(M2);
        else % size(Mz,1) < size(Mx,1)
            M1 = Mz;
            M2 = Mu;
            e = eu(M2);
        end
        % Find largest indicators from the larger set
        M2 = sortrows([M2 e],size([M2 e],2),'descend');
        M2 = M2(1:size(M1,1),1:size(M1,2));
        % Initial set of edges/elements (depends on markedgelem)
        Mset = union(M1, M2, 'rows');
    elseif rule == 3
        if size(Mu,1) < size(Mz,1)
            Mset = Mu;
        else
            Mset = Mz;
        end
    end
    % Overall set of marked elements and edges (fixes hanging nodes etc.)
    % Mset1 = Elements, Mset2 = Edges
    [Mset1, Mset2] = get_all_marked_elem(Ybasis, evtY, Mset, mel);
   
% Parametric refinement
elseif nargin == 6 
% Rule 1 - simple union of primal/dual sets
    if rule == 1 
        Mset1 = union(Mu, Mz, 'rows');
% Rule 2 - union of smaller set and subset of larger set
    elseif rule == 2 
        % Define larger/smaller sets
        if size(Mu,1) <= size(Mz,1)
            M1 = Mu;
            M2 = Mz;
            e = ez;
        else % size(Mz,1) < size(Mx,1)
            M1 = Mz;
            M2 = Mu;
            e = eu;
        end
        % Find largest indicators from the larger set
        M2 = sortrows([M2 e],size([M2 e],2),'descend');
        M2 = M2(1:size(M1,1),1:size(M1,2));
        % Final set of marked collocation points
        % (Note: If the reduced margin is not used in the main driver, this
        %        may need to be modified to force downward-closedness.)
        Mset1 = union(M1, M2, 'rows');
    elseif rule == 3
        if size(Mu,1) < size(Mz,1)
            Mset1 = Mu;
        else
            Mset1 = Mz;
        end
    end
    
end