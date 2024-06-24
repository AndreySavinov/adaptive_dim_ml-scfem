function [rxy,revt,rbound,rinterior,reboundt] = uniform_refinement(xy,evt,bound,eboundt,refine_type,iplot)
%UNIFORM_REFINEMENT uniformly refines the triangular mesh using either red or bisec3 refinement
%
% [rxy,revt,rbound,rinterior,reboundt] = uniform_refinement(xy,evt,bound,...
%                                        eboundt,refine_type,iplot)
%
% input:
%               xy     vertex coordinate vector  
%              evt     element vertex mapping matrix
%            bound     boundary vertex vector
%          eboundt     element boundary mapping matrix
%      refine_type     red or bisec3 refinement type
%            iplot     (optional) grid plotting switch
%
% output:
%              rxy     refined vertex-coordinate vector
%             revt     refined element-vertex-mapping matrix
%           rbound     refined boundary-vertex vector
%        rinterior     refined interior-vertex vector
%         reboundt     refined element-boundary-mapping matrix
%
% Function(s) called:  p2_grid_generator
%                      plot_mesh
%
%   TIFISS function: LR; 22 June 2018
% Copyright (c) 2018 A. Bespalov, L. Rocchi

  if nargin < 6
      iplot = 0;
  end
    
% Calling the p2 grid generator  
  [p2xy,p2evt,p2bound] = p2_grid_generator(xy,evt,bound);
 
% Save new rxy and rbound matrices
  rxy    = p2xy;
  rbound = p2bound;

% Allocate memory 
  nel  = size(evt,1);      % number of "old" elements
  revt = zeros(4*nel,3);   % number of new elements in the refined mesh
  
  if refine_type == 1
      % -----------------------------------------------------------
      % Uniform RED refinement
      % -----------------------------------------------------------
      revt(1:4:(4*nel),:) = p2evt(:,[1,6,5]);
      revt(2:4:(4*nel),:) = p2evt(:,[6,2,4]);  
      revt(3:4:(4*nel),:) = p2evt(:,[5,4,3]);
      revt(4:4:(4*nel),:) = p2evt(:,[4,5,6]);
      %
      % Data structures for the refined reboundt matrix
      %     
      % Extract boundary elements according to the boundary edge's number
      bedgeOne   = eboundt(eboundt(:,2) == 1,:);
      bedgeTwo   = eboundt(eboundt(:,2) == 2,:);
      bedgeThree = eboundt(eboundt(:,2) == 3,:);
      %
      % For a given (old) boundary element in the eboundt matrix, 
      % there will be only 2 (out of the 4) new boundary element-children: 
      % according to the old element's edge number, 1, 2, or 3, the 
      % boundary children will be respectively the [2nd, 3rd], [1st, 3rd] 
      % or [1st 2nd] elements created
      matElOne   = [4*(bedgeOne(:,1)-1),   4*(bedgeOne(:,1)-1)]    +  repmat([2 3],size(bedgeOne,1),1);
      matElTwo   = [4*(bedgeTwo(:,1)-1),   4*(bedgeTwo(:,1)-1)]    +  repmat([1 3],size(bedgeTwo,1),1);
      matElThree = [4*(bedgeThree(:,1)-1), 4*(bedgeThree(:,1)-1)]  +  repmat([1 2],size(bedgeThree,1),1);
      %
      % Adding the corresponding edge'number per each child boundary
      % element: according to the old element's edge number, 1, 2, or 3, 
      % the boundary edges will be [1 1], [2 2], or [3 3] 
      matElOne   = [matElOne(:,1),   ones(size(matElOne,1),1),       matElOne(:,2),   ones(size(matElOne,1),1)];
      matElTwo   = [matElTwo(:,1),   repmat(2,size(matElTwo,1),1),   matElTwo(:,2),   repmat(2,size(matElTwo,1),1)];
      matElThree = [matElThree(:,1), repmat(3,size(matElThree,1),1), matElThree(:,2), repmat(3,size(matElThree,1),1)];
      
  else% refine_type == 2     
      % -----------------------------------------------------------
      % Uniform BISEC3 refinement
      % -----------------------------------------------------------
      revt(1:4:(4*nel),:) = p2evt(:,[1,6,5]);
      revt(2:4:(4*nel),:) = p2evt(:,[5,6,2]);
      revt(3:4:(4*nel),:) = p2evt(:,[2,4,5]);
      revt(4:4:(4*nel),:) = p2evt(:,[5,4,3]);
      %
      % Data structures for the refined reboundt matrix
      %      
      % Extract boundary elements according to the boundary edge's number
      bedgeOne   = eboundt(eboundt(:,2) == 1,:);
      bedgeTwo   = eboundt(eboundt(:,2) == 2,:);
      bedgeThree = eboundt(eboundt(:,2) == 3,:);
      %
      % For a given (old) boundary element in the eboundt matrix, 
      % there will be only 2 (out of the 4) new boundary element-children: 
      % according to the old element's edge number, 1, 2, or 3, the 
      % boundary children will be respectively the [3rd 4th], [1st, 4th] 
      % or [1st 2nd] elements created
      matElOne   = [4*(bedgeOne(:,1)-1),   4*(bedgeOne(:,1)-1)]    +  repmat([3 4],size(bedgeOne,1),1);
      matElTwo   = [4*(bedgeTwo(:,1)-1),   4*(bedgeTwo(:,1)-1)]    +  repmat([1 4],size(bedgeTwo,1),1);
      matElThree = [4*(bedgeThree(:,1)-1), 4*(bedgeThree(:,1)-1)]  +  repmat([1 2],size(bedgeThree,1),1);
      %
      % Adding the corresponding edge'number per each child boundary
      % element: according to the old element's edge number, 1, 2, or 3, 
      % the boundary edges will be [3 1], [2 2], or [3 1] 
      matElOne   = [matElOne(:,1),   repmat(3,size(matElOne,1),1),   matElOne(:,2),   ones(size(matElOne,1),1)];
      matElTwo   = [matElTwo(:,1),   repmat(2,size(matElTwo,1),1),   matElTwo(:,2),   repmat(2,size(matElTwo,1),1)];
      matElThree = [matElThree(:,1), repmat(3,size(matElThree,1),1), matElThree(:,2), ones(size(matElThree,1),1)];
          
  end
  
% --------------------------------------------------------------------------
% Creating the reboundt matrix
% --------------------------------------------------------------------------
  reboundt = [ matElOne(:,[1,2]);   matElOne(:,[3,4]); ...
               matElTwo(:,[1,2]);   matElTwo(:,[3,4]); ...
               matElThree(:,[1,2]); matElThree(:,[3,4] )];               
% NOTE that reboundt does not need to be sorted. If sorting is wanted, 
% sortrows has to be called, but this will slow down (a little) the computation
% reboundt = sortrows(reboundt,1); 
  
% Update the interior nodes
  totalnodes = 1:size(rxy,1);
  rinterior = totalnodes(~ismember(totalnodes,rbound));
  
  if iplot == 1
      % Plot the uniform refined mesh
      plot_mesh(revt,p2xy,'Uniform refined mesh');
  end

end % end function