function goafem_sc_plotsol_p1(dom_type,sol_u,sol_z,var_u,var_z,evt,xy)
% GOAFEM_SC_PLOTSOL_P1 plots expectation and variance of both primal and
% dual solutions.
% input:
%          dom_type   type of spatial domain
%             sol_u   primal solution
%             sol_z   dual solution
%             var_u   variance of primal solution
%             var_z   variance of dual solution
%               evt   indexing of mesh points
%                xy   coordinates of mesh vertices
%
% SEE ALSO: goafem_singlelevelSC
% NOTE that the mean-field and the variance are interpolated on a square grid 
% [X,Y] in order to plot isolines (contour plot); Matlab does not provide a 
% countour function for mesh-based functions.
%
% TR; 28 September 2022

    
  nvtx = size(xy,1);    % Number of vertices
 
% Refine grid and get cartesian product mesh
  npoints = 150; 
  x = linspace(min(xy(:,1)),max(xy(:,1)),npoints);
  y = linspace(min(xy(:,2)),max(xy(:,2)),npoints); 
  [X,Y] = meshgrid(x,y);
  
% -----------------------------------------------------------------------------  
% Expectation and variance of the solution
% -----------------------------------------------------------------------------
  expsol_u = griddata(xy(:,1),xy(:,2),sol_u(1:nvtx),X,Y);
  expsol_z = griddata(xy(:,1),xy(:,2),sol_z(1:nvtx),X,Y);
  varsol_u = griddata(xy(:,1),xy(:,2),var_u(1:nvtx),X,Y);
  varsol_z = griddata(xy(:,1),xy(:,2),var_z(1:nvtx),X,Y);
% Fix to zero (eventual very small) negative values of the variance
% (due to griddata interpolation)

        
% Get rid of points outside the domain for X-Y grid
  [expsol_u,varsol_u] = recover_domain(dom_type,expsol_u,varsol_u,X,Y);  
  [expsol_z,varsol_z] = recover_domain(dom_type,expsol_z,varsol_z,X,Y);  
 
% Plot mean value and variance
  title1_x = 'Expectation of the Primal Solution';
  title2_x = 'Variance of the Primal Solution';
  plot_stuff(X,Y,expsol_u,varsol_u,sol_u(1:nvtx),var_u,xy,evt,dom_type,title1_x,title2_x,8);
  
  title1_z = 'Expectation of the Dual Solution';
  title2_z = 'Variance of the Dual Solution';
  plot_stuff(X,Y,expsol_z,varsol_z,sol_z(1:nvtx),var_z,xy,evt,dom_type,title1_z,title2_z,9);
  
 end % end function  


% -----------------------------------------------------------------------------
% Child function
% -----------------------------------------------------------------------------
function [func1,func2] = recover_domain(dom_type,func1,func2,X,Y)
% Get rid of those points of the cartesian grid X-Y outside the domain
  
% Nothing to do for convex domains (square):
% - dom_type == 1 || dom_type == 2 
  if dom_type == 2
      % L-shaped domain
      func1(X<0 & Y<0) = nan;
      func2(X<0 & Y<0) = nan;
  end
  
end % end child function


% -----------------------------------------------------------------------------
% Child function
% -----------------------------------------------------------------------------
function outline_domain(dom_type)
% Outline the corresponding domain 
  if dom_type == 1
      unitsquare;
  elseif dom_type == 2
      ellx;
  elseif dom_type == 3
      large_square;
  elseif dom_type == 4
      largecrack;
  elseif dom_type == 5
      largesquare;
  end
end % end child function


% -----------------------------------------------------------------------------
% Child function
% -----------------------------------------------------------------------------
function plot_stuff(X,Y,func1,func2,expmesh,varmesh,xy,evt,dom_type,title1,title2,figno)
% In this plot function both mean value and variance are plotted using trimesh 

  fontSize = 14;

% Plot functions 
  figure(figno);
  
  subplot(221)
  contour(X,Y,func1,20);
  axis square;  axis off; 
  outline_domain(dom_type);
  title(title1); 
  
  subplot(222)
  trimesh(evt,xy(:,1),xy(:,2),expmesh); %mesh(X,Y,func1);          
  axis square;
  view(330,30);
  
  subplot(223)
  contour(X,Y,func2,20);    
  axis square;  axis off; 
  outline_domain(dom_type);
  title(title2);
   
  subplot(224)
  trimesh(evt,xy(:,1),xy(:,2),varmesh); % mesh(X,Y,func2);          
  axis square;
  view(330,30);
    
  set(findall(gcf,'-property','Fontsize'),'Fontsize',fontSize);
  
end % end child function

