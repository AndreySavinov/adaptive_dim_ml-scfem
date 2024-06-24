function [psi,dpsidx,dpsidy] = vtqderiv(s,t,xl,yl)
%VTQDERIV  derivatives of super-quadratic shape functions for vector s,t
%   [psi,dpsidx,dpsidy] = vtqderiv(s,t,xl,yl);
%   input
%          s         reference element x coordinate   
%          t         reference element y coordinate
%          xl        physical element x vertex coordinates 
%          yl        physical element y vertex coordinates
%   output
%          psi       shape function
%          dpsidx    x derivative of psi
%          dpsidy    y derivative of psi
%
%    TIFISS function: DJS; 16 January 2016.
% Copyright (c) 2016 D.J. Silvester and Qifeng Liao
      nel=length(xl(:,1));
      zero_v = zeros(nel,1);
      one_v = ones(nel,1);
%
% evaluate shape functions 
     [xi,dxids,dxidt] = vtshape(s,t);
     [xqi,dxqids,dxqidt] = vtqshape(s,t);
%
      dxds = zero_v;
      dxdt = zero_v;
      dyds = zero_v;
      dydt = zero_v;
%	   jac = zero_v;
%    invjac = zero_v; 
%
      for ivtx = 1:3
         dxds(:) = dxds(:) + xl(:,ivtx) .* one_v.*dxids(:,ivtx);
         dxdt(:) = dxdt(:) + xl(:,ivtx) .* one_v.*dxidt(:,ivtx);
         dyds(:) = dyds(:) + yl(:,ivtx) .* one_v.*dxids(:,ivtx);
         dydt(:) = dydt(:) + yl(:,ivtx) .* one_v.*dxidt(:,ivtx);
      end
%
%      jac(:) = dxds(:).*dydt(:) - dxdt(:).*dyds(:);
%      invjac(:) = one_v ./ jac(:);
%
      for ivtx = 1:7
         psi(:,ivtx) = xqi(:,ivtx).*one_v;
		 dpsidx(:,ivtx) = ( dxqids(:,ivtx).*dydt(:) - dxqidt(:,ivtx).*dyds(:));
         dpsidy(:,ivtx) = (-dxqids(:,ivtx).*dxdt(:) + dxqidt(:,ivtx).*dxds(:));
      end
      return