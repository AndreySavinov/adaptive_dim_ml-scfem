function [dpsi2dx2,dpsi2dy2] = vtqderiv_2(s,t,xl,yl)
%VTQDERIV_2 second derivatives of quadratic shape functions for vector s, t
%   [dpsi2dx2,dpsi2dy2] = vtqderiv_2(s,t,xl,yl);
%   input
%          s         reference element x coordinate   
%          t         reference element y coordinate
%          xl        physical element x vertex coordinates 
%          yl        physical element y vertex coordinates  
%   output
%          dpsi2dx2    x second derivatives of psi
%          dpsi2dy2    y second derivatives of psi
%
%   TIFISS function: QL; 17 April 2011.
% Copyright (c) 2011 D.J. Silvester and Qifeng Liao
      nel=length(xl(:,1));
      zero_v = zeros(nel,1);
      one_v = ones(nel,1);
%
% evaluate bilinear shape functions 
     [phi_e,dphids,dphidt] = vtshape(s,t);
% evaluate second direvatives for biquadratic shape funtions
     [dpsi2ds2,dpsi2dt2,dpsi2dsdt] = vtqshape_2(s,t);
% local derivatives
      dxds = zero_v;
      dxdt = zero_v;
      dyds = zero_v;
      dydt = zero_v;
%
      for ivtx = 1:3
         dxds(:) = dxds(:) + xl(:,ivtx) .* one_v.*dphids(:,ivtx);
         dxdt(:) = dxdt(:) + xl(:,ivtx) .* one_v.*dphidt(:,ivtx);
         dyds(:) = dyds(:) + yl(:,ivtx) .* one_v.*dphids(:,ivtx);
         dydt(:) = dydt(:) + yl(:,ivtx) .* one_v.*dphidt(:,ivtx);
      end
%
     for ivtx = 1:6
		 dpsi2dx2(:,ivtx) = ( dpsi2ds2(:,ivtx).*dydt(:).^2 - 2*dpsi2dsdt(:,ivtx).*dydt(:).*dyds(:) + dpsi2dt2(:,ivtx).*dyds(:).^2);
         dpsi2dy2(:,ivtx) = ( dpsi2ds2(:,ivtx).*dxdt(:).^2 - 2*dpsi2dsdt(:,ivtx).*dxdt(:).*dxds(:) + dpsi2dt2(:,ivtx).*dxds(:).^2);
     end
        
   return   