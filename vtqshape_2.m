function [dpsi2ds2,dpsi2dt2,dpsi2dsdt] = vtqshape_2(s,t)
%VTQSHAPE_2  second derivatives on reference element for vector s, t
%   [dpsi2ds2,dpsi2dt2,dpsi2dsdt] = vtqshape_2(s,t);
%   input
%          s         first triangle coordinate   
%          t         second triangle coordinate  
%  
%   output
%          dpsi2ds2     second derivative in s
%          dpsi2dt2     cross derivative
%          dpsi2dsdt    second derivative in t
%
%   TIFISS function: DJS; 16 January 2016.
% Copyright (c) 2016 D.J. Silvester and Qifeng Liao
      n=length(s);
      one = 1.0e0*ones(n,1); zero=0.0e0*one;
%      xi(:,1) = one-s-t;
%      xi(:,2) = s;
%      xi(:,3) = t;
%      dxids(:,1) = -one;
%      dxids(:,2) = one;
%      dxids(:,3) = zero;
%      dxidt(:,1) = -one;
%      dxidt(:,2) = zero;
%      dxidt(:,3) = one;
%
% 	  psi(:,1) = (2*xi(:,1)-1).*xi(:,1);
%	  psi(:,2) = (2*xi(:,2)-1).*xi(:,2);
%	  psi(:,3) = (2*xi(:,3)-1).*xi(:,3);
%	  psi(:,4) = 4*xi(:,2).*xi(:,3);
%	  psi(:,5) = 4*xi(:,1).*xi(:,3);
%	  psi(:,6) = 4*xi(:,1).*xi(:,2);
%     psi(:,7) = 27*xi(:,1).*xi(:,2).*xi(:,3);

      dpsi2ds2(:,1)=4*one;
      dpsi2ds2(:,2)=4*one;
      dpsi2ds2(:,3)=zero;
      dpsi2ds2(:,4)=zero;
      dpsi2ds2(:,5)=zero;
      dpsi2ds2(:,6)=-8*one;
      
      dpsi2dt2(:,1)=4*one;
      dpsi2dt2(:,2)=zero;
      dpsi2dt2(:,3)=4*one;
      dpsi2dt2(:,4)=zero;
      dpsi2dt2(:,5)=-8*one;
      dpsi2dt2(:,6)=zero;
      
      dpsi2dsdt(:,1)=4*one;
      dpsi2dsdt(:,2)=zero;
      dpsi2dsdt(:,3)=zero;
      dpsi2dsdt(:,4)=4*one;
      dpsi2dsdt(:,5)=-4*one;
      dpsi2dsdt(:,6)=-4*one;
  return
     