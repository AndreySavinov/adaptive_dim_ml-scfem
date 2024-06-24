function [xi,dxids,dxidt] = vtshape(s,t)
%VTSHAPE  linear shape functions for vector s, t
%   [xi,dxids,dxidt] = vtshape(s,t);
%   input
%          s         first triangle coordinate   
%          t         second triangle coordinate  
%  
%   output
%          xi        shape function
%          dxids     s derivative of xi
%          dxidt     t derivative of xi
%
%    TIFISS function: QL; 17 April 2011.
% Copyright (c) 2011 D.J. Silvester and Qifeng Liao
      n=length(s);
      one = 1.0e0*ones(n,1); zero=0.0e0*one;
      xi(:,1) = one-s-t;
      xi(:,2) = s;
	  xi(:,3) = t;
      dxids(:,1) = -one;
      dxids(:,2) = one;
      dxids(:,3) = zero;
      dxidt(:,1) = -one;
      dxidt(:,2) = zero;
      dxidt(:,3) = one;
      return