function [pm]=patchfault(m,i,j);
%PATCHFAULT    [pm]=patchfault(m,i,j)
%
%This function discretizes a fault model into i*j distinct patches.
%
%INPUTS:
%
%    m = 1x10 vector defined as follows
%
%       m(1) = fault length along the strike direction (km)
%       m(2) = fault width in dip direction (km)
%       m(3) = if sin(d) > 0, m(3) = depth to lower edge of fault (km)
%              if sin(d) < 0, m(3) = depth to upper edge of fault (km)
%              where d = dip angle
%       m(4) = dip angle, from the horizontal (degrees) 
%       m(5) = strike, clockwise from N (degrees)
%       m(6) = East offset of midpoint of lower edge from origin (km)
%       m(7) = North offset of midpoint of lower edge from origin (km)
%       m(8) = strike slip (negative for right lateral, 0 if N/A)
%       m(9) = dip slip:    if sin(2d) < 0, then m(9) > 0 --> normal slip
%                                       and m(9) < 0 --> reverse slip
%                           if sin(2d) > 0, then m(9) > 0 --> reverse
%                                       and m(9) < 0 --> normal
%               where d = dip angle
%               m(9) 0 if dip-slip N/A
%       m(10) = tensile motion (positive for opening, 0 if N/A)
%
%   Note: See Okada (1985) for details on geometry parameterization.
%
%    i = number of patches along fault length
%
%    j = number of patches along fault width
%
%OUTPUTS:
%
%    pm = Nsx10 matrix of patch models, where Ns=i*j
%
%    NOTE: For backwards compatibility, if m has length 7 instead of
%         10 then the pm returned will be Ns x 7 array.
%
%Examples:
% ---> right lateral strike-slip <---
% m = [100, 10, 10, 90, 0, 0, 0, -1, 0, 0]
% i = 2
% j = 2
% 
% pm(1,:) = [50, 5, 5, 90, 0, 0, -25, -1, 0, 0]
% pm(2,:) = [50, 5, 10, 90, 0, 0,-25, -1, 0, 0]
% pm(3,:) = [50, 5, 5, 90, 0, 0, 25, -1, 0, 0]
% pm(4,:) = [50, 5, 10, 90, 0, 0, 25, -1, 0, 0]
%
% ---> normal slip <---
% m = [100, 10, 10, 45, 75, 10, 10, 0, -1, 0]
% i = 2
% j = 2
% 
% pm(1,:) = [50, 5, 6.4645, 45, 75, -15.0632, 6.9446, 0, -1, 0]
% pm(2,:) = [50, 5, 10, 45, 75, -14.1481, 3.5295, 0, -1, 0]
% pm(3,:) = [50, 5, 6.4645, 45, 75, 33.2331, 19.8855, 0, -1, 0]
% pm(4,:) = [50, 5, 10, 45, 75, 34.1481, 16.4705, 0, -1, 0]

%---------------------------------------------------------------------
% Date modified         by          comments
%---------------------------------------------------------------------
%   ???                 ???         original code
%   Feb. 18, 2003       JRM         Included slip components in input
%                                   model vector and output patch 
%                                   models; expanded help message to
%                                   explain dip, depth, and dip slip;
%                                   also corrected order of input for
%                                   repmat
%   Apr. 21, 2003	JRM         Added backwards compatibility for
%				    input which only gives m vector
%				    of length 7.
%---------------------------------------------------------------------
 
%Set constants

	dip=m(4)*pi/180;
	strike=-m(5)*pi/180;
	sin_dip=sin(dip);
	cos_dip=cos(dip);
	iw=m(1)/i;
	jw=m(2)/j;
	is=(1:i);
	js=(1:j)';
	n=i*j;
	c1=-m(2)*cos_dip;
	c2=0.5*(m(1)+iw);
	c3=m(3)-j*jw*sin_dip;
	
%Calculate midpoints, depths of each patch
	
	p=[cos_dip*(jw*js-m(2))]*ones(1,i);
        q=ones(j,1)*[(iw*is)-0.5*(m(1)+iw)];
	r=[m(3)-jw*sin_dip*(j-js)]*ones(1,i);
	mp=[p(:),q(:),r(:)];

%Adjust midpoints for strike
	
	R=[cos(strike),-sin(strike),0;sin(strike),cos(strike),0;0,0,1];
	mp=mp*R';

%Adjust midpoints for offset from origin

	mp(:,1)=mp(:,1)+m(6);
	mp(:,2)=mp(:,2)+m(7);

%Form patch-models

	pm(:,1)=ones(n,1)*iw;
	pm(:,2)=ones(n,1)*jw;
	pm(:,3)=mp(:,3);
	pm(:,4:5)=ones(n,1)*m(4:5);
	pm(:,6:7)=mp(:,1:2);
    
%Append slip component information

    if length(m) == 10
        pm(:,8:10)=repmat(m(8:10),n,1);
    end


