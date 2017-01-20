function [xnew1,xnew11,phi1,lambda,theta1] = Link3D(N,a,X0,Xn);

%Calculate final link vector in chromosome path
%    Copyright (C) 2010  Christopher Penfold

%    This program is free software: you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation, either version 3 of the License, or
%    (at your option) any later version.
%
%    This program is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with this program.  If not, see <http://www.gnu.org/licenses/>.
x0 = X0(1,1);
y0 = X0(1,2);
z0 = X0(1,3);
xn = Xn(1,1);
yn = Xn(1,2);
zn = Xn(1,3);

dx = x0;
dy = y0;
dz = z0;

x0 = x0-dx;
y0 = y0-dy;
z0 = z0-dz;
xn = xn-dx;
yn = yn-dy;
zn = zn-dz;

r = sqrt((xn-x0)^2 + (yn-y0)^2 + (zn-z0)^2);

phi   = acos((zn-z0)/r);
X     = (xn-x0)/(r*sin(phi));
Y     = (yn-y0)/(r*sin(phi));
theta = atan2(Y,X);

phi1  = phi;
lambda=-theta;
theta1= theta;

RX = [cos(theta1),sin(theta1),0;-sin(theta1),cos(theta1),0;0,0,1];
RZ = [cos(lambda),sin(lambda),0;-sin(lambda),cos(lambda),0;0,0,1];
RY = [cos(phi1),0,-sin(phi1);0,1,0;sin(phi1),0,cos(phi1)];
Rot = RZ*RY*RX;
xnew = Rot*[X0;Xn]';
xnew = xnew';

a1 = r/2;
h  = sqrt(a^2 - a1^2);

th = rand(1,N)*2*pi;
zz = ones(1,N)*r/2;
yy = h*cos(th);
xx = h*sin(th);

xnew1 = inv(Rot)*[xx',yy',zz']';
xnew11 = inv(Rot)*xnew';

xnew1(3,:) = xnew1(3,:) + dz;
xnew1(2,:) = xnew1(2,:) + dy;
xnew1(1,:) = xnew1(1,:) + dx;