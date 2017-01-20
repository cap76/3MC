function [X,Y,Z] = RandFlightCon(X0,XN,a,len);

%Generates a Random Flight with fixed ends using Monte Carlo rejection.
%Note that this is relativelty fast when |X0-Xn|<<a*len but slow otherwise.
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

%Start positions.
x0 = X0(1,1);
y0 = X0(1,2);
z0 = X0(1,3);

%End positions.
xn = XN(1,1);
yn = XN(1,2);
zn = XN(1,3);

Dist = Inf;

%Empty array for storing configurations.
X = zeros(1,len);
Y = zeros(1,len);
Z = zeros(1,len);

%Begin loop to generate random samples until one is found with the end and
%second from the end very close to one another.
while Dist>(2*a)^2
[rx1,ry1,rz1,theta,phi] = fSpherePoint(len-3,a);
xx = [x0,cumsum(rx1)+x0];
yy = [y0,cumsum(ry1)+y0];
zz = [z0,cumsum(rz1)+z0];
Dist = (xx(1,len-2)-xn)^2 + (yy(1,len-2)-yn)^2 + (zz(1,len-2)-zn)^2;
end

X(1,1:len-2) = xx;
Y(1,1:len-2) = yy;
Z(1,1:len-2) = zz;

%Calculate the missing vector.
[xnew1,xnew11] = Link3D(1,a,[X(1,len-2),Y(1,len-2),Z(1,len-2)],[xn,yn,zn]);

X(1,len-1) = xnew1(1,1);
Y(1,len-1) = xnew1(2,1);
Z(1,len-1) = xnew1(3,1);

X(1,len) = xn;
Y(1,len) = yn;
Z(1,len) = zn;