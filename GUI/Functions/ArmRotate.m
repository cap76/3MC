function [X,Y,Z] = ArmRotate(X0,Y0,Z0);

%[X,Y,Z] = ARMROTATE(X0,Y0,Z0) rotates a chromosome arm about a
%chosen bead whilst confining the telomere to the surface of the nuclear
%periphery.
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
X               = X0;
Y               = Y0;
Z               = Z0;
len             = length(X0);
bead            = randint2(1,1,len-1) + 1;
dx              = X0(1,bead);
dy              = Y0(1,bead);
dz              = Z0(1,bead);
Xsub            = [X0(1,bead:len)'-dx, Y0(1,bead:len)'-dy, Z0(1,bead:len)'-dz];
usub            = [X0(1,bead),Y0(1,bead),Z0(1,bead)]/sqrt(X0(1,bead).^2+Y0(1,bead).^2+Z0(1,bead).^2);
[Pr]            = Quaternion3(rand(1,1)*2*pi,usub,Xsub);
X(1,bead:len)   = Pr(:,1)' + dx;
Y(1,bead:len)   = Pr(:,2)' + dy;
Z(1,bead:len)   = Pr(:,3)' + dz;