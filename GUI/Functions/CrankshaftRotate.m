function [X,Y,Z] = CrankshaftRotate(X0,Y0,Z0);

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

lenB            = randint2(1,1,len-2)+2;
maxi            = len - lenB;
Bond            = randint2(1,1,maxi)+1;
Bond2           = Bond + lenB;

dx      = X(1,Bond);
dy      = Y(1,Bond);
dz      = Z(1,Bond);
Xsub    = [X(1,Bond+1:Bond2)'-dx,Y(1,Bond+1:Bond2)'-dy,Z(1,Bond+1:Bond2)'-dz];
usub    = [X(1,Bond2) - X(1,Bond),Y(1,Bond2) - Y(1,Bond),Z(1,Bond2) - Z(1,Bond)];
udiv    = sqrt((X(1,Bond2) - X(1,Bond)).^2 + (Y(1,Bond2) - Y(1,Bond)).^2 + (Z(1,Bond2) - Z(1,Bond)).^2);
[Pr]    = Quaternion3(rand(1,1)*2*pi,usub/udiv,Xsub);


X(1,Bond+1:Bond2)   = Pr(:,1)' + dx;
Y(1,Bond+1:Bond2)   = Pr(:,2)' + dy;
Z(1,Bond+1:Bond2)   = Pr(:,3)' + dz;