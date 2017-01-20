function [x,y,z,theta,phi] = Points_on_Sphere(radius,number)

% Generates positions on the surface of a sphere of radius R, centered at
% [X,Y,Z] = [0,0,0]. Returns [X,Y,Z,Theta,Phi] coordinates.
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
z = rand(1,number)*2*radius - radius;
phi = rand(1,number)*2*pi;
theta =  asin(z/radius);
x = radius*cos(theta).*cos(phi);
y = radius*cos(theta).*sin(phi);