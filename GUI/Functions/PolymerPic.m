function PolymerPic(Parameters,X1,Y1,Z1,Width)

%Graphical plot of an example chromosome path
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

Radius = Parameters(1,4);
c  = 2;
l = size(X1,1);
ll = randint(1,1,l)+1;

x0      = [X1(ll,:)]/Radius;
y0      = [Y1(ll,:)]/Radius;
z0      = [Z1(ll,:)]/Radius;

Width = Width/(Radius);

x1 = x0; y1 = y0; z1 = -c*ones(1,100);
x2 = x0; y2 = c*ones(1,100); z2 = z0;
x3 = c*ones(1,100); y3 = y0; z3 = z0;

x0      = [x0];
y0      = [y0];
z0      = [z0];

%Take a guess at the width of the chromosome. Unless specified take this to
%be the distance of the link vectors.
delta = Width;%sqrt(abs(x0(1,1)-x0(1,2)).^2 + abs(y0(1,1)-y0(1,2)).^2 + abs(z0(1,1)-z0(1,2)).^2);

lighting gouraud
colormap autumn
[x y z] = sphere(40);
for i = 1:length(x0)
s = surf(x*delta+x0(1,i),y*delta+y0(1,i),z*delta+z0(1,i),'LineStyle','None');
hold on
daspect([1 1 1])
end

camlight(40,38)
%return
%hold on
%[X_val,Y_val,Z_val] = SphereLines(1+Width,0,0,0,55);
%hold on
%plot3(x0,y0,-1.5*ones(1,100),'.-')
%plot3(x0,1.5*ones(1,100),z0,'.-')
%plot3(1.5*ones(1,100),y0,z0,'.-')


%Plot encasing sphere for comparison.
%a = linspace(-1,1,1000);
%b = sqrt(1.^2 - a.^2);
%aa = [a,-a];
%bb = [b,-b];

%size(aa)
%size(bb)

%hold on
%plot3(aa,bb,-c*ones(1,2000),'-');
%plot3(aa,c*ones(1,2000),bb,'-');
%plot3(c*ones(1,2000),aa,bb,'-');
%axis on
%zlim([-c c])
%ylim([-c c])
%xlim([-c c])