function [X,Y,Z_val] = SphereLines(Radius,x0,y0,z0,N);
%
% Generates grid to lie on the surface of a sphere of defined radius.
% Copyright (C) 2010  Christopher Penfold

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
Z = linspace(z0-Radius,z0+Radius,N);

for i = 1:N
    Z_val = ones(1,100).*Z(1,i);
    theta = linspace(0,1,100)*pi;
    %X = linspace(-Radius*sin(acos(Z(1,i)/Radius)),Radius*sin(acos(Z(1,i)/Radius)),100);
    X = linspace(-Radius*sin(acos(Z(1,i)/Radius)) - x0,x0,100);
    Y = real(sqrt(Radius^2 - Z_val.^2 - X.^2));
    
    hold on
    plot3(X,Y,Z_val,'-')
    plot3(X,Y,-Z_val,'-')
    plot3(X,-Y,-Z_val,'-')
    plot3(-X,-Y,Z_val,'-')
    plot3(-X,Y,-Z_val,'-')
    plot3(-X,-Y,-Z_val,'-')
  
   

    xlabel('x')
    ylabel('y')
    zlabel('z')
    
    
end
  X1 = zeros(1,100);
  Z1 = linspace(-Radius,Radius,100);
  Y1 = real(sqrt((Radius^2 - X1.^2 - Z1.^2)));
  
  X2 = real(sqrt(0.5*Radius^2 - 0.5*Z1.^2));
  Y2 = X2;
plot3(X1,Y1,Z1,'b-')
plot3(X1,-Y1,Z1,'b-')
plot3(Y1,X1,Z1,'b-')
plot3(-Y1,X1,Z1,'b-')
plot3(X2,-Y2,Z1,'b-')
plot3(-X2,-Y2,Z1,'b-')
plot3(-X2,Y2,Z1,'b-')
plot3(X2,Y2,-Z1,'b-')
plot3(-X1,Y1,-Z1,'b-')