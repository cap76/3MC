function [Pos] = BBRandFlight(No,D,X0,X)

%Function for generating a `Brownian Bridge' led random flight.
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
% Define origin.
xor = X0(1,1);
yor = X0(1,2);
zor = X0(1,3);

% Define end.
xn = X(1,1)-xor;
yn = X(1,2)-yor;
zn = X(1,3)-zor;

% Physical constants.
Dis = D;
N   = No;
NT  = No;
%n = N - 2;

% Vectorise steps.
Pos      = zeros(N,3);
Pos(N,:) = [xn,yn,zn];
Pos(1,:) = [0,0,0];
condition = (4*Dis);
valind = 0;
while condition > (2*Dis)
valind = valind +1;
    for i = 2:N-2
O = Pos(i-1,:);
E = Pos(N,:);
NT = N-i-1-valind;%35;

   MEA = O + ((1-0)/(max(NT,1) - 0))*(E - O);
   VAR = (((1 - 0)*(max(NT,1) - 1))/(max(NT,1)-0))*diag(ones(1,3));
   X_n = gsamp(MEA,VAR,1);
   Pos(i,:) = X_n;

   PosRF = sqrt(sum((Pos(i,:)-Pos(i-1,:)).^2));
   RF(i,:) = Dis*((Pos(i,:)-Pos(i-1,:))/PosRF);
   
   RFDist(1:i,:) = cumsum(RF(1:i,:),1);
   Pos(i,:) = RFDist(i,:);
end
   condition = sqrt(sum(([xn,yn,zn] - Pos(N-2,:)).^2));
%disp([num2str(condition -(2*Dis))])
%keyboard
end

[xnew1,xnew11,phi1,lambda,theta1] = Link3D(1,Dis,Pos(N-2,:),[xn,yn,zn]);
Pos(N-1,:) = xnew1;
Pos(:,1) = Pos(:,1) + xor;
Pos(:,2) = Pos(:,2) + yor;
Pos(:,3) = Pos(:,3) + zor;
%plot3([Pos(:,1)],[Pos(:,2)],[Pos(:,3)],'.-')