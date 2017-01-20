function [Energy,Dis,Distances] = ExcludedVolumeEnergy(X1,Y1,Z1,X2,Y2,Z2,Radius,a,b);

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
len1 = length(X1);
len2 = length(X2);

Distances = sqrt((repmat(X1',1,len2) - repmat(X2,len1,1)).^2 + ...
            (repmat(Y1',1,len2) - repmat(Y2,len1,1)).^2 + ...
            (repmat(Z1',1,len2) - repmat(Z2,len1,1)).^2);

Dis = zeros(len1,len2);

Dis(find(Distances>Radius)) = 0;
Dis(find(Distances<=Radius)) = 10000;

if abs(a-b)==0
    D = [];
for i = 1:len1+1
    D = [D;diag(Dis,i)];
end
Energy = sum(D);
else
    Energy = sum(sum(Dis));
end

Energy;
