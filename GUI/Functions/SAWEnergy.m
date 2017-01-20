function [Energy] = SAWEnergy(X,Y,Z,X1,Y1,Z1,Const,a,b);

% Lennard-Jones potential evaluated over chromosome sets
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

% Const = [ve,vh,sigma,lb];
ve    = Const(1,1); % Electrostatic interaction [see Spakowitz & Wang (2005)]
vh    = Const(1,2); % Hardcore potential: ~10000.
sigma = Const(1,3); % Hydrodynamics radius of DNA (~30 nm for chromatin).
lD    = Const(1,4); % Debye length.

% Set 
D  = [];
le = length(X);
lo = length(X1);

mini = min(le,lo);
maxi = max(le,lo);

% Calculate pairwise distance between all beads.
for m = 1:mini(1,1) 
    for l = 1:maxi(1,1)
 Distances(l,m) = sqrt((X(1,l)-X1(1,m)).^2+(Y(1,l)-Y1(1,m)).^2+ ... 
                  (Z(1,l)-Z1(1,m)).^2);
    end
end

% Criterior for wether or not the interaction is upon a single chromosome
% or between two chromosomes.
val = abs(a-b)==0;

% Remove self-self pairiwse distances and ignore doubling up of pairwise
% distances.
for i = (1+val):lo
    D = [D;diag(Distances,i-1)];
end


Dhc = D(D<=sigma);
Dhc = Dhc(Dhc~=0);
Dsigma = D(D>sigma);

% Calculate potentials.
Vhc = (vh/12)*[(sigma./Dhc).^12 - 2*(sigma./Dhc).^6 + 1] + ...
    ve*exp(-Dhc/lD)./Dhc;
Vsigma = ve*exp(-Dsigma/lD)./Dsigma;

% Calculate energy.
Energy = 0.5*(sum(Vhc) + sum(Vsigma));
%toc