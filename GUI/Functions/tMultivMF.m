function[x,y,z,theta,phi] = tMultivMF(ui,uj,uk,kappa,Radius,loops);

% Application of rejection sampling to generate samples from a von
% Mises-Fisher distribition embedded in R3. Returns user defined number of
% samples.
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

% Initialise variables. Some dummy variables for initial loop.
count       = 0;
total       = 0;
Xold        = []; 
Yold        = []; 
Zold        = []; 
THETAold    = []; 
PHIold      = [];

% Begin loop.
while total < loops
% Generate random positions on a sphere under uniform distribution.
[X_1,Y_1,Z_1,THETA_1,PHI_1]     = fSpherePoint(8000,Radius);
%X   = [(sin(THETA_1).*sin(PHI_1))',(cos(THETA_1).*sin(PHI_1))', ...
%    (cos(THETA_1))']';
X   = [(sin(PHI_1).*sin(THETA_1))',(sin(PHI_1).*cos(THETA_1))', ...
    (cos(PHI_1))']';
mu  = [ui,uj,uk]';
%[X_1,Y_1,Z_1,THETA_1,PHI_1]     = Points_on_Sphere(Radius,8000);
%X   = [(sin(THETA_1).*sin(PHI_1))',(cos(THETA_1).*sin(PHI_1))', ...
%    (cos(PHI_1))']';
%mu  = [ui,uj,uk]';


% Generate the maximal probability under the vMF distribution.
max_vMF         = (kappa^0.5 * exp(kappa * mu'*mu))/((2*pi)^(3/2) ...
    * besseli(0.5,kappa));
vMisesFisher    = (kappa^0.5 * exp (kappa * mu'*X))/((2*pi)^(3/2) ...
    * besseli(0.5,kappa));

% Rejection Sampling Criterion.
Condition   = rand(1,8000).*max_vMF;
X_11        = X_1(find(Condition < vMisesFisher));
total       = total + length(X_11);
Xold        = [Xold,X_1(find(Condition < vMisesFisher))];
Yold        = [Yold,Y_1(find(Condition < vMisesFisher))];
Zold        = [Zold,Z_1(find(Condition < vMisesFisher))];
THETAold    = [THETAold,THETA_1(find(Condition < vMisesFisher))];
PHIold      = [PHIold,PHI_1(find(Condition < vMisesFisher))];
end

x       = Xold(1,1:loops);
y       = Yold(1:loops);
z       = Zold(1:loops);
theta   = THETAold(1:loops);
phi     = PHIold(1:loops);