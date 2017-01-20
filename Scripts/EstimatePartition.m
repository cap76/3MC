function [Z] = EstimatePartition(Parameters);

No     = Parameters(1,1);
len    = Parameters(1,2);
Radius = Parameters(1,3);
U = Parameters(1,4);
a = Parameters(1,5);
k = Parameters(1,6:end);

J = k/a;

if max(J)>600 %Function besseli is unstable for J>700 or thereabouts, 
    %therefore use symbolic toolbox after this point, which is slower.

n       = sym(0.5);
x       = sym(J);
logZest = double(log((1./sqrt(2*pi*J))) + (log(vpa(besseli(n,x)))));
Z       = sum(logZest);

else

Zest = (1./sqrt(2*pi*J)).* besseli(0.5,J);
Z    = sum(log(Zest));

end

return