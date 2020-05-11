function output=EMB_Cont_Diff(a);

global n0 x0 L D x t k
n0 = 1; % initial concentration
x0 = 1; % initial distance in which signal is contained
L = 10; % total length of the array
D = 1; % diffusion (is this the efficient diffusion?)
x = L/2; % position of interest
t = 1; % time point
a = 10;

y = (n0/2)*symsum(erf((x0+2*k*L-x)/sqrt(D*t)) + erf((x0-2*k*L+x)/sqrt(D*t)),k,-a,a)
output=double(y);
end

