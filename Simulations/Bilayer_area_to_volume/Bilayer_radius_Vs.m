function [x1, x2, r]=Bilayer_radius_Vs(a,b,c)

tic 

global mu sigmarel sigma runs 
mu = a; % mean of distribution
sigmarel = b; % CV of distribution
sigma = a*b; % stdev of distribution
runs = c; % number of runs with random valus

r = mu + sigma.*randn(runs,1); % random values within normal distribution

for i=1:runs;
    if r(i) <= 0;
        while r(i) <=0;
            r(i) = mu + sigma.*randn(1,1);
        end
    end
end

for i=1:runs;
    [xc, xd]=Bilayer_radius_Vs_function(r(i));
    x1(i,1)=xc;
    x2(i,1)=xd;
end
    
toc

end


function [x1, x2]=Bilayer_radius_Vs_function(a)

global A B C V l_s theta ls_s Vs
V = 1.2e-12; % m3 average volume of compartments, 1 nL
%l = 150.3e-6; % m, average length of compartments
l_s = 2*(4*V./(3*pi)).^(1./3); % m, average length of compartments assuming
%sphere
theta = 60*pi/180; % ?, interface angle

Vs = a;
ls_s = 2*(4*Vs./(3*pi)).^(1./3); % m, average length of sender compartment
A=ls_s./2;
B=l_s./2;
C=sin(2*theta);

x1 = sqrt((A^4*B^2*C^2 + A^2*B^4*C^2 - 2*sqrt(-A^6*B^6*C^4*(C^2-1)))/(A^4 + 4*A^2*B^2*C^2 - 2*A^2*B^2 + B^4));
x2 = sqrt((A^4*B^2*C^2 + A^2*B^4*C^2 + 2*sqrt(-A^6*B^6*C^4*(C^2-1)))/(A^4 + 4*A^2*B^2*C^2 - 2*A^2*B^2 + B^4));

A1 = pi*x1^2;
A2 = pi*x2^2;

end

