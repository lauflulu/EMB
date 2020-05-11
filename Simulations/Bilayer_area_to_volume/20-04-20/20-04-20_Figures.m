%% 20-04-20 Figures

%% Figure 1: Bilayer area with different theta
mu=1.2e-12;
sigmarel=0.35;
Amu=8.6e-9;
Asigma=1.6e-9;
rA= Amu + Asigma.*randn(10000,1);
[x_60, A_60, r1_60, r2_60]=Bilayer_radius_2(mu, sigmarel, 10000); % theta = 60 ?C
[x_50, A_50, r1_50, r2_50]=Bilayer_radius_2(mu, sigmarel, 10000); % theta = 50 ?C
[x_40, A_40, r1_40, r2_40]=Bilayer_radius_2(mu, sigmarel, 10000); % theta = 40 ?C

[x, A, r1, r2]=Bilayer_radius_2(a,b,c)

mu = a; % mean of distribution
sigmarel = b; % CV of distribution
sigma = a*b; % stdev of distribution
runs = c; % number of runs with random valus

r1 = mu + sigma.*randn(runs,1); % random values within normal distribution
r2 = mu + sigma.*randn(runs,1); % random values within normal distribution


for i=1:runs;
    if r1(i) <= 0;
        while r1(i) <=0;
            r1(i) = mu + sigma.*randn(1,1);
        end
    end
    if r2(i) <= 0;
        while r2(i) <=0;
            r2(i) = mu + sigma.*randn(1,1);
        end
    end
end

for i=1:runs;
    [x1, A1]=Bilayer_radius_2_function(r1(i), r2(i));
    x(i,1)=x1;
    A(i,1)=A1;
end

[x1, A1]=Bilayer_radius_2_function(a, b)

global A B C D V l_s theta ls_s Vs
V = a; % m3 average volume of compartments, 1 nL
%l = 150.3e-6; % m, average length of compartments
l_s = 2*(4*V./(3*pi)).^(1./3); % m, average length of compartments assuming
%sphere
theta = 40*pi/180; % ?, interface angle % or 50, 60

Vs = b;
ls_s = 2*(4*Vs./(3*pi)).^(1./3); % m, average length of sender compartment

A=ls_s./2;
B=l_s./2;
C=cos(2*theta);

D = sqrt(A^2+B^2+2*A*B*C);

x1 = sqrt(-D^4+2*D^2*A^2+2*D^2*B^2-A^4+2*A^2*B^2-B^4)/(2*D);

A1 = pi*x1.^2;

%% Figure 2: same but with 40, 41, 42, 43
% theta is 42?