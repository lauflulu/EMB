%% EMB Geometrical parameters

%% parameters from Aurore

V1 = 1e-12; % m3
A1 = 2e-8; % m2
%% parameters from Lukas
V2 = 1.2e-12; % m3 average volume of compartments, 1 nL
A2 = 8.6e-9; % m2 bilayer area, identical to area of unspecific diffusion

%% droplet diameter from volume assuming spherical droplet
l1_sphere = 2*(4*V1./(3*pi)).^(1./3); % m, average length of compartments
l2_sphere = 2*(4*V2./(3*pi)).^(1./3); % m, average length of compartments

%% droplet diameter from volume assuming truncated sphere droplet
% the non-truncated droplet with volume Vtot has radius (l/2)+r_bil
% the truncation corresponds to 2 half-spheres with radius r_bil
% the total volume is V2 + (3*pi/4)*(r_bil)^3
% r_bil is sqrt(A2/pi)
% from this l is: 
l1_tsphere = 2.*((((4*V1)./(3*pi))+(sqrt(A1./pi)).^3).^(1/3)-sqrt(A1./pi));
l2_tsphere = 2.*((((4*V2)./(3*pi))+(sqrt(A2./pi)).^3).^(1/3)-sqrt(A2./pi));
