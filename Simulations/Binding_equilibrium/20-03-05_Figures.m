%% Demonstrating that DNA/repressor protein/inducer binding is fast compared to TX-TL

% Figure 1: in the lowest concentration of LacI 1 nM (externally added) and
% DNA 1 nM, equilibration takes 50s. 
function [t,outcome]=Binding_equilibrium(a,b,c,d,e)

global kp km nspecies
nspecies=3;
kp=a;
km=b;


%initialization
for ispecies=1:1:nspecies;
    x0(ispecies)=0;
end

%initial values < > 0
x0(1)=c; %M, reservoir
x0(2)=d;
x0(3)=e;

% duration of experiment
dt=1;
dx(1)=-kp*y(1)*y(2) + km*y(3);
dx(2)=-kp*y(1)*y(2) + km*y(3);
dx(3)=kp*y(1)*y(2) - km*y(3);

a=5.1e6;b=3.7e-1;c=1e-9;d=1e-9;e=0;
[t,outcome]=Binding_equilibrium(a,b,c,d,e);
plot(t, outcome);

%Figure 2: same but with high concentrations 5 nM DNA and 6 uM LacI

%Figure 3: same but for LacI binding IPTG (1 nM LacI and 1 nM
%IPTG as little amounts) with 
a=1.2e5;
b=2.1e-1;
c=100e-9;
d=1e-9;
e=0; and same script
