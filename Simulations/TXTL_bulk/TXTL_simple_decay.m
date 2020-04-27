function [t,outcome]= TXTL_simple_decay(a,b)

global  nspecies alphamax tRNA tmat t1 t2
nspecies = 2; % P, Pmat
alphamax = 0.0036; % expression factor, determined from bulk and droplet fitting
tRNA = 15; % min, lifetime of RNA, from Karzbrun et al.
tmat = 7; % min, maturation time of GFP, from Iizuka et al.
t1 = a; 
t2 = b; 


%initialization
for ispecies=1:1:nspecies;
    x0(ispecies, 1)=0;
end


% duration of experiment
dt=10;
duration=480; % min, duration of experiment: 8h

% simulation run
tspan=[0:dt:duration];

options=odeset('AbsTol',1e-9,'RelTol',1e-5); % set tolerances
[t,y]=ode23s(@TXTL_decay_function,tspan,x0,options);

%output
t=t;
outcome=y(:,:);

end


function ydot = TXTL_decay_function(t, y, flags)

global  nspecies alphamax tRNA tmat t1 t2

if t < t1;
    dy(1,1)=alphamax*(1-exp(-t/tRNA))-y(1,1)/tmat;
else;
    dy(1,1)=alphamax*(1-exp(-t/tRNA))*exp(-(t-t1)/t2)-y(1,1)/tmat;
end
dy(2,1)=y(1,1)/tmat;

ydot=dy;

end