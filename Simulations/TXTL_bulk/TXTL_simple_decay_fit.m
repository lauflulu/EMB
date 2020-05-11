function Gmat=TXTL_simple_decay_fit(alpha,beta,alphamax,tRNA,tmat, t1, t2,t0,t)

%% definitions of initial conditions and parameter

% intial concentration vector
c0=zeros(1,2);
k=[alphamax,tRNA,tmat, t1, t2];
%% solve ODE
[~,y]=ode23s(@(t,c)TXTL_simple_decay_fit_define(t,c,k), [0;(t+t0)], c0);

%% readout of simulated data
Gmat = alpha*y(2:end, 2)+beta; %%out of the c vector, you choose which value you want to look at

end

function ydot=TXTL_simple_decay_fit_define(t,c,k);
%% variables
G=c(1); %unmature protein
Gmat=c(2);

%% rate constants
alphamax=k(1); tRNA=k(2); tmat=k(3); t1=k(4); t2=k(5); 

%% equations

if t < t1;
    Gdot=alphamax*(1-exp(-t/tRNA))-G./tmat;
else;
    Gdot=alphamax*(1-exp(-t/tRNA))*exp(-(t-t1)/t2)-G./tmat;
end
Gmatdot=G./tmat;


%% calculation
ydot = [Gdot;Gmatdot];
end