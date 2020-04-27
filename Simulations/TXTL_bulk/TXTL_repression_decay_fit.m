function Gmat=TXTL_repression_decay_fit(alpha,beta,Kd,n,alphamax,tRNA,tmat, t1, t2,R0,t0,t)

%% definitions of initial conditions and parameter

% intial concentration vector
c0=zeros(1,2);
k=[Kd,n,alphamax,tRNA,tmat, t1, t2, R0];
%% solve ODE
[~,y]=ode23s(@(t,c)TXTL_repression_decay_fit_define(t,c,k), [0;(t+t0)], c0);

%% readout of simulated data
Gmat = alpha*y(2:end, 2)+beta; %%out of the c vector, you choose which value you want to look at

end

function ydot=TXTL_repression_decay_fit_define(t,c,k);
%% variables
G=c(1); %unmature protein
Gmat=c(2);

%% rate constants
Kd= k(1); n=k(2); alphamax=k(3); tRNA=k(4); tmat=k(5); t1=k(6); t2=k(7); R0=k(8);

%% equations

if t < t1;
    Gdot = alphamax*(1-exp(-t/tRNA))*(Kd.^n./(Kd.^n+R0.^n))-G./tmat;
else;
    Gdot=alphamax*(1-exp(-t/tRNA))*(Kd.^n./(Kd.^n+R0.^n))*exp(-(t-t1)/t2)-G./tmat;
end
Gmatdot=G./tmat;


%% calculation
ydot = [Gdot;Gmatdot];
end