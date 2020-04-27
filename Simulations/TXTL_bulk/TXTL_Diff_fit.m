function Gmat = TXTL_Diff_fit(alpha,beta,K,V,Kd,n,Kdtot,alphamax,t2,tmat,I0,t0,t)

%% definitions of initial conditions and parameter
% intial concentration vector
c0=zeros(1,6);
c0(1)=I0;
k=[K,V,Kd,n,Kdtot,alphamax,t2,tmat];
%% solve ODE
[~,y]=ode23s(@(t,c)TXTL_Diff_fit_define(t,c,k), [0;(t+t0)], c0);

%% readout of simulated data
Gmat = alpha*y(2:end,6)+beta; %%out of the c vector, you choose which value you want to look at

end

function ydot = TXTL_Diff_fit_define(t, c, k)

%% variables
Is=c(1); Ib=c(2); Ir=c(3); %inducer
Pr=c(4); %active promoter
Gr=c(5); %unmature protein
Gmatr=c(6);

%% rate constants
K=k(1); V=k(2); Kd=k(3); n=k(4); Kdtot=k(5); alphamax=k(6); t2=k(7);tmat=k(8);

%% equations
Isdot=-K*(Is-Ib)./V;
Ibdot=-K*(2*Ib-Is-Ir)/V;
Irdot=-K*(Ir-Ib)/V-Ir./(Kd+Ir);
Prdot=Ir./(Kd+Ir);
Grdot=alphamax*(1-exp(-t/t2))*Pr.^n./(Kdtot.^n+Pr.^n)-Gr/tmat;
Gmatrdot=Gr/tmat;


%% calculation
ydot = [Isdot;Ibdot;Irdot;Prdot;Grdot;Gmatrdot];
end