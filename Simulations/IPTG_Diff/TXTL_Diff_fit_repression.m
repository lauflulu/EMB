function Gmat = TXTL_Diff_fit_repression(alpha,beta,K,V,kon,koff,n,Kdtot,alphamax,t2,tmat,I0,R0,t0,t)

%% definitions of initial conditions and parameter

% intial concentration vector
c0=zeros(1,7);
c0(1)=I0;
c0(4)=R0;
k=[K,V,kon,koff,n,Kdtot,alphamax,t2,tmat];
%% solve ODE
[~,y]=ode23s(@(t,c)TXTL_Diff_fit_define(t,c,k), [0;(t+t0)], c0);

%% readout of simulated data
Gmat = alpha*y(2:end,7)+beta; %%out of the c vector, you choose which value you want to look at

end

function ydot = TXTL_Diff_fit_define(t, c, k)

%% variables
Is=c(1); Ib=c(2); Ir=c(3); %inducer
Rr=c(4);
IRr=c(5); %active promoter
Gr=c(6); %unmature protein
Gmatr=c(7);

%% rate constants
K=k(1); V=k(2); kon=k(3); koff=k(4); n=k(5); Kdtot=k(6); alphamax=k(7); t2=k(8); tmat=k(9);

%% equations
Isdot=-K*(Is-Ib)./V;
Ibdot=-K*(2*Ib-Is-Ir)/V;
Irdot=-K*(Ir-Ib)/V-kon*Ir*Rr+koff*IRr;
Rrdot=-kon*Ir*Rr+koff*IRr;
IRrdot=kon*Ir*Rr-koff*IRr;
Grdot=alphamax*(1-exp(-t/t2))*Kdtot.^n./(Kdtot.^n+Rr.^n)-Gr/tmat;
Gmatrdot=Gr/tmat;


%% calculation
ydot = [Isdot;Ibdot;Irdot;Rrdot;IRrdot;Grdot;Gmatrdot];
end