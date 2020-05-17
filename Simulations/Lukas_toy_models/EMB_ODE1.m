function ydot = EMB_ODE1(t, y, I, X, A, V, P, alphamax, alpha, kHill, k, tmat, deg)
%% re-k
aLacI=alpha(1); aRFP=alpha(2); aTetR=alpha(3); 
aYFP=alpha(4); tau=alpha(5); tRNA=alpha(6); a0=alpha(7);
KdTetR=kHill(1); nTetR=kHill(2); KdLacI=kHill(3); nLacI=kHill(4);
kon=k(1); koff=k(2);
tmat1=tmat(1); tmat2=tmat(2);
kdeg=deg(1); Kdeg=deg(2);

%% reshape y
g=reshape(y,I,X);
dg=zeros(I,X);

%% IPTG diffusion bilayer limited (P*l/D=2e-4<<1)
dg(1,1)=-P/V(1,1)*A(1,1)*(g(1,1)-g(1,2));
for x=2:X-1
    dg(1,x)=-P/V(1,x)*(A(1,x-1)*(g(1,x)-g(1,x-1))+A(1,x)*(g(1,x)-g(1,x+1)));  
end
dg(1,X)=-P/V(1,X)*A(1,X-1)*(g(1,X)-g(1,X-1));

%% reactions 
for x=2:X
    %% gene expression
    dg(2,x)=dg(2,x)+alphamax(1,x)*aLacI*(1-exp(-t/tRNA))*(a0+1./(1+(g(6,x)./KdTetR).^nTetR))*EMB_decay_v1(t,tau);% lacI
    dg(4,x)=dg(4,x)+alphamax(1,x)*aRFP*(1-exp(-t/tRNA))*(a0+1./(1+(g(6,x)./KdTetR).^nTetR))*EMB_decay_v1(t,tau);% RFP
    dg(6,x)=dg(6,x)+alphamax(1,x)*aTetR*(1-exp(-t/tRNA))*(a0+1./(1+(g(2,x)./KdLacI).^nLacI))*EMB_decay_v1(t,tau);% tetR
    dg(7,x)=dg(7,x)+alphamax(1,x)*aYFP*(1-exp(-t/tRNA))*(a0+1./(1+(g(2,x)./KdLacI).^nLacI))*EMB_decay_v1(t,tau);% YFP
    
    %% maturation
    dg(4,x)=dg(4,x)-g(4,x)/tmat1;
    dg(5,x)=dg(5,x)+g(4,x)/tmat1; % RFP mat
    dg(7,x)=dg(7,x)-g(7,x)/tmat2;
    dg(8,x)=dg(8,x)+g(7,x)/tmat2; % YFP mat
    
    %% other reactions
    dg(1,x)=dg(1,x)-kon*g(1,x)*g(2,x)+koff*g(3,x);
    dg(2,x)=dg(2,x)-kon*g(1,x)*g(2,x)+koff*g(3,x);
    dg(3,x)=dg(3,x)+kon*g(1,x)*g(2,x)-koff*g(3,x);
    
    %% degradation
    dg(1,x)=dg(1,x)-EMB_degrade_v1(kdeg,Kdeg,g(3,x))*EMB_decay_v1(t,tau); % degradation of lacI*IPTG frees IPTG
    dg(2,x)=dg(2,x)+EMB_degrade_v1(kdeg,Kdeg,g(2,x))*EMB_decay_v1(t,tau);
    dg(3,x)=dg(3,x)+EMB_degrade_v1(kdeg,Kdeg,g(3,x))*EMB_decay_v1(t,tau);
    dg(4,x)=dg(4,x)+EMB_degrade_v1(kdeg,Kdeg,g(4,x))*EMB_decay_v1(t,tau);
    dg(5,x)=dg(5,x)+EMB_degrade_v1(kdeg,Kdeg,g(5,x))*EMB_decay_v1(t,tau);
end

%% reshape back
icount=0;
for x=1:1:X
    for ispecies=1:1:I
        icount=icount+1;
        dy(icount)=dg(ispecies,x);
    end
end
ydot=dy';

end