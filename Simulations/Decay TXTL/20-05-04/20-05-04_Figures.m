%% 20-05-04 Figures

%% Figure 1: Trying all 3 topologies with a kdeg = 1 nM/min

IPTG_th=[0,logspace(-1, 4)];
    
 i=1:51;
    [t,outcome1]=EMB_1_circuit_decay(IPTG_th(i));
    [t,outcome2]=EMB_2_circuit_decay(IPTG_th(i));
    [t,outcome3]=EMB_3_circuit_decay(IPTG_th(i));
    y_R_1{i,1}=outcome1(:,:);
    y_R_2{i,1}=outcome2(:,:);
    y_R_3{i,1}=outcome3(:,:);
    YFP_th(i,1)=y_R_1{i,1}(end,8);
    YFP_th(i,2)=y_R_2{i,1}(end,8);
    YFP_th(i,3)=y_R_3{i,1}(end,7);
    RFP_th(i,1)=y_R_1{i,1}(end,5);
    RFP_th(i,2)=y_R_2{i,1}(end,5);
    RFP_th(i,3)=y_R_3{i,1}(end,5);


        [t,outcome]= EMB_1_circuit_decay(a) % and other according to the topology
    nspecies = 8; % IPTG, free LacI, IPTG-LacI, RFP, RFPmat, TetR, YFP, YFPmat
    alphamax = 100e-6*60; % uM/min expression factor, around 100pM/s and up, from Schwarz-Schilling, Aufinger et al.

    aLacI = 3*1000/1077; % rNTPs length LacI, factor of expression normalized, and plasmid ratio
    aRFP = 3*1000/693; % rNTPs length RFP
    aTetR = 1000/618; % rNTPs length TetR
    aYFP = 1000/717; % rNTPs length YFP

    tRNA = 15; % min, lifetime of RNA, from Karzbrun et al.
    tmat1 = 33; % min, maturation time of mScarletI, from Balleza et al. https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5765880/pdf/nihms916443.pdf
    tmat2 = 12; % min, maturation time of YPet, from Balleza et al. https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5765880/pdf/nihms916443.pdf
    t1 = 150; % min, delay before cell-extract degrades
    t2 = 170; % min, characteristic time of cell-extract degradation

    KdLacI = 3.6e-2; % uM, LacI binding to pLacO, determined from bulk fitting
    nLacI = 1.35; % Hill coefficient, LacI binding to pLacO, determined from bulk fitting
    KdTetR = 1.3e-1; % uM, determined from bulk fitting
    nTetR = 4.3; % Hill coefficient, TetR binding to pTetO, determined from bulk fitting

    kon = 7.2; % /uM/min, on-rate of IPTG binding LacI, from Xu et al.
    koff = 12.6; % /min, off-rate of IPTG binding LacI, from Xu et al.

    kdeg = 1e-3; % uM/min, rate of 0th-order protein degradation, from Karzbrun et al.
    Kdeg = 1e-3; % uM, deducted from Karzbrun et al.



    %initialization
    for ispecies=1:1:nspecies;
        x0(ispecies, 1)=0;
    end

    %initial values < > 0
    x0(1,1)=a; %M, IPTG sender

    % duration of experiment
    dt=10;
    duration=600; % min, duration of experiment: 8h

    % simulation run
    tspan=[0:dt:duration];

    options=odeset('AbsTol',1e-9,'RelTol',1e-5); % set tolerances
    [t,y]=ode23s(@EMB_1_circuit_decay_function,tspan,x0,options);

    %output
    t=t;
    outcome=y(:,:);

        ydot = EMB_1_circuit_decay_function(t, y, flags)
    dy(1,1)=-kon*y(1,1)*y(2,1)+koff*y(3,1); % IPTG

    if t < t1; 
        dy(2,1)=alphamax*aLacI*(1-exp(-t/tRNA))*(1./(1+(y(6,1)./KdTetR).^nTetR))-kdeg/(1+Kdeg/y(2,1))-kon*y(1,1)*y(2,1)+koff*y(3,1); % LacI
        dy(3,1)=kon*y(1,1)*y(2,1)-koff*y(3,1)-kdeg/(1+Kdeg/y(3,1)); % IPTG-LacI
        dy(4,1)=alphamax*aRFP*(1-exp(-t/tRNA))*(1./(1+(y(6,1)./KdTetR).^nTetR))-y(4,1)/tmat1-kdeg/(1+Kdeg/y(4,1)); % RFP
        dy(5,1)=y(4,1)/tmat1-kdeg/(1+Kdeg/y(5,1)); % RFPmat
        dy(6,1)=alphamax*aTetR*(1-exp(-t/tRNA))*(1./(1+(y(2,1)./KdLacI).^nLacI)); % TetR
        dy(7,1)=alphamax*aYFP*(1-exp(-t/tRNA))*(1./(1+(y(2,1)./KdLacI).^nLacI))-y(7,1)/tmat2; % YFP
    else;
        dy(2,1)=(alphamax*aLacI*(1-exp(-t/tRNA))*(1./(1+(y(6,1)./KdTetR).^nTetR))-kdeg/(1+Kdeg/y(2,1)))*exp(-(t-t1)/t2)-kon*y(1,1)*y(2,1)+koff*y(3,1);
        dy(3,1)=kon*y(1,1)*y(2,1)-koff*y(3,1)-kdeg/(1+Kdeg/y(3,1))*exp(-(t-t1)/t2); 
        dy(4,1)=(alphamax*aRFP*(1-exp(-t/tRNA))*(1./(1+(y(6,1)./KdTetR).^nTetR))-kdeg/(1+Kdeg/y(4,1)))*exp(-(t-t1)/t2)-y(4,1)/tmat1;
        dy(5,1)=y(4,1)/tmat1-kdeg/(1+Kdeg/y(5,1))*exp(-(t-t1)/t2); 
        dy(6,1)=alphamax*aTetR*(1-exp(-t/tRNA))*(1./(1+(y(2,1)./KdLacI).^nLacI))*exp(-(t-t1)/t2);
        dy(7,1)=alphamax*aYFP*(1-exp(-t/tRNA))*(1./(1+(y(2,1)./KdLacI).^nLacI))*exp(-(t-t1)/t2)-y(7,1)/tmat2;
    end

    dy(8,1)=y(7,1)/tmat2; % YFPmat

    ydot=dy;


plot(IPTG_th, YFP_th);
hold on; axis square;
plot(IPTG_th, RFP_th);
legend('boxoff');


%% Figure 2: Redo with more data points and compare with XP.
IPTG_th=[0, logspace(-1, 4, 100)];

for i=1:101;
    [t,outcome1]=EMB_1_circuit_decay(IPTG_th(i));
    [t,outcome2]=EMB_2_circuit_decay(IPTG_th(i));
    [t,outcome3]=EMB_3_circuit_decay(IPTG_th(i));
    y_R_1{i,1}=outcome1(:,:);
    y_R_2{i,1}=outcome2(:,:);
    y_R_3{i,1}=outcome3(:,:);
    YFP_th(i,1)=y_R_1{i,1}(end,8);
    YFP_th(i,2)=y_R_2{i,1}(end,8);
    YFP_th(i,3)=y_R_3{i,1}(end,7);
    RFP_th(i,1)=y_R_1{i,1}(end,5);
    RFP_th(i,2)=y_R_2{i,1}(end,5);
    RFP_th(i,3)=y_R_3{i,1}(end,5);
end

YFP_factor=YFP_th(101,1)/YFP_exp(6,1);
RFP_factor=RFP_th(1,2)/RFP_exp(1,2);

IPTG_th_plot=[1e-2, logspace(-1,4,100)]; 
subplot(2,3,1); hold on; axis square; plot(IPTG_exp, YFP_exp(:,1)*YFP_factor*1e-6); plot(IPTG_th_plot*1e-6, YFP_th(:,1)*1e-6);
subplot(2,3,2); hold on; axis square; plot(IPTG_exp, YFP_exp(:,2)*YFP_factor*1e-6); plot(IPTG_th_plot*1e-6, YFP_th(:,2)*1e-6);
subplot(2,3,3); hold on; axis square; plot(IPTG_exp, YFP_exp(:,3)*YFP_factor*1e-6); plot(IPTG_th_plot*1e-6, YFP_th(:,3)*1e-6);
subplot(2,3,4); hold on; axis square; plot(IPTG_exp, RFP_exp(:,1)*RFP_factor*1e-6); plot(IPTG_th_plot*1e-6, RFP_th(:,1)*1e-6);
subplot(2,3,5); hold on; axis square; plot(IPTG_exp, RFP_exp(:,2)*RFP_factor*1e-6); plot(IPTG_th_plot*1e-6, RFP_th(:,2)*1e-6);
subplot(2,3,6); hold on; axis square; plot(IPTG_exp, RFP_exp(:,3)*RFP_factor*1e-6); plot(IPTG_th_plot*1e-6, RFP_th(:,3)*1e-6);

%% Figure 3: finding the hill curves of the simulations and plotting them together

% with the data from Figure 1, fitted with y=a*(x^n/(Kd^n+x^n))+b (YFP) or y=a*(Kd^n/(Kd^n+x^n));




