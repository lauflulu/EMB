%% 20-05-03 Figures

%% Figure 1: comparing the M/min and uM/min calculations on topology 1, with 1st order degradation but 0th order rate value

IPTG=[0, 1, 10e3];
kdeg=[15e-5, 15e-4, 15e-3];

for i=1:3;
    for j=1:3;
        [t,outcome]=EMB_1_circuit_decay1(IPTG(j),kdeg(i));
        RFP_1{j,i}=outcome(:,5);
        YFP_1{j,i}=outcome(:,8);
    end
end
% computing time: 9 s


        [t,outcome]= EMB_1_circuit_decay1(a,b)
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
    koff = 12.6; %/min, off-rate of IPTG binding LacI, from Xu et al.

    kdeg = b; % 15e-3 uM/min, rate of 0th-order protein degradation, from Karzbrun et al.


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
    [t,y]=ode23s(@EMB_1_circuit_decay1_function,tspan,x0,options);

    %output
    t=t;
    outcome=y(:,:);


        ydot = EMB_1_circuit_decay1_function(t, y, flags)

    dy(1,1)=-kon*y(1,1)*y(2,1)+koff*y(3,1); % IPTG

    if t < t1; 
        dy(2,1)=alphamax*aLacI*(1-exp(-t/tRNA))*(1./(1+(y(6,1)./KdTetR).^nTetR))-kdeg*y(2,1)-kon*y(1,1)*y(2,1)+koff*y(3,1); % LacI
        dy(4,1)=alphamax*aRFP*(1-exp(-t/tRNA))*(1./(1+(y(6,1)./KdTetR).^nTetR))-kdeg*y(4,1)-y(4,1)/tmat1; % RFP
        dy(6,1)=alphamax*aTetR*(1-exp(-t/tRNA))*(1./(1+(y(2,1)./KdLacI).^nLacI)); % TetR
        dy(7,1)=alphamax*aYFP*(1-exp(-t/tRNA))*(1./(1+(y(2,1)./KdLacI).^nLacI))-y(7,1)/tmat2; % YFP
    else;
        dy(2,1)=alphamax*aLacI*(1-exp(-t/tRNA))*(1./(1+(y(6,1)./KdTetR).^nTetR))*exp(-(t-t1)/t2)-kdeg*y(2,1)-kon*y(1,1)*y(2,1)+koff*y(3,1);
        dy(4,1)=alphamax*aRFP*(1-exp(-t/tRNA))*(1./(1+(y(6,1)./KdTetR).^nTetR))*exp(-(t-t1)/t2)-kdeg*y(4,1)-y(4,1)/tmat1;
        dy(6,1)=alphamax*aTetR*(1-exp(-t/tRNA))*(1./(1+(y(2,1)./KdLacI).^nLacI))*exp(-(t-t1)/t2);
        dy(7,1)=alphamax*aYFP*(1-exp(-t/tRNA))*(1./(1+(y(2,1)./KdLacI).^nLacI))*exp(-(t-t1)/t2)-y(7,1)/tmat2;
    end

    dy(3,1)=kon*y(1,1)*y(2,1)-koff*y(3,1)-kdeg*y(3,1); % IPTG-LacI

    dy(5,1)=y(4,1)/tmat1-kdeg*y(5,1); % RFPmat

    dy(8,1)=y(7,1)/tmat2; % YFPmat

    ydot=dy;


for i=1:3;
    [t,outcome]=EMB_1_circuit_decay(IPTG(i));
    RFP_ref{i,1}=outcome(:,5);
    YFP_ref{i,1}=outcome(:,8);
end
% computing time 0.8s


        [t,outcome]= EMB_1_circuit_decay(a)
    nspecies = 8; % IPTG, free LacI, IPTG-LacI, RFP, RFPmat, TetR, YFP, YFPmat

    % alphamax = 0.0036; % expression factor, determined from bulk and droplet fitting
    alphamax = 100e-12*60; % M/min expression factor, around 100pM/s and up, from Schwarz-Schilling, Aufinger et al.

    aLacI = 3*1000/1077; % rNTPs length LacI, factor of expression normalized, and plasmid ratio
    aRFP = 3*1000/693; % rNTPs length RFP
    aTetR = 1000/618; % rNTPs length TetR
    aYFP = 1000/717; % rNTPs length YFP

    tRNA = 15; % min, lifetime of RNA, from Karzbrun et al.
    tmat1 = 33; % min, maturation time of mScarletI, from Balleza et al. https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5765880/pdf/nihms916443.pdf
    tmat2 = 12; % min, maturation time of YPet, from Balleza et al. https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5765880/pdf/nihms916443.pdf
    t1 = 150; % min, delay before cell-extract degrades
    t2 = 170; % min, characteristic time of cell-extract degradation

    KdLacI = 3.6e-8; % M, LacI binding to pLacO, determined from bulk fitting
    nLacI = 1.35; % Hill coefficient, LacI binding to pLacO, determined from bulk fitting
    KdTetR = 1.3e-7; % M, determined from bulk fitting
    nTetR = 4.3; % Hill coefficient, TetR binding to pTetO, determined from bulk fitting

    kon = 7.2e6; % /M/min, on-rate of IPTG binding LacI, from Xu et al.
    koff = 12.6; % /min, off-rate of IPTG binding LacI, from Xu et al.

    kdeg = 15e-9; % M/min, rate of 0th-order protein degradation, from Karzbrun et al.


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
        dy(2,1)=alphamax*aLacI*(1-exp(-t/tRNA))*(1./(1+(y(6,1)./KdTetR).^nTetR))-kdeg*y(2,1)-kon*y(1,1)*y(2,1)+koff*y(3,1); % LacI
        dy(4,1)=alphamax*aRFP*(1-exp(-t/tRNA))*(1./(1+(y(6,1)./KdTetR).^nTetR))-y(4,1)/tmat1-kdeg*y(4,1); % RFP
        dy(6,1)=alphamax*aTetR*(1-exp(-t/tRNA))*(1./(1+(y(2,1)./KdLacI).^nLacI)); % TetR
        dy(7,1)=alphamax*aYFP*(1-exp(-t/tRNA))*(1./(1+(y(2,1)./KdLacI).^nLacI))-y(7,1)/tmat2; % YFP
    else;
        dy(2,1)=alphamax*aLacI*(1-exp(-t/tRNA))*(1./(1+(y(6,1)./KdTetR).^nTetR))*exp(-(t-t1)/t2)-kdeg*y(2,1)-kon*y(1,1)*y(2,1)+koff*y(3,1);
        dy(4,1)=alphamax*aRFP*(1-exp(-t/tRNA))*(1./(1+(y(6,1)./KdTetR).^nTetR))*exp(-(t-t1)/t2)-y(4,1)/tmat1-kdeg*y(4,1);
        dy(6,1)=alphamax*aTetR*(1-exp(-t/tRNA))*(1./(1+(y(2,1)./KdLacI).^nLacI))*exp(-(t-t1)/t2);
        dy(7,1)=alphamax*aYFP*(1-exp(-t/tRNA))*(1./(1+(y(2,1)./KdLacI).^nLacI))*exp(-(t-t1)/t2)-y(7,1)/tmat2;
    end

    dy(3,1)=kon*y(1,1)*y(2,1)-koff*y(3,1)-kdeg*y(3,1); % IPTG-LacI

    dy(5,1)=y(4,1)/tmat1-kdeg*y(5,1); % RFPmat

    dy(8,1)=y(7,1)/tmat2; % YFPmat

    ydot=dy;

subplot(2,2,1); hold on; axis square;
plot(t, [YFP_ref{1,1}, YFP_ref{2,1}, YFP_ref{3,1}]*1e6);
subplot(2,2,2); hold on; axis square; plot(t, [RFP_ref{1,1}, RFP_ref{2,1}, RFP_ref{3,1}]*1e6);
subplot(2,2,3); hold on; axis square; plot(t, [YFP_1{1,3}, YFP_1{2,3}, YFP_1{3,3}]);
subplot(2,2,4); hold on; axis square; plot(t, [RFP_1{1,3}, RFP_1{2,3}, RFP_1{3,3}]);
legend('boxoff');

%% Figure 2: trying the constant decay and the transfer function decay

IPTG=[0, 1, 10e3];
kdeg=[15e-5, 15e-4, 15e-3];
for i=1:3;
    for j=1:3;
        [t,outcome]=EMB_1_circuit_decay2(IPTG(j),kdeg(i));
        RFP_2{j,i}=outcome(:,5);
        YFP_2{j,i}=outcome(:,8);
    end
end

% this took 30 min so far and still bugged.. I stopped it

for i=1:3;
    for j=1:3;
        [t,outcome]=EMB_1_circuit_decay3(IPTG(j),kdeg(i));
        RFP_3{j,i}=outcome(:,5);
        YFP_3{j,i}=outcome(:,8);
    end
end
% 45 s!

        [t,outcome]= EMB_1_circuit_decay3(a,b)
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
    koff = 12.6; %/min, off-rate of IPTG binding LacI, from Xu et al.

    kdeg = b; % 15e-3 uM/min, rate of 0th-order protein degradation, from Karzbrun et al.
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
    [t,y]=ode23s(@EMB_1_circuit_decay3_function,tspan,x0,options);

    %output
    t=t;
    outcome=y(:,:);



        ydot = EMB_1_circuit_decay3_function(t, y, flags)
    dy(1,1)=-kon*y(1,1)*y(2,1)+koff*y(3,1); % IPTG

    if t < t1; 
        dy(2,1)=alphamax*aLacI*(1-exp(-t/tRNA))*(1./(1+(y(6,1)./KdTetR).^nTetR))-kdeg*(1/(1+Kdeg/y(2,1)))-kon*y(1,1)*y(2,1)+koff*y(3,1); % LacI
        dy(4,1)=alphamax*aRFP*(1-exp(-t/tRNA))*(1./(1+(y(6,1)./KdTetR).^nTetR))-kdeg*(1/(1+Kdeg/y(4,1)))-y(4,1)/tmat1; % RFP
        dy(6,1)=alphamax*aTetR*(1-exp(-t/tRNA))*(1./(1+(y(2,1)./KdLacI).^nLacI)); % TetR
        dy(7,1)=alphamax*aYFP*(1-exp(-t/tRNA))*(1./(1+(y(2,1)./KdLacI).^nLacI))-y(7,1)/tmat2; % YFP
    else;
        dy(2,1)=alphamax*aLacI*(1-exp(-t/tRNA))*(1./(1+(y(6,1)./KdTetR).^nTetR))*exp(-(t-t1)/t2)-kdeg*(1/(1+Kdeg/y(2,1)))-kon*y(1,1)*y(2,1)+koff*y(3,1);
        dy(4,1)=alphamax*aRFP*(1-exp(-t/tRNA))*(1./(1+(y(6,1)./KdTetR).^nTetR))*exp(-(t-t1)/t2)-kdeg*(1/(1+Kdeg/y(4,1)))-y(4,1)/tmat1;
        dy(6,1)=alphamax*aTetR*(1-exp(-t/tRNA))*(1./(1+(y(2,1)./KdLacI).^nLacI))*exp(-(t-t1)/t2);
        dy(7,1)=alphamax*aYFP*(1-exp(-t/tRNA))*(1./(1+(y(2,1)./KdLacI).^nLacI))*exp(-(t-t1)/t2)-y(7,1)/tmat2;
    end

    dy(3,1)=kon*y(1,1)*y(2,1)-koff*y(3,1)-kdeg*(1/(1+Kdeg/y(3,1))); % IPTG-LacI

    dy(5,1)=y(4,1)/tmat1-kdeg*(1/(1+Kdeg/y(5,1))); % RFPmat

    dy(8,1)=y(7,1)/tmat2; % YFPmat

    ydot=dy;

subplot(1,2,1); hold on; axis square
for i=1:3;
    for j=1:3;
        plot(t/60, YFP_3{j,i});
    end
end
subplot(1,2,2); hold on; axis square;
for i=1:3;
    for j=1:3;
        plot(t/60, RFP_3{j,i});
    end
end

% 15e-3 uM/min seems really too strong to have RFP win in the absence of
% IPTG. Let's go back to the data to see if we witness clear degradation. 

%% Figure 3

kdeg=[5e-3,10e-3,15e-3];

for i=1:3;
    for j=1:3;
        [t,outcome]=EMB_1_circuit_decay3(IPTG(j),kdeg(i));
        RFP_3{j,i}=outcome(:,5);
        YFP_3{j,i}=outcome(:,8);
    end
end
% 21 s


subplot(1,2,1); hold on; axis square;
for i=1:3, 
    for j=1:3;
        plot(t/60, YFP_3{j,i});
    end
end
subplot(1,2,2); hold on; axis square;
for i=1:3, 
    for j=1:3;
        plot(t/60, RFP_3{j,i});
    end
end

% a degradation rate of 15nM / min is really huge, especially with such a
% low expression (6 nM / min max). 

%% Figure 4: changing expression rate

alphamax=[100e-6*60,200e-6*60,300e-6*60];
IPTG=[0, 1, 10e3];

for i=1:3;
    for j=1:3;
        [t,outcome]=EMB_1_circuit_decay3(IPTG(j),alphamax(i));
        RFP_3{j,i}=outcome(:,5);
        YFP_3{j,i}=outcome(:,8);
    end
end

        [t,outcome]= EMB_1_circuit_decay3(a,b)
    nspecies = 8; % IPTG, free LacI, IPTG-LacI, RFP, RFPmat, TetR, YFP, YFPmat

    alphamax = b; % 100e-6*60 uM/min expression factor, around 100pM/s and up, from Schwarz-Schilling, Aufinger et al. 

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
    koff = 12.6; %/min, off-rate of IPTG binding LacI, from Xu et al.

    kdeg = 15e-3; % uM/min, rate of 0th-order protein degradation, from Karzbrun et al.
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
    [t,y]=ode23s(@EMB_1_circuit_decay3_function,tspan,x0,options);

    %output
    t=t;
    outcome=y(:,:);


        ydot = EMB_1_circuit_decay3_function(t, y, flags)
    dy(1,1)=-kon*y(1,1)*y(2,1)+koff*y(3,1); % IPTG

    if t < t1; 
        dy(2,1)=alphamax*aLacI*(1-exp(-t/tRNA))*(1./(1+(y(6,1)./KdTetR).^nTetR))-kdeg*(1/(1+Kdeg/y(2,1)))-kon*y(1,1)*y(2,1)+koff*y(3,1); % LacI
        dy(4,1)=alphamax*aRFP*(1-exp(-t/tRNA))*(1./(1+(y(6,1)./KdTetR).^nTetR))-kdeg*(1/(1+Kdeg/y(4,1)))-y(4,1)/tmat1; % RFP
        dy(6,1)=alphamax*aTetR*(1-exp(-t/tRNA))*(1./(1+(y(2,1)./KdLacI).^nLacI)); % TetR
        dy(7,1)=alphamax*aYFP*(1-exp(-t/tRNA))*(1./(1+(y(2,1)./KdLacI).^nLacI))-y(7,1)/tmat2; % YFP
    else;
        dy(2,1)=alphamax*aLacI*(1-exp(-t/tRNA))*(1./(1+(y(6,1)./KdTetR).^nTetR))*exp(-(t-t1)/t2)-kdeg*(1/(1+Kdeg/y(2,1)))-kon*y(1,1)*y(2,1)+koff*y(3,1);
        dy(4,1)=alphamax*aRFP*(1-exp(-t/tRNA))*(1./(1+(y(6,1)./KdTetR).^nTetR))*exp(-(t-t1)/t2)-kdeg*(1/(1+Kdeg/y(4,1)))-y(4,1)/tmat1;
        dy(6,1)=alphamax*aTetR*(1-exp(-t/tRNA))*(1./(1+(y(2,1)./KdLacI).^nLacI))*exp(-(t-t1)/t2);
        dy(7,1)=alphamax*aYFP*(1-exp(-t/tRNA))*(1./(1+(y(2,1)./KdLacI).^nLacI))*exp(-(t-t1)/t2)-y(7,1)/tmat2;
    end

    dy(3,1)=kon*y(1,1)*y(2,1)-koff*y(3,1)-kdeg*(1/(1+Kdeg/y(3,1))); % IPTG-LacI

    dy(5,1)=y(4,1)/tmat1-kdeg*(1/(1+Kdeg/y(5,1))); % RFPmat

    dy(8,1)=y(7,1)/tmat2; % YFPmat

    ydot=dy;


subplot(1,2,1); hold on; axis square
for i=1:3;
for j=1:3;
plot(t/60, YFP_3{j,i});
end
end
subplot(1,2,2); hold on; axis square;
for i=1:3;
for j=1:3;
plot(t/60, RFP_3{j,i});
end
end
% improving the expression rate does not help: it affects everything
% homogeneously while degradation only affects LacI/RFP, so it makes
% TetR/YFP win.

%% Figure 5: screening kdeg again

IPTG=[0, 1, 10e3];
kdeg=[5e-4, 1e-3, 5e-3];

for i=1:3;
    for j=1:3;
        [t,outcome]=EMB_1_circuit_decay3(IPTG(j),kdeg(i));
        RFP_3{j,i}=outcome(:,5);
        YFP_3{j,i}=outcome(:,8);
    end
end


        [t,outcome]= EMB_1_circuit_decay3(a,b)

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
    koff = 12.6; %/min, off-rate of IPTG binding LacI, from Xu et al.

    kdeg = b; % 15e-3 uM/min, rate of 0th-order protein degradation, from Karzbrun et al.
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
    [t,y]=ode23s(@EMB_1_circuit_decay3_function,tspan,x0,options);

    %output
    t=t;
    outcome=y(:,:);


        ydot = EMB_1_circuit_decay3_function(t, y, flags)
    dy(1,1)=-kon*y(1,1)*y(2,1)+koff*y(3,1); % IPTG

    if t < t1; 
        dy(2,1)=alphamax*aLacI*(1-exp(-t/tRNA))*(1./(1+(y(6,1)./KdTetR).^nTetR))-kdeg*(1/(1+Kdeg/y(2,1)))-kon*y(1,1)*y(2,1)+koff*y(3,1); % LacI
        dy(4,1)=alphamax*aRFP*(1-exp(-t/tRNA))*(1./(1+(y(6,1)./KdTetR).^nTetR))-kdeg*(1/(1+Kdeg/y(4,1)))-y(4,1)/tmat1; % RFP
        dy(6,1)=alphamax*aTetR*(1-exp(-t/tRNA))*(1./(1+(y(2,1)./KdLacI).^nLacI)); % TetR
        dy(7,1)=alphamax*aYFP*(1-exp(-t/tRNA))*(1./(1+(y(2,1)./KdLacI).^nLacI))-y(7,1)/tmat2; % YFP
    else;
        dy(2,1)=alphamax*aLacI*(1-exp(-t/tRNA))*(1./(1+(y(6,1)./KdTetR).^nTetR))*exp(-(t-t1)/t2)-kdeg*(1/(1+Kdeg/y(2,1)))-kon*y(1,1)*y(2,1)+koff*y(3,1);
        dy(4,1)=alphamax*aRFP*(1-exp(-t/tRNA))*(1./(1+(y(6,1)./KdTetR).^nTetR))*exp(-(t-t1)/t2)-kdeg*(1/(1+Kdeg/y(4,1)))-y(4,1)/tmat1;
        dy(6,1)=alphamax*aTetR*(1-exp(-t/tRNA))*(1./(1+(y(2,1)./KdLacI).^nLacI))*exp(-(t-t1)/t2);
        dy(7,1)=alphamax*aYFP*(1-exp(-t/tRNA))*(1./(1+(y(2,1)./KdLacI).^nLacI))*exp(-(t-t1)/t2)-y(7,1)/tmat2;
    end

    dy(3,1)=kon*y(1,1)*y(2,1)-koff*y(3,1)-kdeg*(1/(1+Kdeg/y(3,1))); % IPTG-LacI

    dy(5,1)=y(4,1)/tmat1-kdeg*(1/(1+Kdeg/y(5,1))); % RFPmat

    dy(8,1)=y(7,1)/tmat2; % YFPmat

    ydot=dy;

%% Figure 6: doing a more systematic screening of the transfer function to understand which kdeg is reasonable

%notes: Looking at Karzbrun paper, in previous work it's been kdeg=ClpX*1
%nM/min. We assume 1nM ClpX, so 1 nM/min is reasonable. 
%other works show a highest kdeg of around 10 nM/min, but
%https://www.nature.com/articles/s41598-018-21739-6.pdf shows no
%degradation almost for low protein concentrations (below 1uM), so deg rate
%can be very low. With ClpX they have 1uM/30min being degraded, so 33
%nM/min for 100nM ClpX. With 1 nM ClpX, we would have 0.33 nM/min

IPTG=logspace(-2, 5, 20);
IPTG=IPTG';
kdeg=linspace(0.1,15,20)*1e-3;

for i=1:20;
    for j=1:20;
        [t,outcome]=EMB_1_circuit_decay3(IPTG(j),kdeg(i));
        transfer_3{1,1}(j,i)=outcome(end,8);
        transfer_3{1,2}(j,i)=outcome(end,5);
    end
end

        [t,outcome]= EMB_1_circuit_decay3(a,b)
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
    koff = 12.6; %/min, off-rate of IPTG binding LacI, from Xu et al.

    kdeg = b; % 15e-3 uM/min, rate of 0th-order protein degradation, from Karzbrun et al.
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
    [t,y]=ode23s(@EMB_1_circuit_decay3_function,tspan,x0,options);

    %output
    t=t;
    outcome=y(:,:);


        ydot = EMB_1_circuit_decay3_function(t, y, flags)
    dy(1,1)=-kon*y(1,1)*y(2,1)+koff*y(3,1); % IPTG

    if t < t1; 
        dy(2,1)=alphamax*aLacI*(1-exp(-t/tRNA))*(1./(1+(y(6,1)./KdTetR).^nTetR))-kdeg*(1/(1+Kdeg/y(2,1)))-kon*y(1,1)*y(2,1)+koff*y(3,1); % LacI
        dy(3,1)=kon*y(1,1)*y(2,1)-koff*y(3,1)-kdeg*(1/(1+Kdeg/y(3,1))); % IPTG-LacI
        dy(4,1)=alphamax*aRFP*(1-exp(-t/tRNA))*(1./(1+(y(6,1)./KdTetR).^nTetR))-kdeg*(1/(1+Kdeg/y(4,1)))-y(4,1)/tmat1; % RFP
        dy(5,1)=y(4,1)/tmat1-kdeg*(1/(1+Kdeg/y(5,1))); % RFPmat
        dy(6,1)=alphamax*aTetR*(1-exp(-t/tRNA))*(1./(1+(y(2,1)./KdLacI).^nLacI)); % TetR
        dy(7,1)=alphamax*aYFP*(1-exp(-t/tRNA))*(1./(1+(y(2,1)./KdLacI).^nLacI))-y(7,1)/tmat2; % YFP
    else;
        dy(2,1)=(alphamax*aLacI*(1-exp(-t/tRNA))*(1./(1+(y(6,1)./KdTetR).^nTetR))-kdeg*(1/(1+Kdeg/y(2,1))))*exp(-(t-t1)/t2)-kon*y(1,1)*y(2,1)+koff*y(3,1);
        dy(3,1)=kon*y(1,1)*y(2,1)-koff*y(3,1)-kdeg*(1/(1+Kdeg/y(3,1)))*exp(-(t-t1)/t2); % IPTG-LacI
        dy(4,1)=(alphamax*aRFP*(1-exp(-t/tRNA))*(1./(1+(y(6,1)./KdTetR).^nTetR))-kdeg*(1/(1+Kdeg/y(4,1))))*exp(-(t-t1)/t2)-y(4,1)/tmat1;
        dy(5,1)=y(4,1)/tmat1-kdeg*(1/(1+Kdeg/y(5,1)))*exp(-(t-t1)/t2); % RFPmat
        dy(6,1)=alphamax*aTetR*(1-exp(-t/tRNA))*(1./(1+(y(2,1)./KdLacI).^nLacI))*exp(-(t-t1)/t2);
        dy(7,1)=alphamax*aYFP*(1-exp(-t/tRNA))*(1./(1+(y(2,1)./KdLacI).^nLacI))*exp(-(t-t1)/t2)-y(7,1)/tmat2;
    end

    dy(8,1)=y(7,1)/tmat2; % YFPmat

    ydot=dy;


for i=1:20;
    plot(IPTG, transfer_3{1,1}(:,i));
end

for i=1:20;
    plot(IPTG, transfer_3{1,2}(:,i));
end


