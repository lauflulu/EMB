clear all
close all

%%
datafiles={'raw_data\data_full_100mM_v2.mat',...
        'raw_data\data_full_10mM_v2.mat',...
        'raw_data\data_full_1mM_v2.mat',...
        'raw_data\data_full_100uM_v2.mat',...
        'raw_data\data_ctrl1_100mM_v2.mat',...
        'raw_data\data_ctrl1_10mM_v2.mat',...
        'raw_data\data_ctrl1_1mM_v2.mat',...
        'raw_data\data_ctrl1_100uM_v2.mat',...
        'raw_data\data_ctrl2_100mM_v2.mat',...
        'raw_data\data_ctrl2_10mM_v2.mat'};

%% PI DIR I_shuffle verification
% because of statistical fluctuations: 
%verify automatically computed B* estimates by looking at the graphs
% a B*=5 is a valid choice for estimating single gene PI for
% all data sets except 100 mM full (N=9) and 1 mM ctrl1 (N=11)
S=length(datafiles);

Bstar=zeros(S,3);
PI=cell(S,1);
stdPI=cell(S,1);

fusionTime=91;
K=100;
tol=0.1;
binNumber=2:20;
B=length(binNumber);

shuffle=true;

for s=1:length(datafiles)
    s
  
    load(datafiles{s});
    g=EMB_data2g(data,fusionTime);
    g=g(:,:,:,fusionTime);
    [N,I,X,T]=size(g);
    
    if N>12
        m=flip(unique(ceil((6:11)/12*N)));
    else
        m=flip(7:N);
    end
    
    PIb=zeros(B-2,3);
    stdPIb=zeros(B-2,3);

    for i=3:B
        [PIb(i-2,:),stdPIb(i-2,:)] = EMB_extrapolatePI(@EMB_g2piDIR,g,K,m,binNumber(1:i),shuffle,false);
    end
    PI{s,1}=PIb; stdPI{s,1}=stdPIb;
    
    figure(3)
        subplot(5,2,s)
        hold all
        errorbar(binNumber(3:B)'*ones(1,3),PI{s,1},stdPI{s,1},'--o')
        plot(binNumber,zeros(1,B),'--k')
        plot(binNumber,0.1*ones(1,B),'--k')
        ylim([-0.1,1.7]); box('on'); xlim([0,21]);
        xlabel('B*');ylabel('PI (bits)');
        title(sprintf('Data set %.0f',s))
    
      
    bb=binNumber(3:B)'*ones(1,3);
    Bmax=bb.*((PIb(:,:)+stdPIb(:,:))>tol); 
    Bmax(Bmax==0)=nan;
    Bmax=min(Bmax,[],1)-1;
    Bmax(isnan(Bmax))=21;
    Bmax(Bmax==3)=nan;
    Bstar(s,:)=Bmax;
end