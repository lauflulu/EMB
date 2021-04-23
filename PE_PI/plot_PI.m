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
    

%% PI DIR/ SGA, extrapolation for t=91
% a B*=5 is a valid choice for estimating single gene PI for
% all data sets except 100 mM full (N=9) and 1 mM ctrl1 (N=12)
for s=1:length(datafiles)
    load(datafiles{s});

    fusionTime=91;
    g=EMB_data2g(data,fusionTime);
    g=g(:,:,:,91);
    [N,I,X,T]=size(g);
    
    if N>12
        m=flip(unique(ceil((6:11)/12*N)));
    else
        m=flip(7:N);
    end
    % number of iterations for each (bin,m)
    K=100;
    [PIdir,stdPIdir] = EMB_extrapolatePIv1(@EMB_g2piDIR,g,K,m,2:6,false,false);
    
    [PIsga,stdPIsga] = EMB_extrapolatePIv1(@EMB_g2piSGA,g,K,m,100,false,false);
    
    figure(11)
        subplot(1,4,4)
        hold all
        set(gca,'ColorOrderIndex',1)
        for p=1:3
            errorbar(PIsga(p,:),PIdir(p,:),...
                stdPIdir(p,:),stdPIdir(p,:),stdPIsga(p,:),stdPIsga(p,:),...
                '.-','MarkerSize',20)
            text(PIsga(p,:),PIdir(p,:),sprintf('%d',s))
        end
        plot(linspace(0,1.2),linspace(0,1.2),'--k')
        box('on'); xlim([-0.1,1.3]);ylim([-0.1,1.3]);
        xlabel('I_{SGA}');ylabel('I_{DIR}')
        legend('g1','g2','joint')
        
        subplot(1,4,1:3)
        hold all
        
        for p=1:3
            set(gca,'ColorOrderIndex',p)
            errorbar(s-0.1,PIsga(p,:),stdPIsga(p,:),'.-','MarkerSize',20)
            set(gca,'ColorOrderIndex',p)
            errorbar(s+0.1,PIdir(p,:),stdPIdir(p,:),'.-','MarkerSize',20)
        end
        box('on'); xlim([0.5,10.5]);ylim([-0.1,1.3]);
        xlabel('Data set number');ylabel('I-DIR, I-SGA')
        legend('g1','g2','joint')
end

%% PI DIR/ SGA, time traces
% a B*=5 is a valid choice for estimating single gene PI for
% all data sets except 100 mM full (N=9) and 1 mM ctrl1 (N=12)
time=0:90; time=time*5/60;

for s=3%1:length(datafiles)
    load(datafiles{s});

    fusionTime=91;
    g=EMB_data2g(data,fusionTime);
    [N,I,X,T]=size(g);
    
    if N>12
        m=flip(unique(ceil((6:11)/12*N)));
    else
        m=flip(7:N);
    end
    % number of iterations for each (bin,m)
    K=100;
    [PIdir,stdPIdir] = EMB_extrapolatePIv1(@EMB_g2piDIR,g,K,m,2:6,false,false);
    
    [PIsga,stdPIsga] = EMB_extrapolatePIv1(@EMB_g2piSGA,g,K,m,100,false,false);
    
    figure(11)
        subplot(1,2,1)
            hold all
            plot(time, PIdir, '-');set(gca,'ColorOrderIndex',1)
            plot(time, PIdir+stdPIdir, '--');set(gca,'ColorOrderIndex',1)
            plot(time, PIdir-stdPIdir, '--');
            box('on'); xlim([0,7.5]);ylim([-0.1,1.2]);
            xlabel('Time (h)');ylabel('I_{DIR}')
            legend('g1','g2','joint')
            xticks([0:2.5:7.5]);
            
        subplot(1,2,2)
            hold all
            plot(time, PIsga{3,1}, '-');set(gca,'ColorOrderIndex',1)
            plot(time, PIsga{3,1}+stdPIsga{3,1}, '--');set(gca,'ColorOrderIndex',1)
            plot(time, PIsga{3,1}-stdPIsga{3,1}, '--');
            box('on'); xlim([0,7.5]);ylim([-0.1,1.2]);
            xlabel('Time (h)');ylabel('I_{SGA}')
            legend('g1','g2','joint');xticks([0:2.5:7.5]);
end

%% PI SGA, time traces
time=0:90; time=time*5/60;
tpoints=1:91;

PIsga=cell(S,1);
stdPIsga=cell(S,1);
PIdir=cell(S,1);
stdPIdir=cell(S,1);

tic;
for s=1:S
    load(datafiles{s});

    fusionTime=91;
    g=EMB_data2g(data,fusionTime);
    g=g(:,:,:,tpoints);
    [N,I,X,T]=size(g);
    
    if N>12
        m=flip(unique(ceil((6:11)/12*N)));
    else
        m=flip(7:N);
    end
    % number of iterations for each (bin,m)
    K=100;
    %[PIsga{s,1},stdPIsga{s,1}] = EMB_g2piSGA(g,100);
    [PIsga{s,1},stdPIsga{s,1}] = EMB_extrapolatePIv1(@EMB_g2piSGA,g,K,m,100,false,false);
    %[PIsga{s,1},stdPIsga{s,1}] = EMB_extrapolatePIv3(@EMB_g2piSGA,g,K,m,100,false,true);
    [PIdir{s,1},stdPIdir{s,1}] = EMB_extrapolatePIv1(@EMB_g2piDIR,g,K,m,2:6,false,false);
    s
    toc
end
toc
% plots SGA
figure(11)
for s=1:S
    subplot(5,2,s)
        hold all
        plot(time, PIsga{s,1}, '-');set(gca,'ColorOrderIndex',1)
        plot(time, PIsga{s,1}+stdPIsga{s,1}, '--');set(gca,'ColorOrderIndex',1)
        plot(time, PIsga{s,1}-stdPIsga{s,1}, '--');
        box('on'); xlim([0,7.5]);ylim([-0.1,1.5]);
        xlabel('Time (h)');ylabel('I_{SGA}')
        legend('g1','g2','joint');xticks([0:2.5:7.5]);
        title(sprintf('N = %d',N));
end
% plots DIR
figure(12)
for s=1:S
    subplot(5,2,s)
        hold all
        plot(time, PIdir{s,1}, '-');set(gca,'ColorOrderIndex',1)
        plot(time, PIdir{s,1}+stdPIdir{s,1}, '--');set(gca,'ColorOrderIndex',1)
        plot(time, PIdir{s,1}-stdPIdir{s,1}, '--');
        box('on'); xlim([0,7.5]);ylim([-0.1,1.5]);
        xlabel('Time (h)');ylabel('I_{DIR}')
        legend('g1','g2','joint');xticks([0:2.5:7.5]);
        title(sprintf('N = %d',N));
end
%% endpoint PI comparison
 %load('analyze_data\PIdir_t91_200430.mat')
 %load('analyze_data\PIsga_200430.mat')
 t=91;
 PIsga_end=zeros(S,3);
 stdPIsga_end=zeros(S,3);
  PIdir_end=zeros(S,3);
 stdPIdir_end=zeros(S,3);
 for s=1:S
     PIsga_end(s,:)=PIsga{s,1}(:,t);
     stdPIsga_end(s,:)=stdPIsga{s,1}(:,t);
     PIdir_end(s,:)=PIdir{s,1}(:,t);
     stdPIdir_end(s,:)=stdPIdir{s,1}(:,t);
 end
 xpoints=[1:10]'*[1,1,1];
 
 figure(20)
    hold all
    errorbar(xpoints-0.1,PIsga_end,stdPIsga_end,'--.');
    errorbar(xpoints+0.1,PIdir_end,stdPIdir_end,'--.')
    plot([1,10],[0,0],'-k')
    xlim([0.5,10.5]); ylim([-0.1,1.5]);
    xlabel('Data Set Number');
    ylabel('PI_{SGA}(7.5 h) (bits)')
    box('on');legend();
    
 figure(22)
    hold all
    errorbar(PIsga_end,PIdir_end,...
        stdPIdir_end,stdPIdir_end,...
        stdPIsga_end,stdPIsga_end,'.','MarkerSize',20)
    plot([0,1.5],[0,1.5],'-k')
    xlim([-0.1,1.5]); ylim([-0.1,1.5]);
    xlabel('PI_{SGA}(7.5 h) (bits)');
    ylabel('PI_{DIR}(7.5 h) (bits)')
    box('on');legend();
    
%%
S=size(PIsga,1);
 t=91;
 PIsga_end=zeros(S,3);
 stdPIsga_end=zeros(S,3);
  PIdir_end=zeros(S,3);
 stdPIdir_end=zeros(S,3);
 color=[[0,0,1];...
     [1,0,0];...
     [0.8,0.8,0.8]];
 for s=1:S
     PIsga_end(s,:)=PIsga{s,1}(:,t);
     stdPIsga_end(s,:)=stdPIsga{s,1}(:,t);
     PIdir_end(s,:)=PIdir{s,1}(:,t);
     stdPIdir_end(s,:)=stdPIdir{s,1}(:,t);
 end
 xpoints=flip(logspace(-1,2,4)'*[1,1,1]);
 figure(21)
    hold all
    for i=1:3
        errorbar(xpoints(1:4,i),PIsga_end(1:4,i),stdPIsga_end(1:4,i),'--.','Color',color(i,:)*3/3)
        errorbar(xpoints(1:4,i),PIsga_end(5:8,i),stdPIsga_end(5:8,i),'--.','Color',color(i,:)*2/3)
        errorbar(xpoints(1:2,i),PIsga_end(9:10,i),stdPIsga_end(9:10,i),'--.','Color',color(i,:)*1/3)
    end
    plot([0.1,100],[0,0],'-k')
    xlim([0.1*0.8,100*1.2]); ylim([-0.1,1.5]);
    xlabel('IPTG (mM)');
    set(gca,'Xscale','log')
    ylabel('PI_{SGA}(7.5 h) (bits)')
    box('on')