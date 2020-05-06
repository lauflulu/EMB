clear all
close all

%%
datafiles={'data\data_full_100mM_v2.mat',...
        'data\data_full_10mM_v2.mat',...
        'data\data_full_1mM_v2.mat',...
        'data\data_full_100uM_v2.mat',...
        'data\data_ctrl1_100mM_v2.mat',...
        'data\data_ctrl1_10mM_v2.mat',...
        'data\data_ctrl1_1mM_v2.mat',...
        'data\data_ctrl1_100uM_v2.mat',...
        'data\data_ctrl2_100mM_v2.mat',...
        'data\data_ctrl2_10mM_v2.mat'};
    
S=length(datafiles);
    
%% PI SGA, time traces
time=0:90; time=time*5/60;
tpoints=1:20;

PIsga=cell(S,1);
stdPIsga=cell(S,1);
tic;
for s=1%:S
    load(datafiles{s});

    fusionTime=91;
    g=EMB_data2g(data,fusionTime);
    g=g(:,:,:,tpoints);
    [N,I,X,T]=size(g);
    
    if N*0.5>8
        m=flip(unique(ceil(N*0.5:N/12:N-1)));
    else
        m=flip(7:N);
    end
    % number of iterations for each (bin,m)
    K=100;
    %[PIsga{s,1},stdPIsga{s,1}] = EMB_g2piSGA(g,100);
    [PIsga{s,1},stdPIsga{s,1}] = EMB_extrapolatePI(@EMB_g2piSGA,g,K,m,100,false,true);
    %[PIsga{s,1},stdPIsga{s,1}] = EMB_extrapolatePIv3(@EMB_g2piSGA,g,K,m,100,false,true);
    %[PIdir{s,1},stdPIdir{s,1}] = EMB_extrapolatePIv3(@EMB_g2piDIR,g,K,m,2:6,false,true);
end
toc
%%
for s=1%:S
    
    load(datafiles{s});

    fusionTime=91;
    g=EMB_data2g(data,fusionTime);
    g=g(:,:,:,tpoints);
    [N,I,X,T]=size(g);
    figure(11)
    
    
        %subplot(5,2,s)
            hold all
            plot(time(tpoints), PIsga{s,1}, '-');set(gca,'ColorOrderIndex',1)
            plot(time(tpoints), PIsga{s,1}+stdPIsga{s,1}, '--');set(gca,'ColorOrderIndex',1)
            plot(time(tpoints), PIsga{s,1}-stdPIsga{s,1}, '--');
            box('on'); xlim([0,7.5]);ylim([-0.1,1.5]);
            xlabel('Time (h)');ylabel('I_{SGA}')
            legend('g1','g2','joint');xticks([0:2.5:7.5]);
            title(sprintf('N = %d',N));
end
%% endpoint PI comparison
 load('analyze_data\PIdir_t91_200430.mat')
 load('analyze_data\PIsga_200430.mat')
 t=91;
 PIsga_end=zeros(S,3);
 stdPIsga_end=zeros(S,3);
  PIdir_end=zeros(S,3);
 stdPIdir_end=zeros(S,3);
 for s=1:S
     PIsga_end(s,:)=PIsga{s,1}(:,t);
     stdPIsga_end(s,:)=stdPIsga{s,1}(:,t);
     PIdir_end(s,:)=PIdir{s,1}(:,1);
     stdPIdir_end(s,:)=stdPIdir{s,1}(:,1);
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
 figure(21)
    hold all
    errorbar(xpoints(1:4,1:2),PIsga_end(1:4,1:2),stdPIsga_end(1:4,1:2),'--.')
    errorbar(xpoints(1:3,1:2),PIsga_end(5:7,1:2),stdPIsga_end(5:7,1:2),'--.')
    errorbar(xpoints(4,1:2),PIsga_end(8,1:2),stdPIsga_end(8,1:2),'--.')
    errorbar(xpoints(1:2,1:2),PIsga_end(9:10,1:2),stdPIsga_end(9:10,1:2),'--.')
    plot([1,4],[0,0],'-k')
    xlim([0.5,4.5]); ylim([-0.1,1.5]);
    xlabel('Data Set Number');
    ylabel('PI_{SGA}(7.5 h) (bits)')
    box('on')