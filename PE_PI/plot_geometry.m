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


%% droplet and bilayer size
fusionTime=91;
s=3;
load(datafiles{s});

time=data(1).t/60; %time=time(1,1:fusionTime);
[Rb,Rd,Ab,Vd,rl,rr,d] = EMB_data2geometry(data,fusionTime);
[T,X,N]=size(d);

meanAb=mean(Ab,2:3,'omitnan');
stdAb=std(Ab,[],2:3,'omitnan');

meanVd=mean(Vd,2:3,'omitnan');
stdVd=std(Vd,[],2:3,'omitnan');
CVVd=stdVd./meanVd;

meanVd_intra=mean(Vd,2,'omitnan');
stdVd_intra=std(Vd,[],2,'omitnan');
CVVd_intra=mean(stdVd_intra./meanVd_intra,3,'omitnan');

meanVd_inter=mean(Vd,2,'omitnan');
stdVd_inter=std(meanVd_inter,[],3,'omitnan');
CVVd_inter=stdVd_inter./mean(meanVd_inter,3,'omitnan');

theta=(asin(rl)+asin(rr))/2/pi*180;
theta2=180-acos(-(d.^2-(Rb.*rl.^(-1)).^2-(Rb.*rr.^(-1)).^2)./(Rb.*rr.^(-1))./(Rb.*rl.^(-1)))/pi*180;
D=Rb.*((rl.^(-2)-1).^0.5 + (rr.^(-2)-1).^0.5);
%%
figure(7)     
    subplot(3,2,1)
        hold all
        plot(time,meanAb,'-r','LineWidth',3)
        plot(time,meanAb+stdAb,'-b','LineWidth',3)
        plot(time,meanAb-stdAb,'-b','LineWidth',3)
        box('on');xlabel('Time (h)');ylabel('DIB Area (?m^2)');
    
    subplot(3,2,2)
        hold all
        plot(time,meanVd,'-r','LineWidth',3)
        plot(time,meanVd+stdVd,'-b','LineWidth',3)
        plot(time,meanVd-stdVd,'-b','LineWidth',3)
        box('on');xlabel('Time (h)');ylabel('Droplet Volume (?m^3)');
        
    subplot(3,2,3)
        hold all
        plot(squeeze(rl(1,:,:)),squeeze(rr(1,:,:)),'ob')
        plot(squeeze(rl(13,:,:)),squeeze(rr(49,:,:)),'or')
        plot(squeeze(rl(49,:,:)),squeeze(rr(91,:,:)),'ok')
        box('on');xlabel('Rb/Rl');ylabel('Rb/Rr');
    subplot(3,2,4)
        hold all
        plot(squeeze(Rb(1,:,:)),squeeze(theta(1,:,:)),'ob');
        plot(squeeze(Rb(13,:,:)),squeeze(theta(13,:,:)),'or');
        plot(squeeze(Rb(49,:,:)),squeeze(theta(49,:,:)),'ok');
        box('on');xlabel('Rb');ylabel('theta');
    subplot(3,2,5)
        hold all
        plot(squeeze(d(1,:,:)),squeeze(D(1,:,:)),'ob')
        plot(squeeze(d(13,:,:)),squeeze(D(49,:,:)),'or')
        plot(squeeze(d(49,:,:)),squeeze(D(91,:,:)),'ok')
        plot(0:140,0:140,'-g','LineWidth',3)
        box('on');xlabel('d');ylabel('D');
    subplot(3,2,6)
        hold all
        plot(squeeze(theta(1,:,:)),squeeze(theta2(1,:,:)),'ob');
        plot(squeeze(theta(13,:,:)),squeeze(theta2(13,:,:)),'or');
        plot(squeeze(theta(49,:,:)),squeeze(theta2(49,:,:)),'ok');
        plot(0:180,0:180,'-g','LineWidth',3)
        box('on');xlabel('Rb');ylabel('theta2');
        
        
%% plots
time=data(1).t/60;
figure(1)
    subplot(2,2,1)
        hold all
        for n=1:N
            p=scatter(reshape(time'*ones(1,5),[],1),...
                reshape(data(n).DIB,[],1),'filled','MarkerFaceColor',[0 0 0]);
            set(p,'MarkerFaceAlpha',.02);
        end
        Mdib=nanmean(cat(2,data.DIB),2);
        STDdib=nanstd(cat(2,data.DIB),[],2);
        plot(time,Mdib,'-r','LineWidth',3)
        plot(time,Mdib+STDdib,'-b','LineWidth',3)
        plot(time,Mdib-STDdib,'-b','LineWidth',3)
        box('on');xlabel('Time (h)');ylabel('DIB Length (px)');
    
    subplot(2,2,2)
        hold all
        for n=1:N
            p=scatter(reshape(time'*ones(1,6),[],1),...
                reshape(data(n).area,[],1),'filled','MarkerFaceColor',[0 0 0]);
            set(p,'MarkerFaceAlpha',.02);
        end
        Marea=nanmean(cat(2,data.area),2);
        STDarea=nanstd(cat(2,data.area),[],2);
        plot(time,Marea,'-r','LineWidth',3)
        plot(time,Marea+STDarea,'-b','LineWidth',3)
        plot(time,Marea-STDarea,'-b','LineWidth',3)
        box('on');xlabel('Time (h)');ylabel('Droplet Area (px)');
        
    subplot(2,2,3)
        hold all
        dibAmicron=(cat(2,data.DIB)*1.3/2).^2*pi;
        MdibArea=nanmean(dibAmicron,2);
        STDdibArea=nanstd(dibAmicron,[],2);
        plot(time,MdibArea,'-r','LineWidth',3)
        plot(time,MdibArea+STDdibArea,'-b','LineWidth',3)
        plot(time,MdibArea-STDdibArea,'-b','LineWidth',3)
        box('on');xlabel('Time (h)');ylabel('DIB Area (?m^2)');
    
    subplot(2,2,4)
        hold all
        dropRmicron=(cat(2,data.area)./pi).^0.5*1.3;
        dropVmicron=4/3*pi*dropRmicron.^3;
        Mvol=nanmean(dropVmicron,2);
        STDvol=nanstd(dropVmicron,[],2);
        plot(time,Mvol,'-r','LineWidth',3)
        plot(time,Mvol+STDvol,'-b','LineWidth',3)
        plot(time,Mvol-STDvol,'-b','LineWidth',3)
        box('on');xlabel('Time (h)');ylabel('Droplet Volume (?m^3)');