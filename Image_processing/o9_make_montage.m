% to generate montage for overview of all datasets
clear all
close all
clc

%%
T=[1,49,49,49,49]; %frame of the montage
datafiles={'190925_full_100mM\data_full_100mM.mat',...
        '190927_full_10mM\data_full_10mM.mat',...
        '190928_full_1mM\data_full_1mM.mat',...
        '190929_full_100uM\data_full_100uM.mat',...
        '190930_ctrl1_100mM\data_ctrl1_100mM.mat',...
        '190930_ctrl1_10mM\data_ctrl1_10mM.mat',...
        '191002_ctrl1_1mM\data_ctrl1_1mM.mat',...
        '191002_ctrl1_100uM\data_ctrl1_100uM.mat',...
        '191010_ctrl2_100mM\data_ctrl2_100mM.mat',...
        '191010_ctrl2_10mM\data_ctrl2_10mM.mat'};

%% colormaps
mymap = cell(5,1);
mymap{1,1}=gray(255);
mymap{2,1}=gray(255);
mymap{3,1}=[zeros(255,1)  linspace(0,1,255)'  linspace(0,1,255)'];
mymap{4,1}=[linspace(0,1,255)'  zeros(255,1)  zeros(255,1)];
mymap{5,1}=gray(255);
%% main loop
for p=1:size(datafiles,2)
    load(datafiles{p})
    folder=split(datafiles{p},'\');
    fileList=dir([folder{1,1},'/o6-final_videos/*.tif']);

    Nfiles=length(fileList);
    nameList=cell(Nfiles,1);
    channelList=zeros(Nfiles,1);
    positionList=cell(Nfiles,1);
    for i=1:Nfiles
        filename=fileList(i).name;
        nameList{i,1}=filename;
        temp = strsplit(filename,'-');
        channel=regexp(temp{2},'\d*','match');
        channelList(i,1)=str2double(channel{1});
        positionList{i,1} = temp{1};
    end

    channelList=unique(channelList);

    %% sort position list
    positionList=unique(positionList); % all cropped images
    positionList=positionList';
    b={data.name}; % all makeData images (T>10)
    dS={data.dataSize}; 
    dS=squeeze(cat(3, dS{:}));
    dS4h=(min(dS,[],1)>49);
    noFusion4h=b(1,dS4h);
    dS7h=(min(dS,[],1)>90);
    noFusion7_5h=b(1,dS7h);
    newPosList=[noFusion7_5h,...
        noFusion4h(~ismember(noFusion4h,noFusion7_5h)),...
        b(~ismember(b,noFusion4h)),...
        positionList(~ismember(positionList,b))];
    positionList=newPosList';
    %% color of box with name identifier, corresponds to sort quality
    boxColor = [ones(length(noFusion7_5h),1)*[0 0.8 0];...
        ones(length(noFusion4h)-length(noFusion7_5h),1)*[0.8 0.8 0];...
        ones(length(b)-length(noFusion4h),1)*[0.8 0.4 0];...
        ones(length(positionList)-length(b),1)*[0.8 0 0]];

    %% load images
    Nchannels = length(channelList);
    Npositions= length(positionList);
    Ncolumns = 5; %number of columns in montage, user input
    Nrows=ceil(Npositions/Ncolumns);


    channelNames=cell(size(positionList));
    images=cell(2,5);
    for i=1:Npositions
        for c=1:Nchannels
            channelNames{i,c}=[folder{1,1},'/o6-final_videos/',positionList{i,1},sprintf('-C%d.tif',c-1)];
            Gray=mat2gray(imread(channelNames{i,c}, T(c)));
            GrayIndex = uint8(floor(Gray * 255));
            images{i,c}      = ind2rgb(GrayIndex, mymap{c,1});
        end
        % mark circles

        id=positionList{i,1};
        idx=find(ismember(b,id));
        if idx
            X=data(idx).x(T(1),:)'; X=X(~isnan(X));
            Y=data(idx).y(T(1),:)'; Y=Y(~isnan(Y));
            R=data(idx).r(T(1),:)'; R=R(~isnan(R));
            D=data(idx).d(T(1),:)'; D=D(~isnan(D));

            if R
%                 images{i,1}=insertShape(images{i,1},'circle',...
%                      [X,Y,R],'LineWidth',2);
                images{i,1}=insertText(images{i,1},[X(D==0),Y(D==0)],'S',...
                 'FontSize',30,'BoxColor','red','BoxOpacity',0,...
                 'TextColor','black','AnchorPoint','Center');
            end
        end

    end
    images=images';        
    %% make montages
    Npages=ceil(Npositions/15);

    for i=1:Npages
        range=(i-1)*15+1:i*15;
        range=range(range<=Npositions);
        Lr=length(range);

            mont=montage(images(:,range),'Size', [Lr Nchannels],'ThumbnailSize',[300,600]);
            [X,Y,~]=size(mont.CData);

            textPos=[zeros(Lr,1),...
                linspace(0,X/Lr*(Lr-1),Lr)'];

            mont=insertText(mont.CData,textPos(1:Lr,:),positionList(range,1),...
                    'FontSize',30,'BoxColor',...
                    boxColor(range,:),'BoxOpacity',1,'TextColor','white');

            imwrite(mont,[folder{1,1},sprintf('_%d.png',i)]);
            %figure(i)
            %imshow(mont);
    end

end