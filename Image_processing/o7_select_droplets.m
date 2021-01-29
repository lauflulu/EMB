%%
clear all
close all

%%
% filter parameters
minDataSize=10;
% minCyclicity=0.8;
% minRadius=20;
% maxRadius=30;

T=5; % time interval in minutes

dropletList=dir('o8-droplets\droplets*.mat');
data=struct;
c=0;
for i=1:size(dropletList,1)
    load(['o8-droplets\',dropletList(i).name])
    filteredDroplets=droplets.select('dataSize',minDataSize);
    % filteredDroplets=filteredDroplets.select('meanCyclicity',minCyclicity);
    % filteredDroplets=filteredDroplets.select('meanRadius',minRadius,'meanRadius',maxRadius);
    
    % manual selection
    if size(filteredDroplets,2)>0
        waitfor(filteredDroplets.dialog);
        filteredDroplets=dropletSelection_1;
    end
    if c==0
        oldDroplets=filteredDroplets;
    end
    
    % write data. only for complete samples and avoid duplicates in case
    % selection is empty
    if size(filteredDroplets,2)==6 && (filteredDroplets(1).p(1,1)~=oldDroplets(1).p(1,1) || c==0)
        c=c+1;
        droplets=filteredDroplets;
        oldDroplets=droplets;
        position=[droplets.p];
        intSum=[droplets.intensitySum];

        data(c).x=position(:,1:2:size(position,2));
        data(c).y=position(:,2:2:size(position,2));
        data(c).r=[droplets.radius];
        data(c).YFP=intSum(:,2:4:size(intSum,2));
        data(c).RFP=intSum(:,3:4:size(intSum,2));
        data(c).Cy5=intSum(:,4:4:size(intSum,2));

        % compute distance and sort
        senderID=findSenderIndex2(data,c);
        data(c).d=distanceFromSender(data,c,senderID);
        
        % filter criteria
        % data(c).fusions=countFusions(data,c);
        % data(c).receivers=numberReceivers(data,c);
        data(c).dataSize=getDataSize(data,c);

        % compute additional measures
        % data(c).hDIB=calculateDIBsize(data,c);
        data(c).N=countNeighbors2(data,c,20);
        
        % other info
        data(c).t=linspace(0,(size(position,1)-1)*5,size(position,1));
        name=droplets.stacks; name=name{1}.file; name=strsplit(name,'\');
        name=strsplit(name{end},'-');
        data(c).name=name{1};
        
        data = sortByDistance(data,c);
    end
end