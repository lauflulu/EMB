clear all
close all

tracking = DropletTracking(); tracking.open();
%selectFolder
%select dummy position
%segment dummy position
%configure tracker
%remove dummy position
%run script
%track
%%
currentDirectory=cd;
tiffDirectory=[currentDirectory,'\o6-final_videos\'];
cd(tiffDirectory);
tiffList=dir('*C0.tif');
NC=4; %=BF+Fluorescence, not binary
%
for p=1:size(tiffList,1)
    for i=0:NC
        tiffName=tiffList(p).name;
        tiffName(end-4)=sprintf('%d',i);
        tracking.positions{1,p}.stacks{1,i+1}=TiffStack([tiffDirectory,tiffName]);
    end
    p
end

cd(currentDirectory);
mkdir('o8-droplets')