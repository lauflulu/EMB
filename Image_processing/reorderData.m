
filename=dir('data*.mat');
load(filename.name);

for i=1:size(data,2)
senderID=findSenderIndex2(data,i);
data(i).d=distanceFromSender(data,i,senderID);
data = sortByDistance(data,i ); %could be replaced by index list
end
fname=string(filename.name);
save(fname,'data');