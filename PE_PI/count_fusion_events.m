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
    
S=length(datafiles);
%%
fused=zeros(S,1);
unfused=zeros(S,1);

fusionTime=91;

for s=1:S
    load(datafiles{s});
    
    a={data.r}; a=cat(3, a{:});
    for n=1:length(data)
        if sum(sum(isnan(a(1:fusionTime,:,n))))==0
            unfused(s,1)=unfused(s,1)+1;
        else
            fused(s,1)=fused(s,1)+1;
        end
    end
end

sum(unfused)/(sum(unfused)+sum(fused))
