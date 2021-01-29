%improved background subtraction
close all
clear all

tic;
main_dir=cd;
dirList=dir(main_dir);
Ndir=length(dirList);
folderList=cell(10,1);
count=0;
for i=1:Ndir
    if dirList(i).isdir==1 && ~isempty(regexp(dirList(i).name,'\d*','match')) %&& dirList(i).name ~= '..'
        count=count+1;
        folderList{count,1}=dirList(i).name;
    end
end

%% big loop
total_tic=tic;
sub_tic=tic;
total_wait=waitbar(0,sprintf('Elapsed: %.0f s,\nEstimated to finish: ',toc(total_tic)));
sub_wait=waitbar(0,[sprintf('Elapsed: %.0f s,\nEstimated to finish: ',toc(sub_tic)),datestr(now)]);

for d=3%1:Ndir
    sub_tic=tic;
    newdir=[main_dir ,'\', folderList{d,1}];
    cd(newdir);
    mkdir('o7-background_videos/');
    fileList=dir('o6-final_videos/*.tif');

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
    positionList=unique(positionList);

    Nchannels = length(channelList);
    Npositions= length(positionList);
    Nframes = size(imfinfo(['o4-cropped_videos/',positionList{1,1},'-C2.tif']),1);

    %% load image
    count=0;
    
    for p=32%1:Npositions % position
        count=count+1;
        for c = 4%2:4 % channel
            for t = 49%:Nframes % timepoint

                filename=[positionList{p,1},sprintf('-C%d.tif',c)];
                flatfile=[positionList{p,1},sprintf('-F%d.tif',c)];
                excludefile=[positionList{p,1},'-B1.tif'];
                
                I0=imread(['o4-cropped_videos/',filename], t);
                F0=imread(['o6-final_videos/',flatfile]);
                B0=imread(['o7-exclude_videos/',excludefile], t);

                x0=meshgrid(1:size(I0,2),1:size(I0,1));
                y0=meshgrid(1:size(I0,1),1:size(I0,2))';

                scale=0.3;
                x_down=imresize(x0,scale);
                y_down=imresize(y0,scale);
                I_down=imresize(I0,scale);
                B_down=imresize(B0,scale);

                x=reshape(x_down,[],1);
                y=reshape(y_down,[],1);
                I=double(reshape(I_down,[],1));
                B=[reshape(B_down,[],1)==0];

                % the fit
                ft = fittype( 'poly22' );
                fo = fitoptions( 'Method', 'LinearLeastSquares','Exclude',B);
                fo.Robust = 'Bisquare';

                % Fit model to data.
                [sf, gof] = fit( [x, y], I, ft, fo);
                
                % correct image
                y0_v=reshape(y0,[],1);
                x0_v=reshape(x0,[],1);
                I0_v=double(reshape(I0,[],1));
                B0_v=[reshape(B0,[],1)==0];
                
                I_corr=single(double(I0_v)-sf(x0_v,y0_v));
                I_corr=reshape(I_corr,size(I0,1),size(I0,2));
                z1=sf(x,y);
                z=reshape(z1,size(I_down,1),size(I_down,2));
                Z=uint16(reshape(sf(x0_v,y0_v),size(I0,1),size(I0,2)));
                I_corr=uint16(I_corr);
                %I_corr=uint16(I_corr./F0);
                imwrite(Z, ['o7-background_videos/',filename],'writemode','append','Compression','none');

                imwrite(I_corr, ['o7-background_videos/',filename],'writemode','append','Compression','none');

            end
        end    
        waitbar(count/Npositions,sub_wait,[sprintf('Elapsed: %.0f s,\nEstimated to finish: ',toc(sub_tic)),datestr(now+toc(sub_tic)/count*(Npositions-count)/24/3600)]);
    end
    waitbar(d/Ndir,total_wait,sprintf('Elapsed: %.0f s',toc(total_tic)));
end

%% display
II=reshape(I_corr,1,[])';
I00=reshape(I0_v,1,[])';

z0=reshape(sf(x0_v,y0_v),size(I0,1),size(I0,2));




figure(2)
        
        subplot(1,2,1)
        hold all
        surf(x0,y0,z0,'edgecolor','none');

        C=(double(I0(~(B0==0)))-min(double(I0(~(B0==0)))))./(max(double(I0(~(B0==0))))-min(double(I0(~(B0==0)))))*0.8;
        scatter3(x0(~(B0==0)),y0(~(B0==0)),I0(~(B0==0)),...
            [],[C C C],'.')
        box('on')
         xlabel('x (px)')
        ylabel('y (px)')
        zlabel('gray value')
        set(gca,'Ydir','reverse')
        view(-15,30)
        xlim([0,800]); ylim([0,400]); zlim([600,1600]);

    subplot(1,2,2)
        hold all
        surf(x0,y0,z0,'edgecolor','none');
        C=(double(I0(~(B0==0)))-min(double(I0(~(B0==0)))))./(max(double(I0(~(B0==0))))-min(double(I0(~(B0==0)))))*0.8;
        scatter3(x0(~(B0==0)),y0(~(B0==0)),I0(~(B0==0)),...
            [],[C C C],'.')
        box('on')
         xlabel('x (px)')
        ylabel('y (px)')
        zlabel('gray value')
        set(gca,'Ydir','reverse')
        view(15,-30)
        xlim([0,800]); ylim([0,400]); zlim([600,1600]);