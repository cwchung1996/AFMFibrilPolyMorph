%% @cwchung1996/AFMFibrilPolyMorph
% Requires 
% HeightSensor.m and the MATLAB Bruker toolkit to be installed beforehand
% ImageDisp.m and scatt.m for creating figures

% Housekeeping
tic;clc;clear;close all;

%Image dimensions (please alter accordingly)
n=[5,2]; %number of images per dataset
umpx=[2 2 2 2 2;
    2 2 0 0 0]; %field of view
px512=512; %image dimension
umpx=umpx/px512; %length (um) per pixel
Dataname=["aS_Day3";"Abeta"]; %Image name to display on graphs
Filename=Dataname; %Image name
Legname=["WT (Day 3)";"A\beta"]; %Image name to appear on legend of figures

% Paramter for thresholding image (please alter accordingly)
thres=[0.18,0.22];
pmax_dd=2000;

% Parameters for finding peaks and troughs (please alter accordingly)
minlength=[420,500];
PeakProminence=[0.2,0.2];
MinSep=[10,1];

%Dummy variables (do not touch)
maxval=2000;
xdist_pk=zeros(maxval,maxval*5);
xdist_trghs=zeros(maxval,maxval*5);
hdist_trghs=zeros(maxval,maxval*5);
xdist=zeros(maxval,maxval*5);
pdist=zeros(maxval,maxval*5);
hdist_pk=zeros(maxval,maxval*5);
hdist=zeros(maxval,maxval*5);
height_save=zeros(maxval,maxval*5);
TF=zeros(maxval*3,maxval*5);
colourmat=zeros(maxval*5,3);
count_save=0;
count=1;
count_fiblength=1;

%% Analysis
for datanum=1:length(Dataname) %For loop for different datasets
    for ndata=1:n(datanum) %For loop for different images within a dataset
        name=sprintf('%s_Test%d.spm',Dataname(datanum),ndata);
        name=fullfile(pwd,name);
        dataout=HeightSensor(name);
        corr=min(dataout,[],'all');
        fig=image((dataout),'CDataMapping','scaled');
        axis('tight', 'square');
        colormap(pink);
        colorbar;
        title(sprintf('%s: %s (Image %d)','Height Sensor',Legname(datanum),ndata));
        axis('tight', 'square');
        drawnow;
        imgname=sprintf('%s%s%d','HeightProf_',Filename(datanum),ndata);
        saveas(fig,imgname,'tif');
        
        %% Masking based on intensity
        datalogi=logical(dataout);
        pmin_dd=thres(datanum)*max(dataout,[],'all');
        mask_P_dd = (dataout>pmin_dd)&(dataout<pmax_dd);
        dataMatrix2_dd = datalogi.*mask_P_dd;
        
        [clusters, nclus]=bwlabeln(dataMatrix2_dd);
        c1=1;c2=length(clusters);
        skelimg=zeros(size(dataout));
        for i=1:nclus
            colourmat(count,:)=rand(1,3);
            [r c]=find(clusters==i);
            if length(r)<minlength(datanum)
                
            else
                clustersedit=clusters;
                idx=clustersedit~=i;
                clustersedit(idx) = 0 ;
                idx=clustersedit==i;
                maskclus(c1:c2,c1:c2)=dataout.*(idx);
                raw_data2_clus=single(maskclus(c1:c2,c1:c2));
                
                %% Skeletonising
                clustersedit = bwmorph(clustersedit,'bridge');
                clustersedit = bwmorph(clustersedit,'clean');
                clustersedit = bwmorph(clustersedit,'close');
                clustersedit = bwmorph(clustersedit,'fill');
                clustersedit=bwmorph(clustersedit, 'thin', Inf) ;
                clustersedit = bwskel(logical(clustersedit),'MinBranchLength',50);
                skelimg=clustersedit+skelimg;
                figure;imshow(clustersedit);
                
                endpointImage = bwmorph(clustersedit,'endpoints');
                
                
                [idx,idy]=find(clustersedit==1);
                co=[idx,idy];
                
                [rows, columns] = find(endpointImage);
                lr=length(rows);
                if lr>2 || lr==0 || lr==1
                    count=count-1;
                else
                    
                    xi(1)=rows(2);yi(1)=columns(2); %endpoint as initial guess
                    xi_save=xi(1);
                    yi_save=yi(1);
                    
                    idx_find=find(co(:,1)==xi_save);
                    idy_find=find(co(:,2)==yi_save);
                    
                    [val,pos]=intersect(idx_find,idy_find);
                    
                    co_save(1,:)=[idx(val),idy(val)];
                    co(val,:)=[];
                    
                    dist_save(1)=0;
                    height(1)=dataout(idx(val),idy(val));
                    % [valy,posy]=intersect(idy_find,idx_find)  %unhash for
                    % checking purposes
                    
                    for iii = 2:length(idx)
                        [k,dist] = dsearchn(co,co_save(iii-1,:));
                        dist_save(iii)=dist+dist_save(iii-1);
                        co_save(iii,:)=co(k,:);
                        x=co(k,1);y=co(k,2);
                        height(iii)=dataout(x,y);
                        co(k,:)=[];

                        
                    end
                    
                    
                    height_save(:,count)=[height,zeros(1,maxval-length(height))]';
                    
                    
                    
                    
                    %% Finding peaks
                    height=smooth(height);height=height';
                    TFmax = islocalmax(height,'MinProminence',PeakProminence(datanum),'MinSeparation',MinSep(datanum),'FlatSelection', 'first');
                    pks=height(TFmax);
                    locspks=dist_save(TFmax);
                    
                    TFmin = islocalmin(height,'MinProminence',PeakProminence(datanum),'MinSeparation',MinSep(datanum),'FlatSelection', 'first');
                    trghs=height(TFmin);
                    locstrghs=dist_save(TFmin);
                    
                    xdist(:,count)=[dist_save*umpx(datanum,ndata),zeros(1,maxval-length(dist_save))]';
                    xdist_pk(:,count)=[locspks*umpx(datanum,ndata),zeros(1,maxval-length(locspks))]';
                    xdist_trghs(:,count)=[locstrghs*umpx(datanum,ndata),zeros(1,maxval-length(locstrghs))]';
                    
                    %% Finding fibril length
                    fibdx0=find(co_save(:,1)==0);
                    fibdx511=find(co_save(:,1)==px512-1);
                    fibdy0=find(co_save(:,2)==0);
                    fibdy511=find(co_save(:,2)==px512-1);
                    if isempty(fibdx0)==1 && isempty(fibdx511)==1 && isempty(fibdy0)==1 && isempty(fibdy511)==1
                        fiblength(:,count_fiblength)=dist_save(end)*umpx(datanum,ndata);
                        count_fiblength=count_fiblength+1;
                    else
                        
                    end
                    
                    
                    
                    if isempty(xdist)==1
                        
                    else
                        hdist(:,count)=[height,zeros(1,maxval-length(height))]';
                        hdist_pk(:,count)=[pks,zeros(1,maxval-length(pks))]';
                        hdist_trghs(:,count)=[trghs,zeros(1,maxval-length(trghs))]';
                        
                        for k=2:length(nonzeros(hdist_pk(:,count)))
                            pdist(k-1,count)=xdist_pk(k,count)-xdist_pk(k-1,count);
                        end
                        scatx(count)=i;
                    end
                    
                end
                
                clear raw_data2_clus
                count=count+1;
                clear dist_save height pks locs co_save co dist
                
                
            end
            %Making skeletal traces
            skelimg=logical(skelimg);
            imwrite(skelimg,fullfile(pwd, sprintf('Skel_%s%d.tif',Filename(datanum),ndata)));
            
            
        end
        clear clusters clustersedit r c
        
    end
    count_save(datanum+1)=count-1;
    count_fiblength_save(datanum+1)=count_fiblength-1;
end
hold off;

%end

%end


close all;


%% Figure for height profile
%maxplot=10;

Filenum1=1;
Filenum=0;

for q=1:length(count_save)-1
    %Plot all data in a matrix
    x1=count_save(q)+1;
    x2=count_save(q+1);
    Filenum=Filenum+1;
    Filenum1=1;
    figure;
    suptitle(Legname(q));
    count_plot=1;
    subpl=round(sqrt(count_save(q+1)));subpl=subpl(end);
    for i=x1:x2
        
        
        subplot(6,5,count_plot);
        
        plotx1=nonzeros(xdist(:,i));plotx1=[0;plotx1];
        ploty1=nonzeros(hdist(:,i));
        if isempty(nonzeros(hdist(:,i)))
            subplot(8,8,count_plot);
            plot(plotx1,ploty1,'color',colourmat(i,:));hold on
            
        else
            plot(plotx1,ploty1,'color',colourmat(i,:));hold on
            plot(nonzeros(xdist_pk(:,i)),nonzeros(hdist_pk(:,i)),'ro');
            plot(nonzeros(xdist_trghs(:,i)),nonzeros(hdist_trghs(:,i)),'r*');
            xlabel('Fibril length [\mum]');
            ylabel('H [nm]');
            Titlename=sprintf('Fibril %d',Filenum1);

        end
        count_plot=count_plot+1;
        Filenum1=Filenum1+1;

    end
end





%% Average for each fibril
for i=1:count-1
    pdist_avg(1,i)=mean(nonzeros(pdist(:,i)));
    pdist_avg_std(1,i)=std(nonzeros(pdist(:,i)));
    hdist_pk_avg(1,i)=mean(nonzeros(hdist_pk(:,i)));
    hdist_pk_avg_std(1,i)=std(nonzeros(hdist_pk(:,i)));
    hdist_trghs_avg(1,i)=mean(nonzeros(hdist_trghs(:,i)));
    hdist_trghs_avg_std(1,i)=std(nonzeros(hdist_trghs(:,i)));
end

idx1=isnan(pdist_avg);pdist_avg(idx1)=0;
idx2=isnan(pdist_avg_std);pdist_avg_std(idx2)=0;

idx3=isnan(hdist_pk_avg);hdist_pk_avg(idx3)=0;
idx4=isnan(hdist_pk_avg_std);hdist_pk_avg_std(idx4)=0;

idx5=isnan(hdist_trghs_avg);hdist_trghs_avg(idx5)=0;
idx6=isnan(hdist_trghs_avg_std);hdist_trghs_avg_std(idx6)=0;

clear x1 x2



%% Creating colour matrix for scatter plot
colourmatrix=rand(datanum,3);
hold on;

%% Scatter plot: peaks
figure('Name','Mean and SD analysis');
subplot(4,1,1);
scatt(n,datanum,Dataname,scatx,count_save,hdist_pk_avg,colourmatrix,'Average per fibril [nm]',Legname,'Peak height')

%% Scatter plot: peak to peak distance
subplot(4,1,2);
scatt(n,datanum,Dataname,scatx,count_save,pdist_avg,colourmatrix,'Average per fibril [\mum]',Legname,'Distance between adjacent peaks')

%% Scatter plot: troughs
subplot(4,1,3);
scatt(n,datanum,Dataname,scatx,count_save,hdist_trghs_avg,colourmatrix,'Average per fibril [nm]',Legname,'Trough height')

%% Scatter plot: fibril length
subplot(4,1,4);
scatt(n,datanum,Dataname,scatx,count_fiblength_save,fiblength,colourmatrix,'Average per fibril [\mum]',Legname,'Fibril length')
hold off;

%% Displaying height profile images
ImageDisp("HeightProf",Dataname,n,Filename,datanum,Legname,"Height Profile")

%% Displaying skeletal images
ImageDisp("Skel",Dataname,n,Filename,datanum,Legname,"Skeleton")


% %% Histograms
% figure;
% %count_scat=1;
% for q=1:length(count_save)-1
%     %Plot all data in a matrix
%     x1=count_save(q)+1;
%     x2=count_save(q+1);
%     subplot(1,datanum,q);
%     histfit(abs(nonzeros(pdist_avg(:,x1:x2))));
%     xlabel('Distance between adjacent peaks [\mum]');
%     ylabel('Frequency');
% 
%     x1=count_save(q+1);
%     x2=count_save(q);
%     zeroval=length(find(hdist_pk_avg(x1:x2)==0));
%     title(sprintf('%s (%d%s)',Legname(q),x1-x2-zeroval,' fibrils'));
% end
% suptitle('Histogram of peak to peak distance');
% 
% 
% clear x1 x2
% figure;
% for q=1:length(count_save)-1
%     %Plot all data in a matrix
%     x1=count_save(q)+1;
%     x2=count_save(q+1);
%     subplot(1,datanum,q);
%     histfit(nonzeros(hdist_pk_avg(:,x1:x2)));
%     xlabel('Height of peaks [nm]');
%     ylabel('Frequency');
%     x1=count_save(q+1);
%     x2=count_save(q);
%     zeroval=length(find(hdist_pk_avg(x1:x2)==0));
%     title(sprintf('%s (%d%s)',Legname(q),x1-x2-zeroval,' fibrils'));
% end
% 
% suptitle('Histogram of peaks');
% clear x1 x2
% figure;
% for q=1:length(count_save)-1
%     %Plot all data in a matrix
%     x1=count_save(q)+1;
%     x2=count_save(q+1);
%     subplot(1,datanum,q);
%     histfit(nonzeros(hdist_trghs_avg(:,x1:x2)));
%     xlabel('Height of troughs [nm]');
%     ylabel('Frequency');
%     x1=count_save(q+1);
%     x2=count_save(q);
%     zeroval=length(find(hdist_pk_avg(x1:x2)==0));
%     title(sprintf('%s (%d%s)',Legname(q),x1-x2-zeroval,' fibrils'));
% end
% 
% suptitle('Histograms of troughs');
% 
% 
% clear x1 x2
% figure;
% %count_scat=1;
% for q=1:length(count_fiblength_save)-1
%     %Plot all data in a matrix
%     x1=count_fiblength_save(q)+1;
%     x2=count_fiblength_save(q+1);
%     subplot(1,datanum,q);
%     histfit(abs(nonzeros(fiblength(:,x1:x2))));
%     xlabel('Length of fibrils [\mum]');
%     ylabel('Frequency');
%     title(sprintf('%s (%d%s)',Legname(q),x2-x1,' fibrils'));
% end
% suptitle('Histogram of fibril length');
