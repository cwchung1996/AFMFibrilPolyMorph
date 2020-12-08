function scatt(n,datanum,Dataname,scatx,count_fiblength_save,fiblength,colourmatrix,ylab,Legname,titl)
x1=1;
x2=n(1);
plotx1=1:datanum;
a=1:length(Dataname);
b=1:length(Dataname);
ad=0.05;
for q=1:length(Dataname)
    
    for i=1:length(scatx)
        aa=a(q)-ad;
        bb=b(q)+ad;
        r(i)=aa+(bb-aa)*rand(1,1);
    end
    
    %Plot all data in a matrix
    x1=count_fiblength_save(q)+1;
    x2=count_fiblength_save(q+1);
    
    meanpeaks=mean(nonzeros(fiblength(:,x1:x2)));
    stdpeaks=std(nonzeros(fiblength(:,x1:x2)));
    
    p(q)=scatter(r(1:length(nonzeros(fiblength(:,x1:x2))))',nonzeros(fiblength(:,x1:x2)),25,colourmatrix(q,:),'o','filled','MarkerFaceAlpha',.5,'MarkerEdgeAlpha',.5);hold on;
    errorbar(plotx1(q),meanpeaks,stdpeaks,stdpeaks,[],[],'Marker','+','LineWidth',0.5,'color','black');
    
    
end

axis([0,datanum+1, -inf,inf]);
ylabel(ylab);
legend([p(1),p(2)],Legname);
title(titl);
set(gca,'xtick',[])
set(gca,'xticklabel',[]);hold off;
end