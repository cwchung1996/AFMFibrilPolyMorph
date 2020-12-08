function ImageDisp(Imgname,Dataname,n,Filename,datanum,Legname,suptitl)
figure('Name',Imgname);

% Creating a dummy variable as figure placeholder
dum1=1:max(n);

if length(Dataname)>1
    for i=2:length(Dataname)
        dum1(i,:)=dum1(i-1,:)+max(n);
    end
else
end


g=1;
while g<sum(n)+1
    for dn=1:length(Dataname)
        name=Filename(dn);
        %% FOR-LOOP 2: Importing series of FLIM data
        for k = 1:n(dn)
            figa=sprintf('%s_%s%d.tif',Imgname,name,k);hold on;
            subplot(datanum,max(n),dum1(dn,k));
            imshow(figa);
            title(sprintf('%s #%d',Legname(dn),k));
            g=g+1;
        end
    end
end
suptitle(suptitl);hold off;