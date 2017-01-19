imageSeq=cell(16,96,21);
for i=4:24
    image=imread(['Stitched_' num2str(i) '_Cy5.tif']);
    dimImage=sort(size(image));
    frame=imresize(image,7500/dimImage(2));
    frame= fliplr(frame);
    frame=flipud(frame);
    clearvars image
    disp(i)
    
    for j=1:96
        imageSeq{1,j,i-3}=frame((data_wellnum(1+(j-1)*16,234)-data_wellnum(1+(j-1)*16,235)):(data_wellnum(1+(j-1)*16,234)+data_wellnum(1+(j-1)*16,235)),(data_wellnum(1+(j-1)*16,233)-data_wellnum(1+(j-1)*16,235)):(data_wellnum(1+(j-1)*16,233)+data_wellnum(1+(j-1)*16,235)));
        imageSeq{2,j,i-3}=frame((data_wellnum(2+(j-1)*16,234)-data_wellnum(2+(j-1)*16,235)):(data_wellnum(2+(j-1)*16,234)+data_wellnum(2+(j-1)*16,235)),(data_wellnum(2+(j-1)*16,233)-data_wellnum(2+(j-1)*16,235)):(data_wellnum(2+(j-1)*16,233)+data_wellnum(2+(j-1)*16,235)));
        imageSeq{3,j,i-3}=frame((data_wellnum(3+(j-1)*16,234)-data_wellnum(3+(j-1)*16,235)):(data_wellnum(3+(j-1)*16,234)+data_wellnum(3+(j-1)*16,235)),(data_wellnum(3+(j-1)*16,233)-data_wellnum(3+(j-1)*16,235)):(data_wellnum(3+(j-1)*16,233)+data_wellnum(3+(j-1)*16,235)));
        imageSeq{4,j,i-3}=frame((data_wellnum(4+(j-1)*16,234)-data_wellnum(4+(j-1)*16,235)):(data_wellnum(4+(j-1)*16,234)+data_wellnum(4+(j-1)*16,235)),(data_wellnum(4+(j-1)*16,233)-data_wellnum(4+(j-1)*16,235)):(data_wellnum(4+(j-1)*16,233)+data_wellnum(4+(j-1)*16,235)));
        imageSeq{5,j,i-3}=frame((data_wellnum(5+(j-1)*16,234)-data_wellnum(5+(j-1)*16,235)):(data_wellnum(5+(j-1)*16,234)+data_wellnum(5+(j-1)*16,235)),(data_wellnum(5+(j-1)*16,233)-data_wellnum(5+(j-1)*16,235)):(data_wellnum(5+(j-1)*16,233)+data_wellnum(5+(j-1)*16,235)));
        imageSeq{6,j,i-3}=frame((data_wellnum(6+(j-1)*16,234)-data_wellnum(6+(j-1)*16,235)):(data_wellnum(6+(j-1)*16,234)+data_wellnum(6+(j-1)*16,235)),(data_wellnum(6+(j-1)*16,233)-data_wellnum(6+(j-1)*16,235)):(data_wellnum(6+(j-1)*16,233)+data_wellnum(6+(j-1)*16,235)));
        imageSeq{7,j,i-3}=frame((data_wellnum(7+(j-1)*16,234)-data_wellnum(7+(j-1)*16,235)):(data_wellnum(7+(j-1)*16,234)+data_wellnum(7+(j-1)*16,235)),(data_wellnum(7+(j-1)*16,233)-data_wellnum(7+(j-1)*16,235)):(data_wellnum(7+(j-1)*16,233)+data_wellnum(7+(j-1)*16,235)));
        imageSeq{8,j,i-3}=frame((data_wellnum(8+(j-1)*16,234)-data_wellnum(8+(j-1)*16,235)):(data_wellnum(8+(j-1)*16,234)+data_wellnum(8+(j-1)*16,235)),(data_wellnum(8+(j-1)*16,233)-data_wellnum(8+(j-1)*16,235)):(data_wellnum(8+(j-1)*16,233)+data_wellnum(8+(j-1)*16,235)));
        imageSeq{9,j,i-3}=frame((data_wellnum(9+(j-1)*16,234)-data_wellnum(9+(j-1)*16,235)):(data_wellnum(9+(j-1)*16,234)+data_wellnum(9+(j-1)*16,235)),(data_wellnum(9+(j-1)*16,233)-data_wellnum(9+(j-1)*16,235)):(data_wellnum(9+(j-1)*16,233)+data_wellnum(9+(j-1)*16,235)));
        imageSeq{10,j,i-3}=frame((data_wellnum(10+(j-1)*16,234)-data_wellnum(10+(j-1)*16,235)):(data_wellnum(10+(j-1)*16,234)+data_wellnum(10+(j-1)*16,235)),(data_wellnum(10+(j-1)*16,233)-data_wellnum(10+(j-1)*16,235)):(data_wellnum(10+(j-1)*16,233)+data_wellnum(10+(j-1)*16,235)));
        imageSeq{11,j,i-3}=frame((data_wellnum(11+(j-1)*16,234)-data_wellnum(11+(j-1)*16,235)):(data_wellnum(11+(j-1)*16,234)+data_wellnum(11+(j-1)*16,235)),(data_wellnum(11+(j-1)*16,233)-data_wellnum(11+(j-1)*16,235)):(data_wellnum(11+(j-1)*16,233)+data_wellnum(11+(j-1)*16,235)));
        imageSeq{12,j,i-3}=frame((data_wellnum(12+(j-1)*16,234)-data_wellnum(12+(j-1)*16,235)):(data_wellnum(12+(j-1)*16,234)+data_wellnum(12+(j-1)*16,235)),(data_wellnum(12+(j-1)*16,233)-data_wellnum(12+(j-1)*16,235)):(data_wellnum(12+(j-1)*16,233)+data_wellnum(12+(j-1)*16,235)));
        imageSeq{13,j,i-3}=frame((data_wellnum(13+(j-1)*16,234)-data_wellnum(13+(j-1)*16,235)):(data_wellnum(13+(j-1)*16,234)+data_wellnum(13+(j-1)*16,235)),(data_wellnum(13+(j-1)*16,233)-data_wellnum(13+(j-1)*16,235)):(data_wellnum(13+(j-1)*16,233)+data_wellnum(13+(j-1)*16,235)));
        imageSeq{14,j,i-3}=frame((data_wellnum(14+(j-1)*16,234)-data_wellnum(14+(j-1)*16,235)):(data_wellnum(14+(j-1)*16,234)+data_wellnum(14+(j-1)*16,235)),(data_wellnum(14+(j-1)*16,233)-data_wellnum(14+(j-1)*16,235)):(data_wellnum(14+(j-1)*16,233)+data_wellnum(14+(j-1)*16,235)));
        imageSeq{15,j,i-3}=frame((data_wellnum(15+(j-1)*16,234)-data_wellnum(15+(j-1)*16,235)):(data_wellnum(15+(j-1)*16,234)+data_wellnum(15+(j-1)*16,235)),(data_wellnum(15+(j-1)*16,233)-data_wellnum(15+(j-1)*16,235)):(data_wellnum(15+(j-1)*16,233)+data_wellnum(15+(j-1)*16,235)));
        imageSeq{16,j,i-3}=frame((data_wellnum(16+(j-1)*16,234)-data_wellnum(16+(j-1)*16,235)):(data_wellnum(16+(j-1)*16,234)+data_wellnum(16+(j-1)*16,235)),(data_wellnum(16+(j-1)*16,233)-data_wellnum(16+(j-1)*16,235)):(data_wellnum(16+(j-1)*16,233)+data_wellnum(16+(j-1)*16,235)));
        
    end
end

save('imageSeq','imageSeq')
mkdir('96 Well Movies')
cd('96 Well Movies')

for i=1:96
    pFig=figure();
    vidObjc=VideoWriter(['MycMITOMI_' num2str(i)],'MPEG-4');
    open(vidObjc)
    for k=1:21
        for j=1:16
            subplot(4,4,j)
            imshow(imageSeq{j,i,k})
            title(['Well' num2str(data_wellnum((i-1)*16+j,1)) ' - Spot ' num2str(data_wellnum((i-1)*16+j,2))])
        end
        M=getframe(pFig);
        writeVideo(vidObjc,M)
    end
    close(vidObjc)
    close(pFig)
end
