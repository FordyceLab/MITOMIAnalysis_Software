function []=MITOMIConcentrationImages()
warning off
imageSeq=cell(16,96,1);
ImName=input('Name of image file: ','s');
FileName=input('Name of data file: ', 's');
load([FileName '.mat'],'data_wellnum');
image=imread([ImName '.tif']);
dimImage=sort(size(image));
frame=imresize(image,7500/dimImage(2));
frame= fliplr(frame);
frame=flipud(frame);
ColorSet=uint16(varycolor(20)*65535);
mkdir('96 Well Concentrations')
cd('96 Well Concentrations')
clearvars image
    
for i=1:96
    
    idx=find(data_wellnum(:,1)==i);
    data=data_wellnum(idx,:);
    
    imageSeq{1,i}=frame((data(16,34)-data(16,35)):(data(16,34)+data(16,35)),(data(16,33)-data(16,35)):(data(16,33)+data(16,35)));
    imageSeq{2,i}=frame((data(15,34)-data(15,35)):(data(15,34)+data(15,35)),(data(15,33)-data(15,35)):(data(15,33)+data(15,35)));
    imageSeq{5,i}=frame((data(14,34)-data(14,35)):(data(14,34)+data(14,35)),(data(14,33)-data(14,35)):(data(14,33)+data(14,35)));
    imageSeq{6,i}=frame((data(13,34)-data(13,35)):(data(13,34)+data(13,35)),(data(13,33)-data(13,35)):(data(13,33)+data(13,35)));
    imageSeq{9,i}=frame((data(12,34)-data(12,35)):(data(12,34)+data(12,35)),(data(12,33)-data(12,35)):(data(12,33)+data(12,35)));
    imageSeq{10,i}=frame((data(11,34)-data(11,35)):(data(11,34)+data(11,35)),(data(11,33)-data(11,35)):(data(11,33)+data(11,35)));
    imageSeq{13,i}=frame((data(10,34)-data(10,35)):(data(10,34)+data(10,35)),(data(10,33)-data(10,35)):(data(10,33)+data(10,35)));
    imageSeq{14,i}=frame((data(9,34)-data(9,35)):(data(9,34)+data(9,35)),(data(9,33)-data(9,35)):(data(9,33)+data(9,35)));
    imageSeq{3,i}=frame((data(8,34)-data(8,35)):(data(8,34)+data(8,35)),(data(8,33)-data(8,35)):(data(8,33)+data(8,35)));
    imageSeq{4,i}=frame((data(7,34)-data(7,35)):(data(7,34)+data(7,35)),(data(7,33)-data(7,35)):(data(7,33)+data(7,35)));
    imageSeq{7,i}=frame((data(6,34)-data(6,35)):(data(6,34)+data(6,35)),(data(6,33)-data(6,35)):(data(6,33)+data(6,35)));
    imageSeq{8,i}=frame((data(5,34)-data(5,35)):(data(5,34)+data(5,35)),(data(5,33)-data(5,35)):(data(5,33)+data(5,35)));
    imageSeq{11,i}=frame((data(4,34)-data(4,35)):(data(4,34)+data(4,35)),(data(4,33)-data(4,35)):(data(4,33)+data(4,35)));
    imageSeq{12,i}=frame((data(3,34)-data(3,35)):(data(3,34)+data(3,35)),(data(3,33)-data(3,35)):(data(3,33)+data(3,35)));
    imageSeq{15,i}=frame((data(2,34)-data(2,35)):(data(2,34)+data(2,35)),(data(2,33)-data(2,35)):(data(2,33)+data(2,35)));
    imageSeq{16,i}=frame((data(1,34)-data(1,35)):(data(1,34)+data(1,35)),(data(1,33)-data(1,35)):(data(1,33)+data(1,35)));
    imageSeq{17,i}=[data(:,50),data(:,39)-data(:,44)];
    
    hFig=figure();
    data2=sortrows(data,-50);
    positions=insertShape(frame,'circle',[data2(1:16,33),data2(1:16,34),data2(1:16,35)],'Color',ColorSet(1:16,:),'LineWidth',10);
    imshow(positions)
    title(['Oligo ' num2str(i) ' Positions'])
    saveas(gcf,['Oligo ' num2str(i) ' Positions'],'png')
    close(hFig);
end

save('ConcentrationData.mat','imageSeq')

for i=1:96
    pFig=figure();
    titleName={'0.936 uM','1.405 uM','2.107 uM','3.160 uM','4.741 uM','7.111 uM','10.667 uM','16 uM'};
    
    for j=1:8
        subplot(4,8,j)
        imshow(imageSeq{(j-1)*2+1,i})
        title(titleName(j))
        subplot(4,8,j+8)
        imshow(imageSeq{(j-1)*2+2,i})
        title(titleName(j))
    end
    
    subplot(4,8,17:32)
    scatter(imageSeq{17,i}(:,1),imageSeq{17,i}(:,2))
    title(['Intensity vs Concentration for Oligo ' num2str(i)])
    ylabel('Net Intensity (a.u.)')
    xlabel('Concentration (uM)')
    saveas(gcf,['WellConcentration_' num2str(i)],'png');
    close(pFig)
end

end