function []=CompileConcentrationCurve_v1(dataName,numFrames)
%Use file" {Name}_AnalysisData_NoHeaders.txt
data=load([dataName '_NoHeaders.txt']);

%% Figure out which spots came from which well
dataindex=data(:,1)~=0;
data_clean=data(dataindex,:);
NumSamples=length(data_clean);

%Identify wells on 384 plate wrt 96 plate source

%96 well plate IDs
Plate(1,:)=[1,1,2,2,3,3,4,4,5,5,6,6,7,7,8,8,9,9,10,10,11,11,12,12];
Plate(2,:)=Plate(1,:);
Plate(3,:)=Plate(1,:)+12;
Plate(4,:)=Plate(1,:)+12;
Plate(5,:)=Plate(1,:)+24;
Plate(6,:)=Plate(1,:)+24;
Plate(7,:)=Plate(1,:)+36;
Plate(8,:)=Plate(1,:)+36;
Plate(9,:)=Plate(1,:)+48;
Plate(10,:)=Plate(1,:)+48;
Plate(11,:)=Plate(1,:)+60;
Plate(12,:)=Plate(1,:)+60;
Plate(13,:)=Plate(1,:)+72;
Plate(14,:)=Plate(1,:)+72;
Plate(15,:)=Plate(1,:)+84;
Plate(16,:)=Plate(1,:)+84;
Plate=[Plate;Plate];
PlateToPin=reshape(Plate',[4,192])';

%Concentration IDs
SeqH=[1,2/3;(2/3)^2,(2/3)^3];
SeqL=[(2/3)^4,(2/3)^5;(2/3)^6,(2/3)^7];
ConcentrationH=[SeqH,SeqH,SeqH,SeqH,SeqH,SeqH,SeqH,SeqH,SeqH,SeqH,SeqH,SeqH];
ConcentrationL=[SeqL,SeqL,SeqL,SeqL,SeqL,SeqL,SeqL,SeqL,SeqL,SeqL,SeqL,SeqL];
Concentration=[ConcentrationH;ConcentrationH;ConcentrationH;ConcentrationH;ConcentrationH;ConcentrationH;ConcentrationH;ConcentrationH;ConcentrationL;ConcentrationL;ConcentrationL;ConcentrationL;ConcentrationL;ConcentrationL;ConcentrationL;ConcentrationL]*16;
ConcentrationToPin=reshape(Concentration',[4,192])';

%384 well plate IDs
for i=1:32
    Plate384(i,:)=[1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24]+(24*(i-1));
end
PlateToPin384=reshape(Plate384',[4,192])';

%Reformatted for data file
for i=1:192
PlateToPinPrint(((i-1)*2+1):((i-1)*2+2),:)=[PlateToPin(i,:);PlateToPin(i,:)];
ConcentrationToPinPrint(((i-1)*2+1):((i-1)*2+2),:)=[ConcentrationToPin(i,:);ConcentrationToPin(i,:)];
PlateToPinPrint384(((i-1)*2+1):((i-1)*2+2),:)=[PlateToPin384(i,:);PlateToPin384(i,:)];
end

PlateToSpot=reshape(PlateToPinPrint,[1,1536])';
ConcentrationToSpot=reshape(ConcentrationToPinPrint,[1,1536])';
PlateToSpot384=reshape(PlateToPinPrint384,[1,1536])';

data_wellnum=sortrows([PlateToSpot,data_clean,PlateToSpot384,ConcentrationToSpot],1);

%% Save data

HeaderFormat={'Well','Index','ColIndx','RowIndx','Removed','Flagged','ButXCoor','ButYCoor','ButRad','BuAreaFG','BuAreaBG','ButAutoF','BNDMedFG','BNDAvgFG','BNDStdFG','BNDSumFG','BNDSatFG','BNDMedBG','BNDAvgBG','BNDStdBG','BNDSumBG','BNDSatBG'};
for z=1:numFrames
    HeaderFormat(22+z)={['CAPMedFG' num2str(z)]};
    HeaderFormat(22+z+numFrames)={['CAPAvgFG' num2str(z)]};
    HeaderFormat(22+z+numFrames*2)={['CAPStdFG' num2str(z)]};
    HeaderFormat(22+z+numFrames*3)={['CAPSumFG' num2str(z)]};
    HeaderFormat(22+z+numFrames*4)={['CAPSatFG' num2str(z)]};
    HeaderFormat(22+z+numFrames*5)={['CAPMedBG' num2str(z)]};
    HeaderFormat(22+z+numFrames*6)={['CAPAvgBG' num2str(z)]};
    HeaderFormat(22+z+numFrames*7)={['CAPStdBG' num2str(z)]};
    HeaderFormat(22+z+numFrames*8)={['CAPSumBG' num2str(z)]};
    HeaderFormat(22+z+numFrames*9)={['CAPSatBG' num2str(z)]};
end

HeaderFormat(end+1:end+6)={'SOLXCoor','SOLYCoor','SOLRad','SOAreaFG','SOAreaBG','SOLAutoF'};
HOffNum=28+10*numFrames;

for SolFrame=1:1
    HeaderFormat(HOffNum+SolFrame)={['SOLMedFG' num2str(SolFrame)]};
    HeaderFormat(HOffNum+SolFrame+1)={['SOLAvgFG' num2str(SolFrame)]};
    HeaderFormat(HOffNum+SolFrame+1*2)={['SOLStdFG' num2str(SolFrame)]};
    HeaderFormat(HOffNum+SolFrame+1*3)={['SOLSumFG' num2str(SolFrame)]};
    HeaderFormat(HOffNum+SolFrame+1*4)={['SOLSatFG' num2str(SolFrame)]};
    HeaderFormat(HOffNum+SolFrame+1*5)={['SOLMedBG' num2str(SolFrame)]};
    HeaderFormat(HOffNum+SolFrame+1*6)={['SOLAvgBG' num2str(SolFrame)]};
    HeaderFormat(HOffNum+SolFrame+1*7)={['SOLStdBG' num2str(SolFrame)]};
    HeaderFormat(HOffNum+SolFrame+1*8)={['SOLSumBG' num2str(SolFrame)]};
    HeaderFormat(HOffNum+SolFrame+1*9)={['SOLSatBG' num2str(SolFrame)]};
end

HeaderFormat(end+1:end+2)={'384','[uM]'};
fName=input('Name files: ' , 's');
DataText=fopen([fName '_sorted.txt'],'w');
fprintf(DataText,'%s\t',HeaderFormat{:} );
fprintf(DataText,'\r\n');

for Z=1:NumSamples
    fprintf(DataText,'%.2f\t',data_wellnum(Z,:));
    fprintf(DataText,'\r\n');
end
fclose(DataText);
save([fName '_sorted.mat'],'data_wellnum');


end