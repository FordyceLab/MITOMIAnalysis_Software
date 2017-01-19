function []=CompileMulti96_v0_0(dataName)
%Use file" {Name}_AnalysisData_NoHeaders.txt
data=load([dataName '_NoHeaders.txt']);

%% Figure out which spots came from which well
dataindex=data(:,1)~=0;
data_clean=data(dataindex,:);
NumSamples=length(data_clean);

%Identify wells on 384 plate wrt 96 plate source
PlateFormat=[1 97; 193 289];
PlateFormat2=[ PlateFormat , PlateFormat + 1 , PlateFormat + 2 , PlateFormat + 3 , PlateFormat + 4 , PlateFormat + 5 , PlateFormat + 6 , PlateFormat + 7 , PlateFormat + 8 , PlateFormat + 9 , PlateFormat + 10 , PlateFormat + 11];
Plate = [PlateFormat2; PlateFormat2 + 12; PlateFormat2 + 24; PlateFormat2 + 36; PlateFormat2 + 48; PlateFormat2 + 60; PlateFormat2 + 72; PlateFormat2 + 84];
PlateToPin=reshape(Plate',[4,96])';

for i=1:16
    Plate384(i,:)=[1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24]+(24*(i-1));
end
PlateToPin384=reshape(Plate384',[4,96])';

for i=1:96
PlateToPinPrint(((i-1)*4+1):((i-1)*4+4),:)=[PlateToPin(i,:);PlateToPin(i,:);PlateToPin(i,:);PlateToPin(i,:)];
PlateToPinPrint384(((i-1)*4+1):((i-1)*4+4),1:4)=[PlateToPin384(i,:);PlateToPin384(i,:);PlateToPin384(i,:);PlateToPin384(i,:)];
end

PlateToSpot=reshape(PlateToPinPrint,[1,1536])';
PlateToSpot384=reshape(PlateToPinPrint384,[1,1536])';
data_wellnum=sortrows([PlateToSpot,data_clean,PlateToSpot384],1);

%% Save data
numFrames=1;
HeaderFormat={'Well','Index','ColIndx','RowIndx','Removed','Flagged','ButXCoor','ButYCoor','ButRad','ButAutoF','BNDMedFG','BNDAvgFG','BNDStdFG','BNDSumFG','BNDSatFG','BNDMedBG','BNDAvgBG','BNDStdBG','BNDSumBG','BNDSatBG'};
for z=1:1
    HeaderFormat(19+z)={['CAPMedFG' num2str(z)]};
    HeaderFormat(19+z+numFrames)={['CAPAvgFG' num2str(z)]};
    HeaderFormat(19+z+numFrames*2)={['CAPStdFG' num2str(z)]};
    HeaderFormat(19+z+numFrames*3)={['CAPSumFG' num2str(z)]};
    HeaderFormat(19+z+numFrames*4)={['CAPSatFG' num2str(z)]};
    HeaderFormat(19+z+numFrames*5)={['CAPMedBG' num2str(z)]};
    HeaderFormat(19+z+numFrames*6)={['CAPAvgBG' num2str(z)]};
    HeaderFormat(19+z+numFrames*7)={['CAPStdBG' num2str(z)]};
    HeaderFormat(19+z+numFrames*8)={['CAPSumBG' num2str(z)]};
    HeaderFormat(19+z+numFrames*9)={['CAPSatBG' num2str(z)]};
end

HeaderFormat(end+1:end+4)={'SOLXCoor','SOLYCoor','SOLRad','SOLAutoF'};
HOffNum=24+10*1;

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

HeaderFormat(end+1)={'384'};
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