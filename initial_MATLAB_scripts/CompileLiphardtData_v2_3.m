function []=CompileLiphardtData_v2_3(dataName,exposure)
%Use file" {Name}_AnalysisData_NoHeaders.txt
data=load([dataName '_NoHeaders.txt']);
% Log=load([dataName '.mat']);
%% Figure out which spots came from which oligo

dataindex=data(:,1)~=0;
data_clean=data(dataindex,:);

for a=1:4
    for b=1:8
        for c=1:40
            DNAConc((a-1)*320+(b-1)*40+c)=b;          
        end
    end
end

DNAConc(find(DNAConc==1))=16;
DNAConc(find(DNAConc==2))=10.7;
DNAConc(find(DNAConc==3))=7.1;
DNAConc(find(DNAConc==4))=4.7;
DNAConc(find(DNAConc==5))=3.2;
DNAConc(find(DNAConc==6))=2.1;
DNAConc(find(DNAConc==7))=1.4;
DNAConc(find(DNAConc==8))=0.94;

data_clean2=[data_clean,DNAConc'];

holder2=[];
holder=[];
ticker=1;
NumSamples=length(data_clean2);

for a=1:20
   for b=1:8
      d=21-a+(b-1)*20;
      holder=horzcat(holder,[d d d d d d d d]);
   end
end


for aa=1:4
    for a=1:8
        for b=1:5
            if lt((aa)*320,641)
                c=10-(b-1)*2;
            else
                c=11-b*2;
            end
            holder2= horzcat(holder2,[c c c c c c c c]);
        end
    end
end


% save('IndexHolder.mat','dataindex')
data_oligonum=horzcat(holder2',data_clean2);
data_oligonum=sortrows(data_oligonum,1);
data_oligonum=[data_oligonum(:,1),holder',data_oligonum(:,2:end)];

%% Save data
fName=input('Name files: ' , 's');
DataText=fopen([fName '_sortedwDNA.txt'],'w');
HeaderFormat={'OligoNum','Well Num' 'Index','ColIndx','RowIndx','Removed','Flagged','ButXCoor','ButYCoor','ButRad','ButAutoF','BNDMedFG','BNDAvgFG','BNDStdFG','BNDSumFG','BNDSatFG','BNDMedBG','BNDAvgBG','BNDStdBG','BNDSumBG','BNDSatFG','CAPMedFG','CAPAvgFG','CAPStdFG','CAPSumFG','CAPSatFG','CAPMedBG','CAPAvgBG','CAPStdBG','CAPSumBG','CAPSatBG','SOLXCoor','SOLYCoor','SOLRad','SOLAutoF','SOLMedFG','SOLAvgFG','SOLStdFG','SOLSumFG','SOLSatFG','SOLMedBG','SOLAvgBG','SOLStdBG','SOLSumBG','SOLSatBG','DNAConc'};
fprintf(DataText,'%s\t',HeaderFormat{:} );
fprintf(DataText,'\r\n');
    for Z=1:NumSamples
        fprintf(DataText,'%.2f\t',data_oligonum(Z,:));
        fprintf(DataText,'\r\n');
    end
fclose(DataText);
save([fName '_sortedwDNA.mat'],'data_oligonum');
%% Compute data
filter=data_oligonum(:,7)==1;
data_oligonum(filter,:)=0;
% filter=(data_oligonum(:,39)-data_oligonum(:,38))<=0 | (data_oligonum(:,14)-data_oligonum(:,18))>=0 | filter==1;

Content.Surface_MedianFG=reshape(data_oligonum(:,12),[128,10])'.*reshape(~filter,[128,10])';
Content.Surface_MedianBG=reshape(data_oligonum(:,17),[128,10])'.*reshape(~filter,[128,10])';
Content.Surface_BGSubMedian=Content.Surface_MedianFG-Content.Surface_MedianBG;

Content.Surface_TotalFG=reshape(data_oligonum(:,15),[128,10])'.*reshape(~filter,[128,10])';
Content.Surface_TotalBG=reshape(data_oligonum(:,20),[128,10])'.*reshape(~filter,[128,10])';
Content.Surface_BGSubTotal=Content.Surface_TotalFG-Content.Surface_TotalBG;

Content.Captured_MedianFG=reshape(data_oligonum(:,22),[128,10])'.*reshape(~filter,[128,10])';
Content.Captured_MedianBG=reshape(data_oligonum(:,27),[128,10])'.*reshape(~filter,[128,10])';
Content.Captured_BGSubMedian=Content.Captured_MedianFG-Content.Captured_MedianBG;

Content.Captured_TotalFG=reshape(data_oligonum(:,25),[128,10])'.*reshape(~filter,[128,10])';
Content.Captured_TotalBG=reshape(data_oligonum(:,30),[128,10])'.*reshape(~filter,[128,10])';
Content.Captured_BGSubTotal=Content.Captured_TotalFG-Content.Captured_TotalBG;

Content.Solubilized_MedianFG=(reshape(data_oligonum(:,36),[128,10]))'./(81.4*exposure/10).*reshape(~filter,[128,10])';
Content.Solubilized_MedianBG=(reshape(data_oligonum(:,41),[128,10]))'./(81.4*exposure/10).*reshape(~filter,[128,10])';
Content.Solubilized_BGSubMedian=Content.Solubilized_MedianFG-Content.Solubilized_MedianBG;

Content.Solubilized_TotalFG=reshape(data_oligonum(:,39),[128,10])'.*reshape(~filter,[128,10])';
Content.Solubilized_TotalBG=reshape(data_oligonum(:,44),[128,10])'.*reshape(~filter,[128,10])';
Content.Solubilized_BGSubTotal=Content.Solubilized_TotalFG-Content.Solubilized_TotalBG;

Content.Ratio_TotalFG=Content.Captured_TotalFG./Content.Surface_TotalFG;
Content.Ratio_TotalBG=Content.Captured_TotalBG./Content.Surface_TotalBG;
Content.Ratio_BGSubTotal=Content.Captured_BGSubTotal./Content.Surface_BGSubTotal;

Content.Ratio_MedianFG=Content.Captured_MedianFG./Content.Surface_MedianFG;
Content.Ratio_MedianBG=Content.Captured_MedianBG./Content.Surface_MedianBG;
Content.Ratio_BGSubMedian=Content.Captured_BGSubMedian./Content.Surface_BGSubMedian;

Content.OligoNames{10}='10 - Classic MAR';
Content.OligoNames{9}='9 - Random DNA Sequence';
Content.OligoNames{1}='1 - (MAR1a) #vl3-4770';
Content.OligoNames{2}='2 - (MAR1b) #vl3-5940';
Content.OligoNames{3}='3 - (MAR2a) #vl3-3898';
Content.OligoNames{4}='4 - (MAR2b) #vl3-5897';
Content.OligoNames{5}='5 - (MAR3a) #10A-15792';
Content.OligoNames{6}='6 - (MAR3b) #10A-12346';
Content.OligoNames{7}='7 - (MAR4a) #10a-31178';
Content.OligoNames{8}='8 - (MAR4b) #10a-30979';
NOTE='Each row is a unique oligo ranging from highest concentration to lowest';
disp(NOTE)
save([fName '_computed.mat'],'Content')

end