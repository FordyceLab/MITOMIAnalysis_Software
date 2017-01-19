function []=FormatDissociationForR(inputname,outputname)
load([inputname '.mat']);

DataMat=[data_wellnum(:,1),data_wellnum(:,245),data_wellnum(:,2),data_wellnum(:,21:41)-data_wellnum(:,126:146)];

HeaderFormat(1:3)={'96WellNum','384WellNum','Index'};

for i=4:25
HeaderFormat{i}=['NetMedian' num2str(i-2)];
end
DataText=fopen([outputname '.txt'],'w');
fprintf(DataText,'%s\t',HeaderFormat{:} );
fprintf(DataText,'\r\n');
for Z=1:1536
fprintf(DataText,'%.0f\t',DataMat(Z,:)');
fprintf(DataText,'\r\n');
end

fclose(DataText);