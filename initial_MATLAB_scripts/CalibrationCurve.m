function [MEAN,STD,CON,P]=CalibrationCurve(CON)
fileList=dir('*.gpr');

for i=1:length(fileList)
    gprHolder=gprread(fileList(i).name);
    if gt(i,5.5)==1
        gprHolder.Data(:,4)=gprHolder.Data(:,4).*10;
    end
    MEAN(i)=mean(gprHolder.Data(:,4));
    STD(i)=std(gprHolder.Data(:,4));
end

errorbar(CON,MEAN,STD,'bo')

P=CON'\MEAN';
X=0:20:100;
hold on
plot(X,X.*P,'r')
axis([-.5 100 0 100000])
xlabel('DNA Concentration (nM)')
ylabel('DNA Intensity (a.u.)')
title('Alexa Fluor-647 Calibration Curve (Setup 2, 3x3 Binning, 200ms)')
stringG=sprintf('Y = m * X + b \n m = %f \n b = %f',P,MEAN(1));
annotation('textbox',[.65 .25 0 0],'String',stringG,'FitBoxToText','on')
end
