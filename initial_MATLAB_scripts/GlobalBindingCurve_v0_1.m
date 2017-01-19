function GlobalBindingCurve_v0_1 (ratio , solubilized)

NumSamples=length(ratio)/16; 
Ydata=reshape(ratio,[16,NumSamples])';

% Initial estimates for [Bmax , Kd1 , Kd2 , ... , KdN]
SampleHolder=ones(1,NumSamples);
p0 = [ 2 , SampleHolder*100 ];

% Estimate parameters 
fn = @(p , solubilized ) objFcn( p , solubilized , NumSamples ) ;% Function handle

lb = [.1 , SampleHolder*0 ]; 
ub = [20 , SampleHolder*10000]; 
options = optimoptions('lsqcurvefit','Display','iter-detailed'); 
options.TolFun = 1e-10;
options.MaxFunEvals = 10000;
[pFit,~,Residual ] = lsqcurvefit( fn , p0 , solubilized , Ydata , lb , ub , options );

SSE=sum(Residual.^2,2);
% assignin('base','pFit',pFit); %save fitted parameters to workspace

mkdir('96 Well Fit')
cd('96 Well Fit')
fid=fopen('fitparameters.csv','w');
fprintf(fid,'%s,%s,%s,%s','Oligo','Saturation Point','Fitted Kd','RSS');

%Generate figures with scattered data and fit
    for i=1:NumSamples
        pFig=figure();
        line(0:5:200,((0:5:200)*pFit(1))./((0:5:200)+pFit(i+1)),'Color','r')
        hold on
        scatter(solubilized((i-1)*16+1:(i-1)*16+16)/(81.4*30/10),ratio((i-1)*16+1:(i-1)*16+16))
        
        xlabel 'DNA Concentration (nM)'
        ylabel 'Binding Ratio'
        TitleName=['Calibration Plate Fit - Well ' num2str(i*2)];
        title(TitleName)
        axis([0 120 0 8])
        grid on
        
        strinG=sprintf('y = (x * Bmax)/ (x + Kd) \n Bmax = %f \n Kd = %f \n SSE = %f',pFit(1),pFit(i+1),SSE(i));
        annotation('textbox',[.65 .3 0 0],'String',strinG,'FitBoxToText','on')
        saveas(gcf,['Fit_CalibrationPlate_' num2str(i*2)],'png');
        close(pFig)
        fprintf(fid,'\n %u,%2f,%2f,%4f',i*2,pFit(1),pFit(i+1),SSE(i));
        
    end
    
    fclose(fid);

end

function yfit = objFcn(p , solubilized , NumSamples )
yfit=zeros(NumSamples,16);
Bmax = p(1) ; % curve independent parameter

    for i=1:NumSamples
        Kd=p(i+1); % curve dependent parameter
        xdata=solubilized((i-1)*16+1:(i-1)*16+16)/(81.4*30/10);
        
        yfit(i,:) = ( xdata * Bmax ) ./ ( xdata + Kd );
    end

end

