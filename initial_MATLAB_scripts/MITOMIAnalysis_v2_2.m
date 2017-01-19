function []=MITOMIAnalysis_v2_2()
% MITOMIAnalysis provides the user with tools for facilitating data
% extraction from the MITOMI microfluidic platform for further analysis.
% MITOMI data analysis has traditionally collected 3 parameters to evaluate
% binding events: localized transcription factor fluorescence, captured DNA
% oligo fluorescence, and initial solubilized DNA fluorescence. 
% 
% To run this function, simply type MITOMIAnalysis() on the command line.
% 
% MITOMIAnalysis has 4 stages:
% 
% 1. User input for loading and setting constraints for source images
% 2. Automated feature finding
% 3. User input for reviewing and correcting the automated process
% 4. Data extraction from features

%% MAIN PROGRAM

[Image,L,ABORT]=ImagePrep();

    if ABORT==1
        cd(L(1).Directory)
        cd ..
        LOGFILENAME=['LOG_MITOMIAnalysis_' datestr(datetime('now')) '.mat'];
        save(LOGFILENAME,'L')
        disp('Log file saved.')
        return
    end
    
[L]=SetCoordinates(Image,L);
[Data,L]=AutomatedFeatureFinding(Image,L);
[Data,L,ABORT]=UserEdit(Image,Data,L);

    if ABORT==1
        cd(L(1).Directory)
        cd ..
        LOGFILENAME=['LOG_MITOMIAnalysis_' datestr(datetime('now')) '.mat'];
        save(LOGFILENAME,'L')
        disp('Log file saved.')
        return
    end
[Data,L]=CompileData(Image,Data,L);
    
fprintMITOMI(L,Data)

%% IMAGEPREP FUNCTION
    function [Image,L,ABORT]=ImagePrep()
        
        %Figure out user's intentions and name dataset
        L.type=questdlg('Select Analysis Type :','MITOMI Analysis','Equilibrium','Dissociation','Custom Multi Input','Equilibrium');
        disp([L.type ' analysis type selected.'])
        
        if strcmp(L.type,'Dissociation')==1
            L.numboundframes=input('Number of time steps starting from t0? : ');
            L.numsolframes=1;
        elseif strcmp(L.type,'Custom Multi Input')==1
            L.numboundframes=input('Number of bound frames? : ');
            L.numsolframes=input('Number of solubilized frames? : ');
        else
            L.numboundframes=1;
            L.numsolframes=1;
        end
        
        L.name=input('Name dataset: ','s');
        disp(['Dataset will be labeled as: ' L.name])
        
        %Select directory containing stitched data images, build file list
        disp('Select directory containing stitched images.')
        cd(uigetdir)
        d=dir('*.tif');
        fileStr={d.name};
        
        %Select images for analysis
        for solframe=1:L.numsolframes
            Sol=listdlg('PromptString',['Select prewash/solubilized molecule image ' num2str(solframe) ' : '],'SelectionMode','single','listsize',[320 300],'ListString',fileStr);
            L.solubilized{solframe}=d(Sol).name;
            solubilized(:,:,solframe)=imread(L.solubilized{solframe});
            disp(['Prewash/solublized molecule image ' num2str(solframe) ' of ' num2str(L.numsolframes) ' loaded: ' L.solubilized{solframe}])
        end

        SB=listdlg('PromptString','Select surface immobilized molecule image: ','SelectionMode','single','listsize',[320 300],'ListString',fileStr);
        L.surface=d(SB).name;
        surface=imread(L.surface);
        disp(['Immobilized molecule image loaded: ' L.surface])
        
        for a=1:L.numboundframes
            C=listdlg('PromptString',['Select captured molecule image ' num2str(a) ' : '],'SelectionMode','single','listsize',[320 300],'ListString',fileStr);
            L.captured{a}=d(C).name;
%             if a==1
%                 captured=imread(L.captured(a));
%             end
            captured(:,:,a)=imread(L.captured{a});
            disp(['Captured molecule image ' num2str(a) ' of ' num2str(L.numboundframes) ' loaded: ' L.captured{a}])
        end

        
        %Scale image size to prevent Matlab from crashing when disp images
        %Scaling small images also aids the auto-detection of features
        dimImage=sort(size(surface));
        Image.surface=imresize(surface,7500/dimImage(2));
        Image.captured=imresize(captured,7500/dimImage(2));
        Image.solubilized=imresize(solubilized,7500/dimImage(2));
        disp('Images rescaled.')

        %Generate preview image and reorient image set
        ImMed=double(median(Image.captured(:)));
        ImPreview=(1-mat2gray(Image.captured(:,:,1), [0 4*ImMed]));
        figOrientation=figure('menubar','none','toolbar','none','numbertitle','off','Name','figOrientation');
        warning('off')
        imshow(ImPreview)
        
        RESPONSE=0;
        
        while RESPONSE~=5
            RESPONSE=menu('Adjust Image Orientation: ','Rotate Clockwise','Rotate Counter Clockwise','Flip Left-to-Right','Flip Top-To-Bottom','>> ACCEPT <<','>> ABORT <<');
            switch RESPONSE
                case 1
                    ImPreview=imrotate(ImPreview,-90);
                    Image= structfun(@(x) (imrotate(x,-90)), Image, 'UniformOutput', false);
                    disp('Image rotated clockwise.')
                case 2
                    ImPreview=imrotate(ImPreview,90);
                    Image= structfun(@(x) (imrotate(x,90)), Image, 'UniformOutput', false);
                    disp('Image rotated counter-clockwise.')
                case 3
                    ImPreview=fliplr(ImPreview);
                    Image= structfun(@fliplr, Image, 'UniformOutput', false);
                    disp('Image flipped left-to-right.')
                case 4
                    ImPreview=flipud(ImPreview);
                    Image= structfun(@flipud, Image, 'UniformOutput', false);
                    disp('Image flipped top-to-bottom.')
                case 5
                    disp('Image orientation accepted.')
                otherwise
                    ABORT=1;
                    L.Error('User aborted program during image ortientation process.');
                    disp(L.Error)
                    return
            end
            imshow(ImPreview)
        end
        
        close(figOrientation)
        warning('on')
        ABORT=0;            

    end

%% SET COORDINATES FUNCTION
    function [L]=SetCoordinates(Image,L)
        %Define grids for automated analysis
        
        %Preset variables
        L.CoorButtons=[];
        L.CoorChamber=[];
        RadiusSample=zeros(1,4);
        SolubilizedRadiusSample=zeros(1,4);
        IntensitySample=zeros(1,4);
        SolubilizedIntSample=zeros(1,4);
        CornerCoor=zeros(4,2);
        SolubilizedCornerCoor=zeros(4,2);
        ImGhost=((1-mat2gray(Image.surface))*1/2 + mat2gray(Image.solubilized(:,:,1))*1/2);
        ImGhostCh=((1-mat2gray(Image.surface))/4 + mat2gray(Image.solubilized(:,:,1))*3/4);
        warning('OFF')
        figButtons=figure('menubar','none','numbertitle','off','toolbar','figure','Name','Array Coordinate Figure');
        figHandles=vertcat(findall(figButtons,'type','uipushtool'),findall(figButtons,'type','uitoggletool'));
        indexHandles=logical(strcmp(get(figHandles,'Tag'),'Exploration.ZoomOut')+strcmp(get(figHandles,'Tag'),'Exploration.ZoomIn')+strcmp(get(figHandles,'Tag'),'Exploration.Pan'));
        notZoomHandles=figHandles(~indexHandles);
        set(notZoomHandles,'enable','off');
        imshow(ImGhost)
        disp('Please enter the following information: ')
        dimensions=inputdlg({'Number of Rows (often 56) :','Number of Columns (often 28) :'},'Array Dimensions',1,{'56','28'});
        L.NumRow=str2double(dimensions{1});
        L.NumCol=str2double(dimensions{2});
        shg
        
        %Define coordinates for button analysis
        for i=1:4
            imshow(ImGhost)
            disp(['GRID CORNER SAMPLE #' num2str(i) ' of 4 - Zoom over a corner then press Enter'])
            title(['GRID CORNER SAMPLE #' num2str(i) ' of 4 - Zoom over a corner then press Enter'])
            zoom on
            pause
            zoom off
            disp('Zoom locked. Click three points on the circumference of the outermost data point.')

            title(['GRID CORNER SAMPLE #' num2str(i) ' of 4 - Click three points on the circumference'])
            [Q,C]=ginput(3);
            p1=[Q(1) C(1) 1]; p2=[Q(2) C(2) 1]; p3=[Q(3) C(3) 1];
            t = p2-p1; u = p3-p1; v = p3-p2;
            w = cross(t,u);
            t2 = sum(t.^2); u2 = sum(u.^2); w2 = sum(w.^2);
            c = p1+(t2*sum(u.*v)*u-u2*sum(t.*v)*t)/(2*w2);
            RadiusSample(i) = (sqrt(t2*u2*sum(v.^2)/w2))/2;
            IntensitySample(i)=double(Image.surface(round(c(1,2)),round(c(1,1))));
            CornerCoor(i,:)=c(1,1:2);
            
            %Define coordinates from solubilized image
            LIMITS=get(gca,{'XLim' 'YLim'});
            imshow(ImGhostCh)
            set(gca,{'XLim' 'YLim'},LIMITS);
            disp('Zoom locked. Click three points on the circumference of the outermost data chamber.')
            title(['GRID CORNER SAMPLE #' num2str(i) ' of 4 - Click three points on the circumference'])
            [S,W]=ginput(3);
            p4=[S(1) W(1) 1]; p5=[S(2) W(2) 1]; p6=[S(3) W(3) 1];
            t = p5-p4; u = p6-p4; v = p6-p5;
            w = cross(t,u);
            t2 = sum(t.^2); u2 = sum(u.^2); w2 = sum(w.^2);
            d = p4+(t2*sum(u.*v)*u-u2*sum(t.*v)*t)/(2*w2);
            SolubilizedRadiusSample(i) = (sqrt(t2*u2*sum(v.^2)/w2))/2;
            SolubilizedCornerCoor(i,:)=d(1,1:2);
            SolubilizedIntSample(i)=double(Image.solubilized(round(d(1,2)),round(d(1,1)),1));
            disp(['Data sampling complete for sample ' num2str(i) '.'])
        end
        
        close('Array Coordinate Figure')
        
        %Reorganize coordinate data collected
        L.approxIntensity=mean(IntensitySample);
        L.Radius=ceil(mean(RadiusSample));        
        vertices=round(sortrows(CornerCoor,2));
        
        verticesSolubilized=round(sortrows(SolubilizedCornerCoor,2));
        L.RadiusSolubilized=ceil(mean(SolubilizedRadiusSample));
        L.approxIntensitySolubilized=max(SolubilizedIntSample);
        
        TopRowXButtons=sort(linspace(vertices(1,1),vertices(2,1),L.NumCol));
        BotRowXButtons=sort(linspace(vertices(3,1),vertices(4,1),L.NumCol));
        TopRowYButtons=interp1(vertices(1:2,1),vertices(1:2,2),TopRowXButtons);
        BotRowYButtons=interp1(vertices(3:4,1),vertices(3:4,2),BotRowXButtons);
        
        TopRowXSolubilized=sort(linspace(verticesSolubilized(1,1),verticesSolubilized(2,1),L.NumCol));
        BotRowXSolubilized=sort(linspace(verticesSolubilized(3,1),verticesSolubilized(4,1),L.NumCol));
        TopRowYSolubilized=interp1(verticesSolubilized(1:2,1),verticesSolubilized(1:2,2),TopRowXSolubilized);
        BotRowYSolubilized=interp1(verticesSolubilized(3:4,1),verticesSolubilized(3:4,2),BotRowXSolubilized);

        %Generate grid from coordinate data
        for j=1:L.NumCol
            ColYValButtons=sort(linspace(TopRowYButtons(j),BotRowYButtons(j),L.NumRow));
            ColXValButtons=interp1([TopRowYButtons(j) BotRowYButtons(j)],[TopRowXButtons(j) BotRowXButtons(j)],ColYValButtons);
            L.CoorButtons=round(vertcat(L.CoorButtons,[ColXValButtons' ColYValButtons']));
            
            ColYValSolubilized=sort(linspace(TopRowYSolubilized(j),BotRowYSolubilized(j),L.NumRow));
            ColXValSolubilized=interp1([TopRowYSolubilized(j) BotRowYSolubilized(j)],[TopRowXSolubilized(j) BotRowXSolubilized(j)],ColYValSolubilized);
            L.CoorChamber=round(vertcat(L.CoorChamber,[ColXValSolubilized' ColYValSolubilized']));
        end
        warning('ON')
        disp('Approximate coordinates for automatation set.')
    end

%% AUTOMATED FEATURE FINDING FUNCTION
    function [Data,L]=AutomatedFeatureFinding(Image,L)
       
        %Adjust radii for mask-making
        L.modRadiusButton=L.Radius+round(L.Radius/2);
        
        buttonFGmask=zeros(L.modRadiusButton*2);
        buttonBGmask=zeros(L.modRadiusButton*2);
        for k=1:L.modRadiusButton*2
            for l=1:L.modRadiusButton*2
                if lt((k-(L.modRadiusButton*2+1)/2)^2+(l-(L.modRadiusButton*2+1)/2)^2,(L.modRadiusButton/2)^2)
                    buttonFGmask(k,l)=1;
                else
                    buttonFGmask(k,l)=0;
                end
                if lt((k-(L.modRadiusButton*2+1)/2)^2+(l-(L.modRadiusButton*2+1)/2)^2,(L.modRadiusButton*2/2)^2) && gt((k-(L.modRadiusButton*2+1)/2)^2+(l-(L.modRadiusButton*2+1)/2)^2,(L.modRadiusButton*1.5/2)^2)
                    buttonBGmask(k,l)=1;
                else
                    buttonBGmask(k,l)=0;
                end
            end
        end
        
        solubilizedFGmask=zeros(L.RadiusSolubilized*2);
        for p=1:L.RadiusSolubilized*2
            for q=1:L.RadiusSolubilized*2
                if lt((p-(L.RadiusSolubilized*2+1)/2)^2+(q-(L.RadiusSolubilized*2+1)/2)^2,L.RadiusSolubilized^2)
                    solubilizedFGmask(p,q)=1;
                end
            end
        end
        
        L.buttonFGmask=buttonFGmask;
        L.buttonBGmask=buttonBGmask;
        L.solubilizedFGmask=solubilizedFGmask;
        
        %Preset variables
        L.BTicker=0;
        L.CTicker=0;
        dimensions=zeros(L.NumCol*L.NumRow,1);
        
        Data.Index=dimensions;
        Data.ColIndex=dimensions;
        Data.RowIndex=dimensions;
        Data.Remove=dimensions;
        Data.Flag=dimensions;
        
        Data.ButtonsXCoor=dimensions;
        Data.ButtonsYCoor=dimensions;
        Data.ButtonsRadius=dimensions;
        Data.ButtonsAreaFG=dimensions;
        Data.ButtonsAreaBG=dimensions;
        Data.AutofindButtons=dimensions;
        
        Data.surfaceMedianFG=dimensions;
        Data.surfaceAverageFG=dimensions;
        Data.surfaceSTDFG=dimensions;
        Data.surfaceTotalFG=dimensions;
        Data.surfaceFractionSaturatedFG=dimensions;
        Data.surfaceMedianBG=dimensions;
        Data.surfaceAverageBG=dimensions;
        Data.surfaceSTDBG=dimensions;
        Data.surfaceTotalBG=dimensions;
        Data.surfaceFractionSaturatedBG=dimensions;
        
        Data.capturedMedianFG=dimensions;
        Data.capturedAverageFG=dimensions;
        Data.capturedSTDFG=dimensions;
        Data.capturedTotalFG=dimensions;
        Data.capturedFractionSaturatedFG=dimensions;
        Data.capturedMedianBG=dimensions;
        Data.capturedAverageBG=dimensions;
        Data.capturedSTDBG=dimensions;
        Data.capturedTotalBG=dimensions;
        Data.capturedFractionSaturatedBG=dimensions;
        
        Data.ChamberXCoor=dimensions;
        Data.ChamberYCoor=dimensions;
        Data.ChamberRadius=dimensions;
        Data.ChamberAreaFG=dimensions;
        Data.ChamberAreaBG=dimensions;
        Data.AutofindChamber=dimensions;
        
        Data.solubilizedMedianFG=dimensions;
        Data.solubilizedMeanFG=dimensions;
        Data.solubilizedSTDFG=dimensions;
        Data.solubilizedTotalFG=dimensions;
        Data.solubilizedFractionSaturatedFG=dimensions;
        Data.solubilizedMedianBG=dimensions;
        Data.solubilizedMeanBG=dimensions;
        Data.solubilizedSTDBG=dimensions;
        Data.solubilizedTotalBG=dimensions;
        Data.solubilizedFractionSaturatedBG=dimensions;

        disp('Variables defined. Beginning automated identification of buttons...')
        
        figButtonGrid=figure('menubar','none','numbertitle','off','toolbar','none','Name','Button Preview');
        WAIT=waitbar(0,'Processing button positions...','Name','Button Positions');
        
        for m=1:length(dimensions);
            
            %Fill in data identity
            Data.ColIndex(m)=ceil(m/L.NumRow);
            Data.RowIndex(m)=m-L.NumRow*(Data.ColIndex(m)-1);
            Data.ButtonsRadius(m,1)=L.Radius;

            CoorX=L.CoorButtons(m,1);
            CoorY=L.CoorButtons(m,2);
            
            %Adjust image intensities and find surface bound spot within 1
            %radius of XY coordinates with slight deviation from defined R
            
            screenSurface=double(Image.surface((CoorY-2*L.modRadiusButton):(CoorY+2*L.modRadiusButton),(CoorX-2*L.modRadiusButton):(CoorX+2*L.modRadiusButton)));
            screenSTD=std(screenSurface(:));
            screenMED=median(screenSurface(:));
            screenSurfaceMod=imadjust(uint16(mat2gray(screenSurface,[screenMED-screenSTD*2 L.approxIntensity+screenSTD*2])*65535));
            [spotLocations,radii]=imfindcircles((screenSurfaceMod),[round(L.modRadiusButton/2.5) round(L.modRadiusButton/1.25)],'ObjectPolarity','bright');
            imshow(screenSurfaceMod)
            
            %if autofind with hough transform finds something, process info
            if isempty(radii)~=1
                %Convert local coordinates to global coordinates
                Data.ButtonsXCoor(m,1)=round(spotLocations(1,1)-L.modRadiusButton*2-1+CoorX);
                Data.ButtonsYCoor(m,1)=round(spotLocations(1,2)-L.modRadiusButton*2-1+CoorY);
                Data.Remove(m)=false;
                Data.AutofindButtons(m)=true;
                
                L.BTicker=L.BTicker+1;
                
            else %Autofind failed, try to find where button is
                surfaceFGdataholder=cell(length(L.modRadiusButton*4+1));
                surfaceBGdataholder=cell(length(L.modRadiusButton*4+1));
                Data.ButtonsRadius(m,1)=L.Radius;

                %Sample data in the phase space
                for n=(CoorX-L.modRadiusButton*2):(CoorX+L.modRadiusButton*2)
                    for o=(CoorY-L.modRadiusButton*2):(CoorY+L.modRadiusButton*2)
                        ExtractImageSurface=double(Image.surface((o-L.modRadiusButton):(o+L.modRadiusButton-1),(n-L.modRadiusButton):(n+L.modRadiusButton-1)));
                        surfaceFGsampletemp=ExtractImageSurface.*L.buttonFGmask;
                        surfaceBGsampletemp=ExtractImageSurface.*L.buttonBGmask;
                        surfaceFGdataholder{n+L.modRadiusButton*2-CoorX+1,o+L.modRadiusButton*2-CoorY+1}=surfaceFGsampletemp(surfaceFGsampletemp>0);
                        surfaceBGdataholder{n+L.modRadiusButton*2-CoorX+1,o+L.modRadiusButton*2-CoorY+1}=surfaceBGsampletemp(surfaceBGsampletemp>0);
                    end
                end
                
                %Find data point w highest net int in phase space datasets
                NetInt=cellfun(@(x,y) (sum(x)-sum(y)),surfaceFGdataholder,surfaceBGdataholder);
                [n_local,o_local]=find(NetInt==max(NetInt(:)));

                %Convert "best data" local coor into global for capt image
                Data.ButtonsXCoor(m)=n_local-1+CoorX-L.modRadiusButton*2;
                Data.ButtonsYCoor(m)=o_local-1+CoorY-L.modRadiusButton*2;   
                Data.Remove(m)=false;
                Data.AutofindButtons(m)=false;
                
            end
            waitbar(m/length(dimensions));
        end
        close(figButtonGrid)
        delete(WAIT)
        disp(['Buttons identified with automation: ' num2str(L.BTicker) ' out of ' num2str(length(dimensions))])            
        disp('Beginning automated identification of chambers...')
            

        figChamberGrid=figure('menubar','none','numbertitle','off','toolbar','none','Name','Chamber Preview');
        WAIT=waitbar(0,'Processing chamber positions...','Name','Chamber Positions');

        for r=1:length(dimensions)

            %Refresh coordinates with chamber positions
            CoorX=L.CoorChamber(r,1);
            CoorY=L.CoorChamber(r,2);
            Data.ChamberRadius(r,1)=L.RadiusSolubilized;

            %Adjust image intensities and find where chamber might be
            screenSol=double(Image.solubilized((CoorY-L.RadiusSolubilized):(CoorY+L.RadiusSolubilized),(CoorX-L.RadiusSolubilized):(CoorX+L.RadiusSolubilized),1));
            screenSolSTD=std(screenSol(:));
            screenSolMED=median(screenSol(:));
            screenSolMod=imadjust(uint16(mat2gray(screenSol,[screenSolMED-screenSolSTD*2 L.approxIntensitySolubilized+screenSolSTD*2])*65535));
            [chamberLocations,chamberradii]=imfindcircles((screenSolMod),[round(L.RadiusSolubilized*.80) round(L.RadiusSolubilized*1.2)],'ObjectPolarity','bright');
            imshow(screenSolMod)

            %if autofind with hough transform finds something, process 
            if isempty(chamberradii)~=1
                %Convert local coordinates to global coordinates
                Data.ChamberXCoor(r)=round(chamberLocations(1,1)-L.RadiusSolubilized-1+CoorX);
                Data.ChamberYCoor(r)=round(chamberLocations(1,2)-L.RadiusSolubilized-1+CoorY);
                Data.AutofindChamber(r)=true;

                L.CTicker=L.CTicker+1;

            else %autofind failed: try, to, find... chamber. so slow.

                solubilizedFGdataholder=cell(length(L.RadiusSolubilized*2+1));

                %sample data for the phase space about the coor
                for s=(CoorX-ceil(L.RadiusSolubilized*7/8)):(CoorX+floor((L.RadiusSolubilized-1)*7/8))
                    for t=(CoorY-ceil(L.RadiusSolubilized*7/8)):(CoorY+floor((L.RadiusSolubilized-1)*7/8))
                        ExtractImageSolubilized=double(Image.solubilized((t-L.RadiusSolubilized):(t+L.RadiusSolubilized-1),(s-L.RadiusSolubilized):(s+L.RadiusSolubilized-1),1));
                        solubilizedFGsampletemp=ExtractImageSolubilized.*L.solubilizedFGmask;
                        solubilizedFGdataholder{s+L.RadiusSolubilized-CoorX+1,t+L.RadiusSolubilized-CoorY+1}=solubilizedFGsampletemp(solubilizedFGsampletemp>0);                             
                    end
                end

                %Find data point w highest net int in phase space datasets
                TotInt=cellfun(@sum,solubilizedFGdataholder);
                [s_local,t_local]=find(TotInt==max(TotInt(:)));

                %Convert "best data" local coor into global for capt image
                Data.ChamberXCoor(r)=s_local(1)-1+CoorX-L.RadiusSolubilized;
                Data.ChamberYCoor(r)=t_local(1)-1+CoorY-L.RadiusSolubilized;
                Data.AutofindChamber(r)=false;
            end

            waitbar(r/length(dimensions));

        end
        close(figChamberGrid)
        delete(WAIT)
        disp(['Chambers identified with automation: ' num2str(L.CTicker) ' out of ' num2str(length(dimensions))])
save('test.mat')
    end
%% USER EDIT FUNCTION
    function [Data,L,ABORT]=UserEdit(Image,Data,L)
        
        %Initialize variables
        ImPreviewButtons=((1-mat2gray(imadjust(Image.surface)))*3/4+mat2gray(Image.solubilized(:,:,1))/4);
        surmenu=0;
        
        while surmenu~=1
            
            %Generate image detailing button locations and detection method
            AutoButtons=find(Data.AutofindButtons==true);
            AutoLength=length(AutoButtons);
            MissButtons=find(Data.AutofindButtons==false);
            MissLength=length(MissButtons);
            CellAuto=cell(AutoLength,1);
            CellFull=cell(AutoLength+MissLength,1);
            iax=cellfun('isempty',CellAuto);
            CellFull(iax)={'green'};
            imx=cellfun('isempty',CellFull);
            CellFull(imx)={'blue'};
            ImWithAutoMiss=insertShape(ImPreviewButtons,'circle',[[Data.ButtonsXCoor(AutoButtons),Data.ButtonsYCoor(AutoButtons),Data.ButtonsRadius(AutoButtons)];[Data.ButtonsXCoor(MissButtons),Data.ButtonsYCoor(MissButtons),Data.ButtonsRadius(MissButtons)]],'Color',CellFull,'LineWidth',3);
            
            %Generate or update graphical feature interface
            warning('OFF')
            if surmenu==0
                apiControl=scrollImage(ImWithAutoMiss);
            else
                apiControl.replaceImage(ImWithAutoMiss,'PreserveView',1);
            end
            warning('ON')
            
            surmenu=menu('Select command : ','Continue (without edits)','Edit Position','ABORT');
        
            switch surmenu
                
                case 1 %Continue to next step
                    
                    disp('Continuing to next stage.')
                    
                case 2 %User chose to manually move button location   

                    %Find object to relocate
                    disp('Click near circle you would like to relocate')
                    h=impoint(gca,[]); 
                    h.setColor('m')
                    [initialXYCoordinates]=h.getPosition();
                    [~,Nearest]=sortrows((Data.ButtonsXCoor-initialXYCoordinates(1)).^2+(Data.ButtonsYCoor-initialXYCoordinates(2)).^2);
                    N=Nearest(1);
                    
                    %Identify where object should be relocated to
                    disp('Click where you would like to center object')
                    g=impoint(gca,[]); 
                    g.setColor('r')
                    [finalXYCoordinates]=round(g.getPosition());
                    
                    Data.ButtonsXCoor(N)=finalXYCoordinates(1);
                    Data.ButtonsYCoor(N)=finalXYCoordinates(2);
                    Data.AutofindButtons(N)=false;
                    
                otherwise %window was closed or user manually aborted
                    ABORT=1;
                    L.Error('User aborted program during User Reposition mode.');
                    disp(L.Error)
                    return
            end
        end
        
        %Remove memory hogs
        close('Graphical Feature Interface')
        clear ImWithAutoMiss ImPreviewButtons
        
        
        %Initialize variables
        bndmenu=0;
        ImPreviewBound=1-mat2gray(imadjust(Image.captured(:,:,1)));

        while bndmenu~=1
            
            %Generate image with active feature locations
            RemainingFeatures=find(Data.Remove==false & Data.Flag==false);
            FlagFeatures=find(Data.Remove==false & Data.Flag==true);
            RemLength=length(RemainingFeatures);
            FlagLength=length(FlagFeatures);
            
            CellRem=cell(RemLength,1);
            CellRF=cell(RemLength+FlagLength,1);
            iax=cellfun('isempty',CellRem);
            CellRF(iax)={'red'};
            irfx=cellfun('isempty',CellRF);
            CellRF(irfx)={'magenta'};

            ImWithFeatures=insertShape(ImPreviewBound,'circle',[[Data.ButtonsXCoor(RemainingFeatures),Data.ButtonsYCoor(RemainingFeatures),Data.ButtonsRadius(RemainingFeatures)];[Data.ButtonsXCoor(FlagFeatures),Data.ButtonsYCoor(FlagFeatures),Data.ButtonsRadius(FlagFeatures)]],'Color',CellRF,'LineWidth',3);
           
            %Generate or update graphical feature interface
            warning('OFF')
            if bndmenu==0
                apiControl=scrollImage(ImWithFeatures);
            else
                apiControl.replaceImage(ImWithFeatures,'PreserveView',1);
            end
            warning('ON')
            
            %User action pane
            bndmenu=menu('Select command : ','Continue (without edits)','Flag','UNDO last flagging','Remove points','UNDO last removal','ABORT');
        
            switch bndmenu
                
                case 1 %User chose to continue
                    disp('Continuing to next stage.')
                    
                case 2
                    disp('Click near a corder of the data you would like to FLAG and drag cursor to opposing corner.')
                    [p1a,p1b]=ginput(1);
                    rbbox; 
                    axes_handle = gca;
                    p2=get(axes_handle,'CurrentPoint');
                    
                    %Sort vertices in preparation of data selection
                    if ( p1a < p2(1,1) )
                       lowX = p1a; highX = p2(1,1);
                    else
                       lowX = p2(1,1); highX = p1a;
                    end

                    if ( p1b < p2(1,2) )
                       lowY = p1b; highY = p2(1,2);
                    else
                       lowY = p2(1,2); highY = p1b;
                    end
                    
                    FlaggedFeatures=find(Data.ButtonsXCoor > lowX & Data.ButtonsXCoor < highX & Data.ButtonsYCoor > lowY & Data.ButtonsYCoor < highY);
                    Data.Flag(FlaggedFeatures)=true;
                    disp('Data points flagged by user.')
                    
                case 3
                    
                    Data.Flag(FlaggedFeatures)=false;
                    disp('Data points unflagged by user.')
                    
                case 4 %User chose to remove points
                    %Identify data points to remove
                    disp('Click near a corner of the data you would like to OMIT and drag cursor to opposing corner.')
                    [p1a,p1b]=ginput(1);
                    rbbox; 
                    axes_handle = gca;
                    p2=get(axes_handle,'CurrentPoint');
                    
                    %Sort vertices in preparation of data selection
                    if ( p1a < p2(1,1) )
                       lowX = p1a; highX = p2(1,1);
                    else
                       lowX = p2(1,1); highX = p1a;
                    end

                    if ( p1b < p2(1,2) )
                       lowY = p1b; highY = p2(1,2);
                    else
                       lowY = p2(1,2); highY = p1b;
                    end
                    
                    RemovedFeatures=find(Data.ButtonsXCoor > lowX & Data.ButtonsXCoor < highX & Data.ButtonsYCoor > lowY & Data.ButtonsYCoor < highY);
                    Data.Remove(RemovedFeatures)=true;
                    disp('Data points removed by user.')
                    
                case 5 %Undo previous removal
                    
                    Data.Remove(RemovedFeatures)=false;
                    disp('Last removal undone by user.')
                    
                otherwise %Window closed or manually aborted
                    ABORT=1;
                    L.Error('User aborted program during User Removal mode.');
                    disp(L.Error)
                    return
            end
        end
        
        %Remove memory hogs
        close('Graphical Feature Interface')
        clear ImWithFeatures ImPreviewBound
        
        solmenu=0;
        ImPreviewChamber=1-mat2gray(imadjust(Image.solubilized(:,:,1)));
        
        while solmenu~=1
            
            %Generate image detailing button locations and detection method
            AutoChamber=find(Data.AutofindChamber==true & Data.Remove==false & Data.Flag==false);
            MissChamber=find(Data.AutofindChamber==false & Data.Remove==false & Data.Flag==false);
            FlagChamber=find(Data.Flag==true & Data.Remove==false);
            
            AuLength=length(AutoChamber);
            MiLength=length(MissChamber);
            FlLength=length(FlagChamber);
            
            CellAu=cell(AuLength,1);
            CellAM=cell(AuLength+MiLength,1);
            CellAMF=cell(AuLength+MiLength+FlLength,1);
            iaux=cellfun('isempty',CellAu);
            CellAM(iaux)={'green'};
            iamx=cellfun('isempty',CellAM);
            CellAMF(iaux)={'green'};
            CellAMF(iamx)={'blue'};
            iamfx=cellfun('isempty',CellAMF);
            CellAMF(iamfx)={'magenta'};
            
            ImWithAMChambers=insertShape(ImPreviewChamber,'circle',[[Data.ChamberXCoor(AutoChamber),Data.ChamberYCoor(AutoChamber),Data.ChamberRadius(AutoChamber)];[Data.ChamberXCoor(MissChamber),Data.ChamberYCoor(MissChamber),Data.ChamberRadius(MissChamber)];[Data.ChamberXCoor(FlagChamber),Data.ChamberYCoor(FlagChamber),Data.ChamberRadius(FlagChamber)]],'Color',CellAMF,'LineWidth',3);
            
            %Generate navigable image and suppress image size warnings
            warning('OFF')
            if solmenu==0
                apiControl=scrollImage(ImWithAMChambers);
            else
                apiControl.replaceImage(ImWithAMChambers,'PreserveView',1);
            end
            warning('ON')
            
            solmenu=menu('Select command : ','Continue (without edits)','Edit Position','ABORT');
        
            switch solmenu
                
                case 1
                    disp('Continuing to next stage.')
                    
                case 2
                    
                    %Find object to relocate
                    disp('Click near chamber you would like to relocate')
                    u=impoint(gca,[]); 
                    u.setColor('m')
                    [initialXYCCoordinates]=u.getPosition();
                    [~,NearestC]=sortrows((Data.ChamberXCoor-initialXYCCoordinates(1)).^2+(Data.ChamberYCoor-initialXYCCoordinates(2)).^2);
                    NC=NearestC(1);
                    
                    %Identify where object should be relocated to
                    disp('Click where you would like to recenter object')
                    v=impoint(gca,[]); 
                    v.setColor('r')
                    [finalXYCCoordinates]=round(v.getPosition());
                    Data.ChamberXCoor(NC)=finalXYCCoordinates(1);
                    Data.ChamberYCoor(NC)=finalXYCCoordinates(2);
                    Data.AutofindChamber(NC)=false;
            
                otherwise %Window closed or manually aborted
                    ABORT=1;
                    L.Error('User aborted program during User Removal mode.');
                    disp(L.Error)
                    return
            end
        end
    close('Graphical Feature Interface')    
    ABORT=0;
    disp('User input stage complete.') 
%     save('tempdata2.mat','Data','Image','L');
    end
%% DATA COMPILIATION FUNCTION

    function [Data,L]=CompileData(Image,Data,L)
        
        index=0;
        L.NumWells=length(Data.Index);
        L.NumSamples=sum(~Data.Remove);
        WAIT=waitbar(0,'Extracting data from positions...','Name','Data Extraction Percent Complete: ');
        disp('Extracting data from positions...')
        window=4*L.RadiusSolubilized+1;
        [MaskX,MaskY]=meshgrid(1:window,1:window);
        
        %Masks are made such that button is always centered
        %Chamber is then inserted into mask relative to button coordinates
        ButtonMask=uint16(sqrt((MaskX-(2*L.RadiusSolubilized+1)).^2+(MaskY-(2*L.RadiusSolubilized+1)).^2)<=L.Radius*.9);
        
        for W=1:L.NumWells
            index=index+double(~Data.Remove(W));
            if Data.Remove(W)==0
            Data.Index(W)=index;
            end
            
            %Generate masks for data extraction
            ChamberBGMask =       uint16(sqrt((MaskX-(Data.ChamberXCoor(W)-Data.ButtonsXCoor(W)+2*L.RadiusSolubilized+1)).^2+(MaskY-(Data.ChamberYCoor(W)-Data.ButtonsYCoor(W)+2*L.RadiusSolubilized+1)).^2)>=L.RadiusSolubilized*1.1 & sqrt((MaskX-(Data.ChamberXCoor(W)-Data.ButtonsXCoor(W)+2*L.RadiusSolubilized+1)).^2+(MaskY-(Data.ChamberYCoor(W)-Data.ButtonsYCoor(W)+2*L.RadiusSolubilized+1)).^2)<=L.RadiusSolubilized*1.3 );
            ChamberNoButtonMask = uint16(sqrt((MaskX-(Data.ChamberXCoor(W)-Data.ButtonsXCoor(W)+2*L.RadiusSolubilized+1)).^2+(MaskY-(Data.ChamberYCoor(W)-Data.ButtonsYCoor(W)+2*L.RadiusSolubilized+1)).^2)<=L.RadiusSolubilized &~ sqrt((MaskX-(2*L.RadiusSolubilized+1)).^2+(MaskY-(2*L.RadiusSolubilized+1)).^2)<=L.Radius*1.1 );
            Data.ButtonsAreaFG(W)=sum(sum(ButtonMask));
            Data.ButtonsAreaBG(W)=sum(sum(ChamberNoButtonMask));
            Data.ChamberAreaFG(W)=sum(sum(ChamberNoButtonMask));
            Data.ChamberAreaBG(W)=sum(sum(ChamberBGMask));
            
            %Collect data from solubilized chambers
            for frameS=1:L.numsolframes
                
                imageS=Image.solubilized((Data.ButtonsYCoor(W)-2*L.RadiusSolubilized):(Data.ButtonsYCoor(W)+2*L.RadiusSolubilized),(Data.ButtonsXCoor(W)-2*L.RadiusSolubilized):(Data.ButtonsXCoor(W)+2*L.RadiusSolubilized),frameS);
                DNAChamber = double(imageS.*ChamberNoButtonMask);
                DNAChamberBG = double(imageS.*ChamberBGMask);

                Data.solubilizedMedianFG(W,frameS)=median(DNAChamber(DNAChamber(:)>0));
                Data.solubilizedMeanFG(W,frameS)=mean(DNAChamber(DNAChamber(:)>0));
                Data.solubilizedSTDFG(W,frameS)=std(DNAChamber(DNAChamber(:)>0));
                Data.solubilizedTotalFG(W,frameS)=sum(DNAChamber(DNAChamber(:)>0));
                Data.solubilizedMedianBG(W,frameS)=median(DNAChamberBG(DNAChamberBG(:)>0));
                Data.solubilizedMeanBG(W,frameS)=mean(DNAChamberBG(DNAChamberBG(:)>0));
                Data.solubilizedSTDBG(W,frameS)=std(DNAChamberBG(DNAChamberBG(:)>0));
                Data.solubilizedTotalBG(W,frameS)=sum(DNAChamberBG(DNAChamberBG(:)>0))*Data.ChamberAreaFG(W)./Data.ChamberAreaBG(W);
                Data.solubilizedFractionSaturatedFG(W,frameS)=length(find(DNAChamber==65535))./length(find(DNAChamber(:)>0));
                Data.solubilizedFractionSaturatedBG(W,frameS)=length(find(DNAChamberBG==65535))./length(find(DNAChamberBG(:)>0));
                
            end
            
            %Collect data from surface immobilized molecules
            imageB=Image.surface((Data.ButtonsYCoor(W)-2*L.RadiusSolubilized):(Data.ButtonsYCoor(W)+2*L.RadiusSolubilized),(Data.ButtonsXCoor(W)-2*L.RadiusSolubilized):(Data.ButtonsXCoor(W)+2*L.RadiusSolubilized));
            SurfaceButton=double(imageB.*ButtonMask);
            SurfaceBG=double(imageB.*ChamberNoButtonMask);
            
            Data.surfaceMedianFG(W)=median(SurfaceButton(SurfaceButton(:)>0));
            Data.surfaceAverageFG(W)=mean(SurfaceButton(SurfaceButton(:)>0));
            Data.surfaceSTDFG(W)=std(SurfaceButton(SurfaceButton(:)>0));
            Data.surfaceTotalFG(W)=sum(SurfaceButton(SurfaceButton(:)>0));
            Data.surfaceMedianBG(W)=median(SurfaceBG(SurfaceBG(:)>0));
            Data.surfaceAverageBG(W)=mean(SurfaceBG(SurfaceBG(:)>0));
            Data.surfaceSTDBG(W)=std(SurfaceBG(SurfaceBG(:)>0));
            Data.surfaceTotalBG(W)=sum(SurfaceBG(SurfaceBG(:)>0))*Data.ButtonsAreaFG(W)./Data.ButtonsAreaBG(W);
            Data.surfaceFractionSaturatedFG(W)=length(find(SurfaceButton==65535))./length(find(SurfaceButton(:)>0));
            Data.surfaceFractionSaturatedBG(W)=length(find(SurfaceBG==65535))./length(find(SurfaceBG(:)>0));
            
            %Collect data from captured molecule images
            for frameC=1:L.numboundframes
                
                imageC=Image.captured((Data.ButtonsYCoor(W)-2*L.RadiusSolubilized):(Data.ButtonsYCoor(W)+2*L.RadiusSolubilized),(Data.ButtonsXCoor(W)-2*L.RadiusSolubilized):(Data.ButtonsXCoor(W)+2*L.RadiusSolubilized),frameC);
                CapturedButton=double(imageC.*ButtonMask);
                CapturedBG=double(imageC.*ChamberNoButtonMask);
                
                Data.capturedMedianFG(W,frameC)=median(CapturedButton(CapturedButton(:)>0));
                Data.capturedAverageFG(W,frameC)=mean(CapturedButton(CapturedButton(:)>0));
                Data.capturedSTDFG(W,frameC)=std(CapturedButton(CapturedButton(:)>0));
                Data.capturedTotalFG(W,frameC)=sum(CapturedButton(CapturedButton(:)>0));
                Data.capturedMedianBG(W,frameC)=median(CapturedBG(CapturedBG(:)>0));
                Data.capturedAverageBG(W,frameC)=mean(CapturedBG(CapturedBG(:)>0));
                Data.capturedSTDBG(W,frameC)=std(CapturedBG(CapturedBG(:)>0));
                Data.capturedTotalBG(W,frameC)=sum(CapturedBG(CapturedBG(:)>0))*Data.ButtonsAreaFG(W)./Data.ButtonsAreaBG(W); 
                Data.capturedFractionSaturatedFG(W,frameC)=length(find(CapturedButton==65535))./length(find(CapturedButton(:)>0));
                Data.capturedFractionSaturatedBG(W,frameC)=length(find(CapturedBG==65535))./length(find(CapturedBG(:)>0));
                
            end
            
            waitbar(W/L.NumWells,WAIT,sprintf('%6.3f',W/L.NumWells*100));

        end

        disp('Data Extraction Complete.')
        disp('Preparing to save data.')
        delete(WAIT)
    end


%% PRINT DATA FUNCTION
    function []=fprintMITOMI(L,Data)
        
        savemat=strcat('editted_',L.name,'_AnalysisData_v2_2.mat');
        save(savemat,'L','Data')

        %create string header for dissociation data
        HeaderFormat={'Index','ColIndx','RowIndx','Removed','Flagged','ButXCoor','ButYCoor','ButRad','BuAreaFG','BuAreaBG','ButAutoF','BNDMedFG','BNDAvgFG','BNDStdFG','BNDSumFG','BNDSatFG','BNDMedBG','BNDAvgBG','BNDStdBG','BNDSumBG','BNDSatBG'};
        for z=1:L.numboundframes
            HeaderFormat(21+z)={['CAPMedFG' num2str(z)]};
            HeaderFormat(21+z+L.numboundframes)={['CAPAvgFG' num2str(z)]};
            HeaderFormat(21+z+L.numboundframes*2)={['CAPStdFG' num2str(z)]};
            HeaderFormat(21+z+L.numboundframes*3)={['CAPSumFG' num2str(z)]};
            HeaderFormat(21+z+L.numboundframes*4)={['CAPSatFG' num2str(z)]};
            HeaderFormat(21+z+L.numboundframes*5)={['CAPMedBG' num2str(z)]};
            HeaderFormat(21+z+L.numboundframes*6)={['CAPAvgBG' num2str(z)]};
            HeaderFormat(21+z+L.numboundframes*7)={['CAPStdBG' num2str(z)]};
            HeaderFormat(21+z+L.numboundframes*8)={['CAPSumBG' num2str(z)]};
            HeaderFormat(21+z+L.numboundframes*9)={['CAPSatBG' num2str(z)]};
        end

        HeaderFormat(end+1:end+6)={'SOLXCoor','SOLYCoor','SOLRad','SOAreaFG','SOAreaBG','SOLAutoF'};
        HOffNum=27+10*L.numboundframes;

        for SolFrame=1:L.numsolframes
            HeaderFormat(HOffNum+SolFrame)={['SOLMedFG' num2str(SolFrame)]};
            HeaderFormat(HOffNum+SolFrame+L.numsolframes)={['SOLAvgFG' num2str(SolFrame)]};
            HeaderFormat(HOffNum+SolFrame+L.numsolframes*2)={['SOLStdFG' num2str(SolFrame)]};
            HeaderFormat(HOffNum+SolFrame+L.numsolframes*3)={['SOLSumFG' num2str(SolFrame)]};
            HeaderFormat(HOffNum+SolFrame+L.numsolframes*4)={['SOLSatFG' num2str(SolFrame)]};
            HeaderFormat(HOffNum+SolFrame+L.numsolframes*5)={['SOLMedBG' num2str(SolFrame)]};
            HeaderFormat(HOffNum+SolFrame+L.numsolframes*6)={['SOLAvgBG' num2str(SolFrame)]};
            HeaderFormat(HOffNum+SolFrame+L.numsolframes*7)={['SOLStdBG' num2str(SolFrame)]};
            HeaderFormat(HOffNum+SolFrame+L.numsolframes*8)={['SOLSumBG' num2str(SolFrame)]};
            HeaderFormat(HOffNum+SolFrame+L.numsolframes*9)={['SOLSatBG' num2str(SolFrame)]};
        end

        DataText=fopen(['editted_' L.name '_AnalysisData_v2_2.txt'],'w');
        IntermediateFormat=struct2cell(Data);
        DataFormat=cat(2,IntermediateFormat{:});
        fprintf(DataText,'%s\t',HeaderFormat{:} );
        fprintf(DataText,'\r\n');

        for Z=1:L.NumWells;
            fprintf(DataText,'%.0f\t',DataFormat(Z,:)');
            fprintf(DataText,'\r\n');
        end
        fclose(DataText);

        savetxt=strcat('editted_',L.name,'_AnalysisData_v2_2_NoHeaders.txt');
        save(savetxt,'DataFormat','-ascii')
        disp('Data saved. Data extraction complete.')    

    end

end

%% NOTEPAD
%{
v2.2 restricts the search area of the solubility chamber by 1/8th of the
length of each direction away from the center point. It also includes
columns for 

v2.1 supports multiple prewash and postwash images via the "Custom Multi
Image" option. All analyses now include saturation percentages for all FG
and BG images. This is an indicator for aggregates in the BG is as well as how
unreliable a FG is.

v2.0 supports dissociation data images as well as equilibrium. GUI controls
updated to be compatible with all screen sizes as well as Windows/Mac OS.
Data extraction now uses whole chamber - button for button BG as well as  
chamber FG. Chamber BG is now a ring outside of the chamber. Data
extraction is now performed after user confirms positions for data
extraction. Dual channel images are now used during the coordinate setting
stage to help the user visualize the extent of their printed DNA array. In
addition, all functions are embedded - therefore there is now no need for
external files.

v1.0 supports equilibrium data extraction. User sets coodinates for the
vertices of their printed array and program attempts to automatically find
features within the expected location. The first attempt is using a hough
transform via imfindcircles. if that fails, the program finds a location
with the brightest intensity and centers around that. Data extraction
occurs after the feature location is determined. The user can then redefine
positions if a feature was incorrectly identified.

Robert Puccinelli
Fordyce Lab 20160309
rpuccinelli@stanford.edu
robert.puccinelli@outlook.com

%}