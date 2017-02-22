function MITOMI_Function=MITOMIAnalysis()
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

%Core Data Structures
global Log
Log=[];
Image=[];
Data=[];
ME=MException('','');

%Core Operations
try
    LogStructureInitialization();
    MITOMIAnalysis_Initialization();
    AnalysisValidation();
    MITOMIAnalysis_FileManager();
    ImagePrep();
    Image=MITOMIAnalysis_ImageManipulation(Image.Captured(:,:,1),Image);
    MITOMIAnalysis_SetCoordinates(Image.Surface);
    
    if ~isempty(Log.BackgroundFrame)
        MITOMIAnalysis_SetCoordinates((Image.Background+Image.Surface)/2);       
    end
    
    DataStructureInitialization();
    AutomatedFeatureFinding();
    MITOMIAnalysis_UserEdit(Log,Image,Data);
    abortMITOMI()
    
catch ME
    
    warning('off','backtrace')
    abortMITOMI()
    warning(ME.message)
    warning('on','backtrace')
    
end

%% LOG INITIALIZATION
    function []=LogStructureInitialization()
        
        Log.LastUpdated='20170207';
        Log.Error=[];
        
        %Progress check
        Log.Initialization=[];
        Log.FileManager=[];
        Log.ImageManipulation=[];
        Log.SetCoordinatesSurface=[];
        Log.SetCoordinatesBackground=[];
        Log.SurfaceFeatureFinding=[];
        Log.BackgroundFeatureFinding=[];
        Log.SurfaceFeatureEvaluation=[];
        Log.CapturedFeatureEvaluation=[];
        Log.BackgroundFeatureEvaluation=[];
        Log.InitiatedFileSave=[];
        
        %Core variables
        Log.Directory=[];
        Log.AnalysisType=[];
        Log.RadiusType=[];
        Log.Background=[];
        Log.NumberFrames=[];
        Log.NameFrames=[];
        Log.NameOutput=[];

        Log.BackgroundFrame=[];
        Log.SurfaceFrame=[];
        Log.CapturedFrames=[];
        
        Log.Rows=[];
        Log.Cols=[];        
        
        
        Log.ApproxGammaSurface=[];
        Log.ApproxButtonRadius=[];
        Log.ApproxBackgroundRadius=[];
        Log.CoorButtons=[]; 
        Log.CoorBackground=[];
        
        Log.SubimageButtonRadius=[];
        Log.ButtonFGmask=[];
        Log.ButtonBGmask=[];
        Log.ButtonTicker=[];

        Log.SolubilizedFGmask=[];
        Log.SolubilizedBGmask=[];
        Log.BackgroundTicker=[];
        
        %Internal GUI tracking variables
        Log.Vertex=0; 
        Log.CoorList=[0,0];
        Log.CurrView=[];
        Log.MagView=[];
        Log.ManipulationAPI=[];
        Log.TempGamma=[];
        Log.Gamma=[];
        Log.RadiusSample=[];

        
    end

%% DATA STRUCTURE INITIALIZATION
    function []=DataStructureInitialization()
        
        dimensions=zeros(Log.Rows*Log.Cols,1);
        
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
        
        Data.SurfaceMedianFG=dimensions;
        Data.SurfaceAverageFG=dimensions;
        Data.SurfaceSTDFG=dimensions;
        Data.SurfaceTotalFG=dimensions;
        Data.SurfaceFractionSaturatedFG=dimensions;
        Data.SurfaceMedianBG=dimensions;
        Data.SurfaceAverageBG=dimensions;
        Data.SurfaceSTDBG=dimensions;
        Data.SurfaceTotalBG=dimensions;
        Data.SurfaceFractionSaturatedBG=dimensions;
        
        Data.CapturedMedianFG=dimensions;
        Data.CapturedAverageFG=dimensions;
        Data.CapturedSTDFG=dimensions;
        Data.CapturedTotalFG=dimensions;
        Data.CapturedFractionSaturatedFG=dimensions;
        Data.CapturedMedianBG=dimensions;
        Data.CapturedAverageBG=dimensions;
        Data.CapturedSTDBG=dimensions;
        Data.CapturedTotalBG=dimensions;
        Data.CapturedFractionSaturatedBG=dimensions;
       
        Data.ChamberXCoor=dimensions;
        Data.ChamberYCoor=dimensions;
        Data.ChamberRadius=dimensions;
        Data.ChamberAreaFG=dimensions;
        Data.ChamberAreaBG=dimensions;
        Data.AutofindChambers=dimensions;
        
        Data.SolubilizedMedianFG=dimensions;
        Data.SolubilizedMeanFG=dimensions;
        Data.SolubilizedSTDFG=dimensions;
        Data.SolubilizedTotalFG=dimensions;
        Data.SolubilizedFractionSaturatedFG=dimensions;
        Data.SolubilizedMedianBG=dimensions;
        Data.SolubilizedMeanBG=dimensions;
        Data.SolubilizedSTDBG=dimensions;
        Data.SolubilizedTotalBG=dimensions;
        Data.SolubilizedFractionSaturatedBG=dimensions;
        
    end
%% ANALYSIS VALIDATION FUNCTION
    function []=AnalysisValidation()

        try
            %Check to see if GUI stage passed parameters
            assert(~isempty(Log.Initialization),'MITOMIAnalysis:MITOMIAnalysis_Initialization:usrCancel','User cancelled GUI operation');

            %Validate passed parameters
            switch [Log.AnalysisType Log.Background]

                case 'EquilibriumYes'
                    assert(Log.NumberFrames==3,'MITOMIAnalysis:ImagePrep:NumberFrames','Expected 3 frames for equilibrium analysis with background, received %i',Log.NumberFrames);
                    Log.BackgroundFrame=1;
                    Log.SurfaceFrame=1;
                    Log.CapturedFrames=1;

                case 'EquilibriumNo'
                    assert(Log.NumberFrames==2,'MITOMIAnalysis:ImagePrep:NumberFrames','Expected 2 frames for equilibrium analysis without background, received %i',Log.NumberFrames);
                    Log.BackgroundFrame=[];
                    Log.SurfaceFrame=1;
                    Log.CapturedFrames=1;

                case 'DissociationYes'
                    assert(Log.NumberFrames>3,'MITOMIAnalysis:ImagePrep:NumberFrames','Expected more than 3 frames for dissociation analysis with background, received %i',Log.NumberFrames);
                    Log.BackgroundFrame=1;
                    Log.SurfaceFrame=1;
                    Log.CapturedFrames=Log.NumberFrames-2;

                case 'DissociationNo'
                    assert(Log.NumberFrames>2,'MITOMIAnalysis:ImagePrep:NumberFrames','Expected more than 2 frames for dissociation analysis without background, received %i',Log.NumberFrames);
                    Log.BackgroundFrame=[];
                    Log.SurfaceFrame=1;
                    Log.CapturedFrames=Log.NumberFrames-1;

            end

        catch ME

            Log.Error=getReport(ME);
            abortMITOMI();
            throw(ME)

        end
    end

%% IMAGEPREP FUNCTION
    function []=ImagePrep()
    
        %Initialize image structure
        Image.Background=[];
        Image.Surface=[];
        Image.Captured=[];
        
        try
            %Check to see if GUI stage passed parameters
            assert(~isempty(Log.FileManager),'MITOMIAnalysis:MITOMIAnalysis_FileManager:usrCancel','User cancelled GUI operation');

            %Load images and update waitbar
            WAIT=waitbar(0,'Loading images into MATLAB...','Name','Loading Images ','CreateCancelBtn','setappdata(gcbf,''canceling'',1);');
            setappdata(WAIT,'canceling',0);
            waitbar(0/Log.NumberFrames,WAIT,sprintf('Fraction Complete: %i / %i',0,Log.NumberFrames));
            try
                Image.Background=imread(Log.NameFrames{Log.BackgroundFrame});
                dimBG=size(Image.Background); 
                Image.Background=imresize(Image.Background,7500/min(dimBG));
                bgFilled=1;
            catch
                Image.Background=[];
                bgFilled=~isempty(Image.Background);
            end

            waitbar(bgFilled/Log.NumberFrames,WAIT,sprintf('Fraction Complete: %i / %i',bgFilled,Log.NumberFrames));

            Image.Surface=imread(Log.NameFrames{~isempty(Log.BackgroundFrame) + Log.SurfaceFrame});
            dimSurface=size(Image.Surface);
            Image.Surface=imresize(Image.Surface,7500/min(dimSurface));


            waitbar((bgFilled+1)/Log.NumberFrames,WAIT,sprintf('Fraction Complete: %i / %i',bgFilled+1,Log.NumberFrames));

            for i = (bgFilled+Log.SurfaceFrame+1):(bgFilled+Log.SurfaceFrame+Log.CapturedFrames)

                %Check to see if cancel button has been pressed, throw error
                if isequal(getappdata(WAIT,'canceling'),1)
                    delete(WAIT)
                    assert(false,'MITOMIAnalysis:ImagePrep:waitCancel','User cancelled image loading operation');
                end

                tempCaptured(:,:)=imread(Log.NameFrames{i});
                dimCaptured=size(tempCaptured);
                Image.Captured(:,:,i-bgFilled-Log.SurfaceFrame)=imresize(tempCaptured,7500/min(dimCaptured));
                tempCaptured=[];

                waitbar(i/Log.NumberFrames,WAIT,sprintf('Fraction Complete: %i / %i',i,Log.NumberFrames));
            end

            delete(WAIT);

            %Check to see if all images are the same dimension
            if bgFilled
                assert(isequal(dimSurface,dimCaptured,dimBG),'MITOMIAnalysis:ImagePrep:dimImages','Image dimensions do not match')
            else
                assert(isequal(dimSurface,dimCaptured),'MITOMIAnalysis:ImagePrep:dimImages','Image dimensions do not match')
                Image.Background=zeros(dimSurface);
            end

            %Scale image size to prevent Matlab from crashing when disp large image  

        catch ME

            Log.Error=getReport(ME);
            abortMITOMI();
            throw(ME)

        end
    end

%% AUTOMATED FEATURE FINDING FUNCTION
    function []=AutomatedFeatureFinding()
       
        buttonMSGBOX=[];
        backgroundMSGBOX=[];
        
        %Adjust radii for subimage mask-making
        Log.SubimageButtonRadius=Log.ApproxButtonRadius+round(Log.ApproxButtonRadius/2);

        %Masks dimensions for surface image
        buttonFGmask=zeros(Log.SubimageButtonRadius*2);
        buttonBGmask=zeros(Log.SubimageButtonRadius*2);
        
        %Generate FG and BG masks
        for i=1:Log.SubimageButtonRadius*2
            for j=1:Log.SubimageButtonRadius*2
                
                if lt((i-(Log.SubimageButtonRadius*2+1)/2)^2+(j-(Log.SubimageButtonRadius*2+1)/2)^2,(Log.SubimageButtonRadius/2)^2)
                    buttonFGmask(i,j)=1;
                else
                    buttonFGmask(i,j)=0;
                end
                
                if lt((i-(Log.SubimageButtonRadius*2+1)/2)^2+(j-(Log.SubimageButtonRadius*2+1)/2)^2,(Log.SubimageButtonRadius*2/2)^2) && gt((i-(Log.SubimageButtonRadius*2+1)/2)^2+(j-(Log.SubimageButtonRadius*2+1)/2)^2,(Log.SubimageButtonRadius*1.5/2)^2)
                    buttonBGmask(i,j)=1;
                else
                    buttonBGmask(i,j)=0;
                end
                
            end
        end
        
        %Save button masks and generate autofind counter
        Log.ButtonFGmask=buttonFGmask;
        Log.ButtonBGmask=buttonBGmask;
        Log.ButtonTicker=0;
        buttonRadius=zeros(length(Data.Index),1);
        
        figButtonGrid=figure('menubar','none','numbertitle','off','toolbar','none','Name','Button Preview');
        WAIT=waitbar(0,'Processing button positions...','Name','Finding Buttons');
        
        for m=1:length(Data.Index); %For each point in the array
            
            %Fill in spot identity
            Data.ColIndex(m)=ceil(m/Log.Rows);
            Data.RowIndex(m)=m-Log.Rows*(Data.ColIndex(m)-1);
            
            %Set approximate X and Y coordinates
            CoorX=Log.CoorButtons(m,1);
            CoorY=Log.CoorButtons(m,2);
            
            %Generate subimage and adjust contrast
            screenSurface=double(Image.Surface((CoorY-2*Log.SubimageButtonRadius):(CoorY+2*Log.SubimageButtonRadius),(CoorX-2*Log.SubimageButtonRadius):(CoorX+2*Log.SubimageButtonRadius)));
            medianSS=median(median(screenSurface));
            stdSS=std(std(screenSurface));
            screenSurfaceMod=uint16(mat2gray(screenSurface,[medianSS+stdSS medianSS+stdSS*4])*65535);
            imshow(screenSurfaceMod)
            
            %Apply Hough transform to find button
            warning('OFF','all') %suppress small radius warning
            [spotLocations,radii]=imfindcircles((screenSurfaceMod),[round(Log.ApproxButtonRadius/2) round(Log.ApproxButtonRadius/0.85)],'ObjectPolarity','bright');
            warning('ON','all')
            
            if ~isempty(radii) %If autofind with hough transform detects something, process it
                
                %Convert local coordinates to global coordinates
                Data.ButtonsXCoor(m,1)=round(spotLocations(1,1)-Log.SubimageButtonRadius*2-1+CoorX);
                Data.ButtonsYCoor(m,1)=round(spotLocations(1,2)-Log.SubimageButtonRadius*2-1+CoorY);
                Data.Remove(m)=false;
                Data.AutofindButtons(m)=true;
                buttonRadius(m,1)=radii(1);
                
                %Increment button autofind ticker
                Log.ButtonTicker=Log.ButtonTicker+1;
                
            else %Autofind failed, try to find where button is with a scan
                surfaceFGdataholder=cell(length(Log.SubimageButtonRadius*4+1));
                surfaceBGdataholder=cell(length(Log.SubimageButtonRadius*4+1));
                
                %Generate map of net local intensities for subimage
                for n=(CoorX-Log.SubimageButtonRadius*2):(CoorX+Log.SubimageButtonRadius*2)
                    for o=(CoorY-Log.SubimageButtonRadius*2):(CoorY+Log.SubimageButtonRadius*2)
                        
                        ExtractImageSurface=double(Image.Surface((o-Log.SubimageButtonRadius):(o+Log.SubimageButtonRadius-1),(n-Log.SubimageButtonRadius):(n+Log.SubimageButtonRadius-1)));
                        surfaceFGsampletemp=ExtractImageSurface.*Log.ButtonFGmask;
                        surfaceBGsampletemp=ExtractImageSurface.*Log.ButtonBGmask;
                        surfaceFGdataholder{n+Log.SubimageButtonRadius*2-CoorX+1,o+Log.SubimageButtonRadius*2-CoorY+1}=surfaceFGsampletemp(surfaceFGsampletemp>0);
                        surfaceBGdataholder{n+Log.SubimageButtonRadius*2-CoorX+1,o+Log.SubimageButtonRadius*2-CoorY+1}=surfaceBGsampletemp(surfaceBGsampletemp>0);
                        
                    end
                end
                
                %Find data point w highest net int in phase space datasets
                NetInt=cellfun(@(x,y) (sum(x)-sum(y)),surfaceFGdataholder,surfaceBGdataholder);
                [n_local,o_local]=find(NetInt==max(NetInt(:)));

                %Convert "best data" local coor into global for surface image
                Data.ButtonsXCoor(m)=n_local-1+CoorX-Log.SubimageButtonRadius*2;
                Data.ButtonsYCoor(m)=o_local-1+CoorY-Log.SubimageButtonRadius*2;   
                Data.Remove(m)=false;
                Data.AutofindButtons(m)=false;
                buttonRadius(m)=Log.ApproxButtonRadius;
               
            end
            waitbar(m/length(Data.Index));
        end
        
        close(figButtonGrid)
        delete(WAIT)
        
        %Set button radius length
        if ~strcmp(Log.RadiusType,'Fixed')
            Data.ButtonsRadius(:,1)=buttonRadius(:);
        else
            Data.ButtonsRadius(:,1)=ones(length(Data.Index),1)*Log.ApproxButtonRadius;
        end
        
        %Display quality of autofind for buttons
        buttonMSGBOX=msgbox(['Buttons identified with automation: ' num2str(Log.ButtonTicker) ' out of ' num2str(length(Data.Index))],'Button Detection Complete');            
        movegui(buttonMSGBOX,'north');

        if ~isempty(Log.BackgroundFrame)

            solubilizedFGmask=zeros(Log.ApproxBackgroundRadius*2);
            solubilizedBGmask=ones(Log.ApproxBackgroundRadius*2);

            %Generate 
            for p=1:Log.ApproxBackgroundRadius*2
                for q=1:Log.ApproxBackgroundRadius*2
                    if lt((p-(Log.ApproxBackgroundRadius*2+1)/2)^2+(q-(Log.ApproxBackgroundRadius*2+1)/2)^2,Log.ApproxBackgroundRadius^2)
                        solubilizedFGmask(p,q)=1;
                    end
                end
            end

            Log.SolubilizedFGmask=solubilizedFGmask;
            Log.SolubilizedBGmask=solubilizedBGmask-solubilizedFGmask;
            Log.BackgroundTicker=0;            

            figChamberGrid=figure('menubar','none','numbertitle','off','toolbar','none','Name','Chamber Preview');
            WAIT=waitbar(0,'Processing chamber positions...','Name','Chamber Positions');
            Log.SurfaceFeatureFinding='Passed';
        
            for n=1:length(Data.Index)

                %Refresh coordinates with chamber positions
                CoorX=Log.CoorBackground(n,1);
                CoorY=Log.CoorBackground(n,2);
                Data.ChamberRadius(n,1)=Log.ApproxBackgroundRadius;

                %Generate subimage and adjust contrast
                screenSol=double(Image.Background((CoorY-Log.ApproxBackgroundRadius):(CoorY+Log.ApproxBackgroundRadius),(CoorX-Log.ApproxBackgroundRadius):(CoorX+Log.ApproxBackgroundRadius),1));
                screenSolSTD=std(screenSol(:));
                screenSolMED=median(screenSol(:));
                screenSolMod=imadjust(uint16(mat2gray(screenSol,[screenSolMED-screenSolSTD*2, screenSolMED+screenSolSTD*2])*65535));
                imshow(screenSolMod)

                %Apply Hough transform to identify chamber
                [backgroundLocations,backgroundRadii]=imfindcircles((screenSolMod),[round(Log.ApproxBackgroundRadius*.80) round(Log.ApproxBackgroundRadius*1.2)],'ObjectPolarity','bright');

                if ~isempty(backgroundRadii) %if autofind with hough transform finds something, process 

                    %Convert local coordinates to global coordinates
                    Data.ChamberXCoor(n)=round(backgroundLocations(1,1)-Log.ApproxBackgroundRadius-1+CoorX);
                    Data.ChamberYCoor(n)=round(backgroundLocations(1,2)-Log.ApproxBackgroundRadius-1+CoorY);
                    Data.AutofindChambers(n)=true;

                    Log.BackgroundTicker=Log.BackgroundTicker+1;

                else %autofind failed
                    solubilizedFGdataholder=cell(length(Log.ApproxBackgroundRadius*2+1));
                    solubilizedBGdataholder=cell(length(Log.ApproxBackgroundRadius*2+1));

                    %Generate map of net local intensities for subimage
                    for s=(CoorX-ceil(Log.ApproxBackgroundRadius*7/8)):(CoorX+floor((Log.ApproxBackgroundRadius-1)*7/8))
                        for t=(CoorY-ceil(Log.ApproxBackgroundRadius*7/8)):(CoorY+floor((Log.ApproxBackgroundRadius-1)*7/8))
                            
                            ExtractImageSolubilized=double(Image.Background((t-Log.ApproxBackgroundRadius):(t+Log.ApproxBackgroundRadius-1),(s-Log.ApproxBackgroundRadius):(s+Log.ApproxBackgroundRadius-1),1));
                            solubilizedFGsampletemp=ExtractImageSolubilized.*Log.SolubilizedFGmask;
                            solubilizedBGsampletemp=ExtractImageSolubilized.*Log.SolubilizedBGmask;
                            solubilizedFGdataholder{s+Log.ApproxBackgroundRadius-CoorX+1,t+Log.ApproxBackgroundRadius-CoorY+1}=solubilizedFGsampletemp(solubilizedFGsampletemp>0);                             
                            solubilizedBGdataholder{s+Log.ApproxBackgroundRadius-CoorX+1,t+Log.ApproxBackgroundRadius-CoorY+1}=solubilizedBGsampletemp(solubilizedBGsampletemp>0);                             

                        end
                    end

                    %Find data point w highest net int in phase space datasets
                    NetInt=cellfun(@(x,y) (sum(x)-sum(y)),solubilizedFGdataholder,solubilizedBGdataholder);
                    [s_local,t_local]=find(NetInt==max(NetInt(:)));

                    %Convert "best data" local coor into global for capt image
                    Data.ChamberXCoor(n)=s_local(1)-1+CoorX-Log.ApproxBackgroundRadius;
                    Data.ChamberYCoor(n)=t_local(1)-1+CoorY-Log.ApproxBackgroundRadius;
                    Data.AutofindChambers(n)=false;
                end

                waitbar(n/length(Data.Index));

            end

            close(figChamberGrid)
            delete(WAIT)
            
            Log.BackgroundFeatureFinding='Passed';

            backgroundMSGBOX=msgbox(['Background chambers identified with automation: ' num2str(Log.BackgroundTicker) ' out of ' num2str(length(Data.Index))],'Background Detection Complete');
            movegui(backgroundMSGBOX,'north');
        end
        
        pause(3.0)
        
        %Close figures if they remain open
        if ishandle(buttonMSGBOX)
            close(buttonMSGBOX)
        end
        
        if ishandle(backgroundMSGBOX)
            close(backgroundMSGBOX)
        end
    end

%% ABORT MITOMI
    function []=abortMITOMI()

    if ~strcmp(ME.identifier,'MITOMIAnalysis:MITOMIAnalysis_Initialization:usrCancel')
        LOGFILENAME=['LOG_MITOMIAnalysis_' datestr(now,30) '.mat'];
        save(LOGFILENAME,'Log')
        disp('Log file saved.')
        clear global
        close()
    elseif strcmp(Log.InitiatedFileSave,'Passed')
        clear global
        disp('Congrats - MITOMI Analysis is complete')
    end
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