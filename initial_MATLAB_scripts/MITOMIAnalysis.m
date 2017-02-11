function []=MITOMIAnalysis()
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
    MITOMIAnalysis_ImageManipulation(Image.Captured(:,:,1),Image);
    MITOMIAnalysis_SetCoordinates(Image.Surface);
    
    if ~isempty(Log.BackgroundFrame)
        MITOMIAnalysis_SetCoordinates((Image.Background+Image.Surface)/2);       
    end
    
    DataStructureInitialization();
    AutomatedFeatureFinding();
    MITOMIAnalysis_UserEdit(Image,Data);
    CompileData();
    fprintMITOMI();
    
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
        Log.ImageHolder=[];
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

            Image.Background=imread(Log.NameFrames{Log.BackgroundFrame});
            bgFilled=~isempty(Image.Background);
            waitbar(bgFilled/Log.NumberFrames,WAIT,sprintf('Fraction Complete: %i / %i',bgFilled,Log.NumberFrames));

            Image.Surface=imread(Log.NameFrames{~isempty(Log.BackgroundFrame) + Log.SurfaceFrame});
            waitbar((bgFilled+1)/Log.NumberFrames,WAIT,sprintf('Fraction Complete: %i / %i',bgFilled+1,Log.NumberFrames));

            for i = (bgFilled+Log.SurfaceFrame+1):(bgFilled+Log.SurfaceFrame+Log.CapturedFrames)

                %Check to see if cancel button has been pressed, throw error
                if isequal(getappdata(WAIT,'canceling'),1)
                    delete(WAIT)
                    assert(false,'MITOMIAnalysis:ImagePrep:waitCancel','User cancelled image loading operation');
                end

                Image.Captured(:,:,i-bgFilled-Log.SurfaceFrame)=imread(Log.NameFrames{i});
                waitbar(i/Log.NumberFrames,WAIT,sprintf('Fraction Complete: %i / %i',i,Log.NumberFrames));
            end

            delete(WAIT);

            %Check to see if all images are the same dimension
            dimSurface=size(Image.Surface);
            if bgFilled
                assert(isequal(dimSurface,size(Image.Captured(:,:,1)),size(Image.Background)),'MITOMIAnalysis:ImagePrep:dimImages','Image dimensions do not match')
            else
                assert(isequal(dimSurface,size(Image.Captured(:,:,1))),'MITOMIAnalysis:ImagePrep:dimImages','Image dimensions do not match')
                Image.Background=zeros(dimSurface);
            end

            %Scale image size to prevent Matlab from crashing when disp large image  
            Image.Background=imresize(Image.Background,7500/min(dimSurface));
            Image.Surface=imresize(Image.Surface,7500/min(dimSurface));
            Image.Captured=imresize(Image.Captured,7500/min(dimSurface));

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
            screenSurface=uint16(Image.Surface((CoorY-2*Log.SubimageButtonRadius):(CoorY+2*Log.SubimageButtonRadius),(CoorX-2*Log.SubimageButtonRadius):(CoorX+2*Log.SubimageButtonRadius)));
            screenSurfaceMod=imadjust(screenSurface,[],[],Log.ApproxGammaSurface);
            imshow(mat2gray(screenSurfaceMod))
            
            %Apply Hough transform to find button
            warning('OFF','all') %suppress small radius warning
            [spotLocations,radii]=imfindcircles((screenSurfaceMod),[round(Log.SubimageButtonRadius/2.5) round(Log.SubimageButtonRadius/1.25)],'ObjectPolarity','bright');
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
            Data.ButtonRadius(:,1)=buttonRadius(:);
        else
            Data.ButtonsRadius(:,1)=ones(length(Data.Index),1)*Log.ApproxButtonRadius;
        end
        
        %Display quality of autofind for buttons
        buttonMSGBOX=msgbox(['Buttons identified with automation: ' num2str(Log.ButtonTicker) ' out of ' num2str(length(Data.Index))],'Button Detection Complete');            
        
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
        
        if ~isempty(Log.BackgroundFrame)

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
                    Data.AutofindChamber(n)=true;

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
                    Data.AutofindChamber(n)=false;
                end

                waitbar(n/length(Data.Index));

            end

            close(figChamberGrid)
            delete(WAIT)
            
            Log.BackgroundFeatureFinding='Passed';

            backgroundMSGBOX=msgbox(['Background chambers identified with automation: ' num2str(Log.BackgroundTicker) ' out of ' num2str(length(Data.Index))],'Background Detection Complete');
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
%% USER EDIT FUNCTION
    function []=UserEdit()
        
        %Initialize variables
        ImPreviewButtons=((1-mat2gray(imadjust(Image.Surface)))*3/4+mat2gray(Image.Background(:,:,1))/4);
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
            warning('OFF','all')
            if surmenu==0
                apiControl=scrollImage(ImWithAutoMiss);
            else
                apiControl.replaceImage(ImWithAutoMiss,'PreserveView',1);
            end
            warning('ON','all')
            
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
        ImPreviewBound=1-mat2gray(imadjust(Image.Captured(:,:,1)));

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
            warning('OFF','all')
            if bndmenu==0
                apiControl=scrollImage(ImWithFeatures);
            else
                apiControl.replaceImage(ImWithFeatures,'PreserveView',1);
            end
            warning('ON','all')
            
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
        ImPreviewChamber=1-mat2gray(imadjust(Image.Background(:,:,1)));
        
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
            warning('OFF','all')
            if solmenu==0
                apiControl=scrollImage(ImWithAMChambers);
            else
                apiControl.replaceImage(ImWithAMChambers,'PreserveView',1);
            end
            warning('ON','all')
            
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
        window=4*L.ApproxBackgroundRadius+1;
        [MaskX,MaskY]=meshgrid(1:window,1:window);
        
        %Masks are made such that button is always centered
        %Chamber is then inserted into mask relative to button coordinates
        ButtonMask=uint16(sqrt((MaskX-(2*L.ApproxBackgroundRadius+1)).^2+(MaskY-(2*L.ApproxBackgroundRadius+1)).^2)<=L.Radius*.9);
        
        for W=1:L.NumWells
            index=index+double(~Data.Remove(W));
            if Data.Remove(W)==0
            Data.Index(W)=index;
            end
            
            %Generate masks for data extraction
            ChamberBGMask =       uint16(sqrt((MaskX-(Data.ChamberXCoor(W)-Data.ButtonsXCoor(W)+2*L.ApproxBackgroundRadius+1)).^2+(MaskY-(Data.ChamberYCoor(W)-Data.ButtonsYCoor(W)+2*L.ApproxBackgroundRadius+1)).^2)>=L.ApproxBackgroundRadius*1.1 & sqrt((MaskX-(Data.ChamberXCoor(W)-Data.ButtonsXCoor(W)+2*L.ApproxBackgroundRadius+1)).^2+(MaskY-(Data.ChamberYCoor(W)-Data.ButtonsYCoor(W)+2*L.ApproxBackgroundRadius+1)).^2)<=L.ApproxBackgroundRadius*1.3 );
            ChamberNoButtonMask = uint16(sqrt((MaskX-(Data.ChamberXCoor(W)-Data.ButtonsXCoor(W)+2*L.ApproxBackgroundRadius+1)).^2+(MaskY-(Data.ChamberYCoor(W)-Data.ButtonsYCoor(W)+2*L.ApproxBackgroundRadius+1)).^2)<=L.ApproxBackgroundRadius &~ sqrt((MaskX-(2*L.ApproxBackgroundRadius+1)).^2+(MaskY-(2*L.ApproxBackgroundRadius+1)).^2)<=L.Radius*1.1 );
            Data.ButtonsAreaFG(W)=sum(sum(ButtonMask));
            Data.ButtonsAreaBG(W)=sum(sum(ChamberNoButtonMask));
            Data.ChamberAreaFG(W)=sum(sum(ChamberNoButtonMask));
            Data.ChamberAreaBG(W)=sum(sum(ChamberBGMask));
            
            %Collect data from solubilized chambers
            for frameS=1:L.numsolframes
                
                imageS=Image.solubilized((Data.ButtonsYCoor(W)-2*L.ApproxBackgroundRadius):(Data.ButtonsYCoor(W)+2*L.ApproxBackgroundRadius),(Data.ButtonsXCoor(W)-2*L.ApproxBackgroundRadius):(Data.ButtonsXCoor(W)+2*L.ApproxBackgroundRadius),frameS);
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
            imageB=Image.surface((Data.ButtonsYCoor(W)-2*L.ApproxBackgroundRadius):(Data.ButtonsYCoor(W)+2*L.ApproxBackgroundRadius),(Data.ButtonsXCoor(W)-2*L.ApproxBackgroundRadius):(Data.ButtonsXCoor(W)+2*L.ApproxBackgroundRadius));
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
                
                imageC=Image.captured((Data.ButtonsYCoor(W)-2*L.ApproxBackgroundRadius):(Data.ButtonsYCoor(W)+2*L.ApproxBackgroundRadius),(Data.ButtonsXCoor(W)-2*L.ApproxBackgroundRadius):(Data.ButtonsXCoor(W)+2*L.ApproxBackgroundRadius),frameC);
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
            WAIT=waitbar(0,'Extracting data from positions...','Name','Data Extraction Percent Complete: ');
        
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

function []=abortMITOMI()

    if ~strcmp(ME.identifier,'MITOMIAnalysis:MITOMIAnalysis_Initialization:usrCancel')
        Log.ImageHolder='Images not saved';
        LOGFILENAME=['LOG_MITOMIAnalysis_' datestr(now,30) '.mat'];
        save(LOGFILENAME,'Log')
        disp('Log file saved.')
        close()
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