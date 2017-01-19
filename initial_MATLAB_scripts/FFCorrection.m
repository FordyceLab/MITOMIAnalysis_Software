function []=FFCorrection(var1,OVERRIDE) 
% FFCorrection is used for flat field correcting image sets from Micromanager 
% before attempting to stitch them identically using the Grid / Collection 
% stitching option in FIJI. Assumes stacked TIFF image output from uM.
% 
% Sequence of events:
% 1. Determine number of channels to be examined
% 2. Open directories containing source images
% 3. Extract metadata from all layers
% 4. Load image, correct it, and save in directory neighboring source images
% 5. Stitch one layer with FIJI plugin for MATLAB 
% 6. Apply that stitching pattern to remainder of image/channel sets
% 
% To use:
% 
% FFCorrection(NumChannels)
% 
% where NumChannels is the designated number of image sets to be stitched in
% sync with each other. It is important that all image sets are the same
% dimensions.
% 
% FFCorrection(NumChannels,OVERRIDE)
% 
% where OVERRIDE is a cell containing the FF channel names to be used in the
% correction. The Nikon microscopes have a tendency to get buggy and forget
% which channel they are using. 
% 
% Example: FFCorrection(3,{'Cy5','eGFP','NONE'})
% 
% Three channels will be FF corrected and stitched. The first image set will 
% be forced to be corrected with Cy5, the second will be forced to be
% corrected with eGFP and the thrid will be "corrected" with a blank image.
% 
% To setup:
% 
% FFCorrection('-setup')
% 
% Before FFCorrection can be used, all of the necessary paths and
% dependencies must be set and system requirements must be met. If MATLAB is
% 32-bit, it will not have sufficient access to RAM to stitch the images. 6Gb 
% of Java Heap Memory is recommended. The FIJI application in Mac cannot be 
% accessed directly through MATLAB, so update the path in the setup function 
% if FIJI is in a different location. The program has not been confirmed to
% work without the mij.jar file, so please add it for the time being.
%
% To remove:
% 
% FFCorrection('-remove')
% 
% When FFCorrection was setup, it modified paths in MATLAB to include the
% directories containing flat field correction images and the MIJI startup
% script in addition to changing the amount of RAM available to MIJI/ImageJ
% when it starts up. To undo all of these settings, run the command listed
% above.
% 
% Version 0.0.160526
%
% Robert Puccinelli
% robert.puccinelli@outlook.com
% Fordyce Lab, Stanford, 2016
%

%% Evaluate inputs

switch nargin
    
    case 0
        NumChannels=input('Number of image/channel sets to be processed: ');
        OVERRIDE=[];
        setPaths='false';
    
    case 1
        if ischar(var1)==1
            NumChannels=0;
            OVERRIDE=[];
            setPaths=var1;
            if strcmp(var1,'-setup')==1
                disp('FFCorrection will be setup in this run.')
            elseif strcmp(var1,'-remove')==1
                disp('FFCorrection paths and lock file will be removed.')
            else
                disp('Unrecognized input. Please try: help FFCorrection')
                return
            end
        else
            NumChannels=var1;
            OVERRIDE=[];
            setPaths='false';
        end
        
        
    case 2
        if ischar(var1)==1 || iscell(OVERRIDE)~=1
            disp('Incorrect input arguments. Please try: help FFCorrection')
            return
        else
            NumChannels=var1;
            setPaths='false';
        end
        
    otherwise
        disp('Too many input arguments. Please try: help FFCorrection')
        return
end

%% Evaluate setPath string

%Eval of setPaths
if strcmp(setPaths,'false')~=1
    
    %Optional setup call
    if strcmp(setPaths,'-setup')==1
        InitialSetupFunction(); 
        return
        
    %Optional uninstall call
    elseif strcmp(setPaths,'-remove')==1
        RemovalFunction();
        return
        
    %Someone might have accidentally changed a line in the input eval
    else
        disp('setPaths string is not recognized. Check input eval section')
        return
    end
end

%% Extract metadata from directories of interest

[M,ABORT]=MetadataExtraction(NumChannels,OVERRIDE);

%If aborted, save log and break program 
if ABORT==1
    cd(M(1).Directory)
    cd ..
    LOGFILENAME=['LOG_FFCorrection_' datestr(datetime('now')) '.mat'];
    save(LOGFILENAME,'M')
    disp('Log file saved.')
    return
end

%% Reformat images for flat field correction

[M,ABORT]=ReformatTIFF(NumChannels,M);

%If aborted, save log and break program 
if ABORT==1
    cd(M(1).Directory)
    cd ..
    LOGFILENAME=['LOG_FFCorrection_' datestr(datetime('now')) '.mat'];
    save(LOGFILENAME,'M')
    disp('Log file saved.')
    return
end

%% Run MIJI for automated image stitching

% [M,ABORT]=MIJIDeployment(M);

%If aborted, save log and break program
if ABORT==1
    cd(M(1).Directory)
    cd ..
    LOGFILENAME=['LOG_FFCorrection_Error' datestr(datetime('now')) '.mat'];
    save(LOGFILENAME,'M')
    disp('Log file saved.')
    return
end
%If successful, save log
cd(M(1).OutputPath)
LOGFILENAME=[ 'LOG_FFCorrection_Successful' datestr(datetime('now')) '.mat'];
save(LOGFILENAME, 'M')
disp('Log file saved.')
cd ..
disp('Congrats - FFCorrection was successfully completed.')

%% METADATA EXTRACTION FUNCTION
    function [M,ABORT]=MetadataExtraction(NumChannels,OVERRIDE)
        NumHold=1;

        %User loops over directories until number of image channels is satisfied
        while NumHold<=NumChannels

            %Open directory for new image set
            cd(uigetdir)
            cdir=pwd;
            M(NumHold).Directory=cdir;
            disp(M(NumHold).Directory)
            
            %Build file list and metadata
            imagelist=dir('*.ome.tif');
            NumFiles=length(imagelist);
            fInfo=imfinfo(imagelist(NumFiles,1).name);
            NumStacks=length(fInfo);
            fInfo=imfinfo(imagelist(1,1).name);
            longS=cell(NumStacks,200);
            disp(['NOTE: Number of stacks detected - ' num2str(NumStacks)])
            waitfor(1)

            %ImageDescription checks asks - Is this the source image file? 
            %If true, pull out relevant metadata for each layer
            if isfield(fInfo,'ImageDescription')==1	
                for a=1:NumStacks
                    meta=fInfo(a).UnknownTags;
                    ChanIdx=NumHold+a-1;
                    splitmeta=[];
                    M(ChanIdx).ImageSetIndex=ChanIdx;
                    M(ChanIdx).ImageLayerLocation=a;

                    if a==1 %&& NumStacks>1
                        splitmeta=strsplit(meta(3).Value,',');
                        longS(a,1:length(splitmeta))=splitmeta;
                    else
                        splitmeta=strsplit(meta(1).Value,',');
                        longS(a,1:length(splitmeta))=splitmeta;
                    end

                    %Find channel metadata for each stack	
                    FilterSearch=strncmp('"TIFilterBlock1-Label"',longS(a,:),22);
                    [rowF, colF]=find(FilterSearch>0);
                    FilterLocation=longS{rowF,colF};
                    FilterString=FilterLocation(25:end);

                    if strcmp('eGFP',FilterString(1:4))==1 
                        M(ChanIdx).Channel='eGFP';
                        disp('NOTE: eGFP detected.')
                    elseif strcmp('Cy5',FilterString(1:3))==1
                        M(ChanIdx).Channel='Cy5';
                        disp('NOTE: Cy5 detected.')
                    elseif strcmp('DAPI',FilterString(1:4))==1
                        M(ChanIdx).Channel='DAPI';
                        disp('NOTE: DAPI detected.')
                    else
                        disp('NOTE: Only eGFP, DAPI and Cy5 are supported for FF corrections.')
                        M(ChanIdx).Channel='NONE';
                    end
                    if isempty(OVERRIDE)==1
                        M(ChanIdx).Override=[];
                    else
                        M(ChanIdx).Override=OVERRIDE{ChanIdx};
                        disp(['NOTE: Channel is being overriden with ' OVERRIDE{ChanIdx} ' by user request.'])
                    end
                    waitfor(1)
                end

                %Find microscope metadata       
                ScopeSearch=strncmp('"TIScope-SoftwareVersion"',longS,25);
                [rowS,colS]=find(ScopeSearch>0);
                ScopeLocation=longS{rowS,colS};
                ScopeString=ScopeLocation(28:36);

%                 if strcmp(ScopeString,'4.4.1.728')==1
%                     disp('NOTE: Setup 1 used - Software Version 4.4.1.728')
%                     Scope='Setup2';
%                 elseif strcmp(ScopeString,'4.4.1.714')==1
%                     disp('NOTE: Setup 2 used - Software Version 4.4.1.714')
                    Scope='Setup2';
%                 else
%                     M(ChanIdx).Error=['Setup could not be determined. Software Version _' ScopeString '_ is not accounted for.'];
%                     disp(M(ChanIdx).Error);
%                     disp('ABORTING . . .')
%                     ABORT=1;
%                     return
%                 end        
%                 waitfor(1)

                %Find objective metadata        
%                 if strcmp(Scope,'Setup1')==1
%                 ObjectiveSearch=strncmp('"TIObjectiveTurret-Label"',longS,25);
%                 [rowO, colO]=find(ObjectiveSearch>0);
%                 ObjectiveLocation=longS{rowO,colO};
%                 ObjectiveString=ObjectiveLocation(28);
%                     if strcmp(ObjectiveString,'1')==1 
%                         Objective='4x';
%                         disp('NOTE: 4x Objective detected.')
%                     elseif strcmp(ObjectiveString,'2')==1
                        Objective='10x';
                        disp('NOTE: 10x Objective detected.')
%                     else
%                         M(ChanIdx).Error=['Objective _' ObjectiveString '_ is not accounted for.'];
%                         disp(M(ChanIdx).Error)
%                         disp('ABORTING . . .')
%                         ABORT=1;
%                         return
%                     end
%                 else
%                     Objective=input('Objective could not be determined for Setup2 - Objective? (4x/10x) :','s');
%                 end
                waitfor(1)

                %Find binning metadata        
                BinningSearch=strncmp('"Binning"',longS,9);
                [rowB,colB]=find(BinningSearch>0);
                BinningLocation=longS{rowB,colB};
                BinningString=BinningLocation(12:14);
                BinningNumber=str2double(BinningString(1,1));
                BinDisp=['NOTE: ',BinningString,' binning detected.'];
                if BinningNumber~=1 && BinningNumber~=2 && BinningNumber~=3
                    M(ChanIdx).Error=['Binning other than 1x1, 2x2 or 3x3 is not supported at this time. ' BinDisp];
                    disp(M(ChanIdx).Error)
                    disp('ABORTING . . .')
                    ABORT=1;
                    return
                end
                disp(BinDisp)
                waitfor(1)

                %Find ROI metadata
                %NOTE: This is for a custom ROI that was predefined
                ROISearch=strncmp('"ROI"',longS,5);
                [rowR,colR]=find(ROISearch>0);
                ROILocation=longS{rowR,colR};
                ROIString=ROILocation(8);
                if strcmp(ROIString,'0')==1 || strcmp(ROIString,'-')==1
                    ROI=0;
                    disp('NOTE: ROI was not detected.')
                else
                    ROI=1;
                    disp('NOTE: ROI was detected')
                end
                waitfor(1)

                %Pass metadata to logging structure
                for b=NumHold:ChanIdx
                    M(b).Scope=Scope;
                    M(b).Objective=Objective;
                    M(b).Binning=BinningNumber;
                    M(b).ROI=ROI;
                    M(b).Directory=cdir;
                    M(b).Error='None';
                end

                ABORT=0;       

            else

                %Program believes these are not the source files.
                M(NumHold).Error='Image metadata has been altered previously.';
                disp(M(NumHold).Error)
                disp('ABORTING . . .')
                ABORT=1;
                return
            end

        %Update counter
        NumHold=NumHold+NumStacks;
        cd ..
        end

        %Function completed
        disp('Metadata successfully extracted.')
        
    end

% %% REFORMAT TIFF FUNCTION
    function [M,ABORT]=ReformatTIFF(NumChannels,M)
        
        %Find and make directories for saving modified images
        disp('Select a directory where the "FFCorrectionOutput" directory can be made');
        waitfor(msgbox('Select a directory where the "FFCorrectionOutput" directory can be made'))
        cd(uigetdir);
        mkdir('FFCorrectionOutput')
        cd('FFCorrectionOutput')
        h=waitbar(0,'Flat Field Correction In Progress','CreateCancelBtn','setappdata(gcbf,''canceling'',1)','name','Flat Field Correction');
        
        for c=1:NumChannels
            tempDir=[];
            tempDirbBreak=[];
            tempDir=M(c).Directory;
            tempDirBreak=strsplit(tempDir,filesep);
            
            %Generate source file list and save dimensions
            fList=dir([M(c).Directory filesep '*.ome.tif']);
            M(c).NumFiles=length(fList);
            M(c).NumCol=length(dir([M(c).Directory filesep '*_000.ome.tif']));
            M(c).NumRow=floor(M(c).NumFiles/M(c).NumCol);
            
            %Load correction image and modify
            disp('Loading and adjusting correction image. . .')
            if ischar(M(c).Override)==1
                M(c).OutputPath=[pwd filesep num2str(M(c).ImageSetIndex) '_' M(c).Override '_' tempDirBreak{end}];
                M(c).CorrectionImageName=['CorrectionImage_' M(c).Scope '_' M(c).Objective '_' M(c).Override '.tif'];
                M(c).tag=M(c).Override;
            else
                M(c).OutputPath=[pwd filesep num2str(M(c).ImageSetIndex) '_' M(c).Channel '_' tempDirBreak{end}];
                M(c).CorrectionImageName=['CorrectionImage_' M(c).Scope '_' M(c).Objective '_' M(c).Channel '.tif']; %%%
                M(c).tag=M(c).Channel;
            end
            ValidateImage=which('-all',M(c).CorrectionImageName);
            if isempty(ValidateImage)==1
                M(c).Error=['Error locating: ' M(c).CorrectionImageName];
                disp(M(c).Error)
                disp('ABORTING . . .')
                ABORT=1;
                return
            end
            FFSource=imread(M(c).CorrectionImageName);
            [FFImage,medianFF]=CorrectionImageAdjustment(FFSource,M(c).Binning,M(c).ROI);
            disp('Correction image successfully loaded and adjusted.')
            mkdir(M(c).OutputPath)
            
            %Load source image, apply flat field, write output image
            for d=1:M(c).NumFiles
                waitbar((d+(c-1)*M(1).NumFiles)/(NumChannels*M(1).NumFiles),h);
                image=imread([M(c).Directory filesep fList(d,1).name],'Index',M(c).ImageLayerLocation);
                imfi=uint16(((double(image))*medianFF./(FFImage)));
                imwrite(imfi,[M(c).OutputPath filesep num2str(d) '.tif'],'TIFF','Compression','none');
            end           
        end
        ABORT=0;
        delete(h)
        pause(.0001)
    end

%% CORRECTION IMAGE ADJUSTMENT FUNCTION    
    function [FFImage,medianFF]=CorrectionImageAdjustment(FFSource, BinningNumber, ROI)
        %Run through binning number and ROI to determine propper image size
        switch BinningNumber
            
            case 1
                
                if ROI==1
                    FFImage=double(FFSource(257:1792,257:1792));
                else
                    FFImage=double(FFSource);
                end

            case 2
                
                if ROI==1
                    FFHolder=zeros(1024);
                    for i=1:1024
                        for j=1:1024
                            FFHolder(i,j)=round((double(FFSource((i-1)*2+2,(j-1)*2+2))+double(FFSource((i-1)*2+1,(j-1)*2+1))+double(FFSource((i-1)*2+1,(j-1)*2+2))+double(FFSource((i-1)*2+2,(j-1)*2+1)))/4);
                        end
                    end        
                    FFImage=FFHolder(129:896,129:896);
                else
                    FFImage=zeros(1024);
                    for i=1:1024
                        for j=1:1024
                            FFImage(i,j)=round((double(FFSource((i-1)*2+2,(j-1)*2+2))+double(FFSource((i-1)*2+1,(j-1)*2+1))+double(FFSource((i-1)*2+1,(j-1)*2+2))+double(FFSource((i-1)*2+2,(j-1)*2+1)))/4);
                        end
                    end
                end
                
            case 3
                
                if ROI==1
                
                    FFHolder=zeros(682);
                    for i=1:682
                        for j=1:682
                            FFHolder(i,j)=round((double(FFSource((i-1)*3+3,(j-1)*3+3))+double(FFSource((i-1)*3+1,(j-1)*3+1))+double(FFSource((i-1)*3+1,(j-1)*3+2))+double(FFSource((i-1)*3+2,(j-1)*3+1))+double(FFSource((i-1)*3+1,(j-1)*3+3))+double(FFSource((i-1)*3+3,(j-1)*3+1))+double(FFSource((i-1)*3+2,(j-1)*3+3))+double(FFSource((i-1)*3+3,(j-1)*3+2))+double(FFSource((i-1)*3+2,(j-1)*3+2)))/9);
                        end
                    end
                    FFImage=FFHolder(85:596,85:596);  
                else
                    FFImage=zeros(682);
                    for i=1:682
                        for j=1:682
                            FFImage(i,j)=round((double(FFSource((i-1)*3+3,(j-1)*3+3))+double(FFSource((i-1)*3+1,(j-1)*3+1))+double(FFSource((i-1)*3+1,(j-1)*3+2))+double(FFSource((i-1)*3+2,(j-1)*3+1))+double(FFSource((i-1)*3+1,(j-1)*3+3))+double(FFSource((i-1)*3+3,(j-1)*3+1))+double(FFSource((i-1)*3+2,(j-1)*3+3))+double(FFSource((i-1)*3+3,(j-1)*3+2))+double(FFSource((i-1)*3+2,(j-1)*3+2)))/9);
                        end
                    end           
                end
        end
        
        medianFF=median(FFImage(:)');
    end

%% MIJI DEPLOYMENT FUNCTION
    function [M,ABORT]=MIJIDeployment(M)

        %Determine if FIJI is part of path - if not, declare and return
        if ~exist('Miji.m')
            M(1).MIJIError='Miji.m was not detected in the path.';
            disp(M(1).MIJIError)
            ABORT=1;
            return
        end
            
        %Boot the MIJ and make a directory for stitched images
        MIJ=Miji('false');
        cd(M(1).OutputPath)
        cd ..
        mkdir([pwd filesep 'MIJI_Stitch'])
        cd('MIJI_Stitch')
        cdir=pwd;
        StitchOutput='img_t1_z1_c1';

        for e=1:length(M)
        
            if e==1
                %Customize string to run MIJI macro
                gridx=[' grid_size_x=' num2str(M(e).NumCol)];
                gridy=[' grid_size_y=' num2str(M(e).NumRow)];
                M(e).MIJIOutputDirectory=cdir;
                M(e).MIJIInput=['type=[Grid: column-by-column] order=[Down & Left]' gridx gridy ' tile_overlap=20 first_file_index_i=1 directory=[' M(e).OutputPath '] file_names={i}.tif output_textfile_name=TileConfiguration.txt fusion_method=[Linear Blending] regression_threshold=0.30 max/avg_displacement_threshold=2.50 absolute_displacement_threshold=3.50 compute_overlap invert_x computation_parameters=[Save computation time (but use more RAM)] image_output=[Write to disk] output_directory=[' M(e).MIJIOutputDirectory ']'];
                
                %Run MIJI macro
                MIJ.run('Grid/Collection stitching', M(e).MIJIInput)
                pause(.0001)

                %Rename output and duplicate TileConfiguration file
                M(e).FinalName=['Stitched_' num2str(M(e).ImageSetIndex) '_' M(e).tag '.tif'];
                movefile(StitchOutput,M(e).FinalName);
                pause(.0001)
                disp(['Stitch 1 out of ' num2str(length(M)) ' complete.'])
            else
                
                disp(['Copying coordinates from master stitch. Beginning stitch ' num2str(e) ' out of ' num2str(length(M)) '.'])
                copyfile([M(1).OutputPath filesep 'TileConfiguration.registered.txt'],[M(e).OutputPath filesep 'TileConfiguration.registered.txt']);
                %Customize string to run MIJI macro
                M(e).MIJIOutputDirectory=cdir;
                M(e).MIJIInput=['type=[Positions from file] order=[Defined by TileConfiguration] directory=[' M(e).OutputPath '] layout_file=TileConfiguration.registered.txt fusion_method=[Linear Blending] regression_threshold=0.30 max/avg_displacement_threshold=2.50 absolute_displacement_threshold=3.50 computation_parameters=[Save computation time (but use more RAM)] image_output=[Write to disk] output_directory=[' M(e).MIJIOutputDirectory ']'];
                
                %Run MIJI macro
                MIJ.run('Grid/Collection stitching', M(e).MIJIInput)
                pause(.0001)
                
                %Rename output
                M(e).FinalName=['Stitched_' num2str(M(e).ImageSetIndex) '_' M(e).tag '.tif'];
                movefile(StitchOutput,M(e).FinalName); 
                pause(.0001)
                disp(['Stitch ' num2str(e) ' out of ' num2str(length(M)) ' complete.'])
            end
        end
        %Close the MIJ
        MIJ.exit;
        ABORT=0;
    end

%% INTIAL SETUP FUNCTION
    function []=InitialSetupFunction()

        %Set FF path, FIJI path and Java  Memory for first system use
        

        %Need 64 bit version of Matlab for image stitching
        bitVer=computer('arch');
        if strcmp(bitVer(end-1:end),'64')==1
            disp('64-bit version of Matlab detected. Proceeding.')
        else
            disp('64-bit version of Matlab not detected. Aborting.')
            return
        end

        %Terminate setup process if it has been completed before        
        optsLocation=[matlabroot filesep 'bin' filesep bitVer filesep];
        LockFile=[optsLocation 'FFCorrectionLOCKED.mat'];
        if exist(LockFile)
            disp('Setup has been completed before. Aborting.')
            return
        end

        
        %Add FF folder to path
        disp('Navigate to the directory containing FF images.')
        waitfor(msgbox('Navigate to the directory containing FF images.'))
        FFpath=uigetdir();
        addpath(FFpath)
        savepath
        disp(['Path for FF images saved to system: ' FFpath])

        %Find FIJI and save Miji.m to path
        if ismac==1
            FIJIpath='/Applications/Fiji.app/scripts';
            addpath(FIJIpath)
        else
            disp('Navigate to Fiji/scripts folder.')
            waitfor(msgbox('Navigate to Fiji/scripts folder.'))
            FIJIpath=uigetdir();
            addpath(FIJIpath)
        end
        savepath
        disp('Miji.m saved to path in MATLAB.')
%         disp('Please update FIJI. Loading MIJI. . .')
%         pause(2)
%         MIJ=Miji('false');
%         MIJ.run('Update Fiji');
%         MIJ.exit();

        %Add 6Gb RAM to ImageJ Java if Java setting is not present.
        cdir=pwd;
        optsFile=[optsLocation 'java.opts'];
        if ~exist(optsFile)
            cd(optsLocation)
            fid=fopen(optsFile);
            fSpec='%s\n';
            fprintf(fid,fSpec,'-Xmx6037m');
            fclose(fid);
            disp('Java max heaps set to 6Gb.')
        else
            cd(optsLocation)
            disp('Java file already exists.')
            optsFileD=importdata('java.opts');
            XmxLocation=find(strncmp(optsFileD,'-Xmx6037m',9));
            if sum(XmxLocation)==0
                optsFileD{end+1}='-Xmx6037m';
                fid=fopen('java.opts','w');
                nrow=length(optsFileD);
                fSpec='%s\n';
                for f=1:nrow
                    fprintf(fid,fSpec,optsFileD{f});
                end
                fclose(fid);
                disp('Java max heaps setting inserted as 6Gb.')
            elseif sum(XmxLocation)==1
                optsFileD{XmxLocation}='-Xmx6037m';
                fid=fopen('java.opts','w');
                nrow=length(optsFileD);
                fSpec='%s\n';
                for g=1:nrow
                    fprintf(fid,fSpec,optsFileD{g});
                end
                fclose(fid);
                disp('Java max heaps setting changed to 6Gb.')
            else
                disp('Unexpected encounter. Please manually change -Xmx setting to -Xmx6037Gb')
                edit(java.opts)
            end
        end
        save('FFCorrectionLOCKED.mat','FFpath','FIJIpath')
        cd(cdir)
        disp('Installation complete.')
    end
    
%% REMOVAL FUNCTION   
    function []=RemovalFunction()
        
        %Terminate removal process if it has never been completed before 
        bitVer=computer('arch');
        optsLocation=[matlabroot filesep 'bin' filesep bitVer filesep];
        LockFile=[optsLocation 'FFCorrectionLOCKED.mat'];
        if ~exist(LockFile)
            disp('Setup has never been completed before. Aborting.')
            return
        end
        
        %Remove paths
        VAR=load(LockFile);
        rmpath(VAR.FFpath)
        disp('Flat field correction image directory removed from path.')
        rmpath(VAR.FIJIpath)
        disp('FIJI/MIJI removed from path.')
        
        %Remove max memory setting and save file
        cdir=pwd;
        cd(optsLocation)
        optsFile=importdata([optsLocation 'java.opts']);
        XmxLocation=find(strncmp(optsFile,'-Xmx6037m',9));
        optsFile(XmxLocation)=[];
        nrows=length(optsFile);
        copyfile('java.opts','java_removed.opts');
        fid=fopen('java.opts','w');
        fSpec='%s\n';
        for h=1:nrows
            fprintf(fid,fSpec,optsFile{h});
        end
        fclose(fid);
        disp('Max Java Heaps memory setting has been removed.')
        delete(LockFile);
        disp('Lock file deleted.')
        cd(cdir)
        disp('Uninstallation complete.')
        savepath
        disp('Path changes saved.')
        return
        
    end

end

%% NOTEPAD
%{

%% Possible upgrades
    %Give the user an option to select which layer is the master stitch
    %Reject input directory if it doesn't meet stacked tiff name structure
    %Or include function that rewrites file name to stacked tiff structure
        %Is the metadata still the same format?
    %Sample function for converting:
%         folderlist=dir('*-*');
%         curdir=pwd;
%         for i=1:length(folderlist)
%             cd(folderlist(i).name);
%             filelist=dir('*tif');
%             folderdir=pwd;
%             movefile([filelist(1).name],['../' folderlist(i).name '.ome.tif']);
%             cd ..
%             folderdir=[];
%             filelist=[];
%         end
    


%% Formal MIJ macros                

    %First image set setting stitching standards (Protein channel most stable, DNA chamber often works too)
    MIJ.run('Grid/Collection stitching', 'type=[Grid: column-by-column] order=[Down & Left] grid_size_x=18 grid_size_y=18 tile_overlap=20 first_file_index_i=1 directory=[/Volumes/NO NAME/20160518_Pho4_RepetitivePlateWithCGC_Postwash_eGFP_3000ms_3x3Binning/renamed] file_names={i}.tif output_textfile_name=TileConfiguration.txt fusion_method=[Linear Blending] regression_threshold=0.30 max/avg_displacement_threshold=2.50 absolute_displacement_threshold=3.50 compute_overlap invert_x computation_parameters=[Save computation time (but use more RAM)] image_output=[Write to disk] output_directory=[/Volumes/NO NAME/20160518_Pho4_RepetitivePlateWithCGC_Postwash_eGFP_3000ms_3x3Binning/renamed]]');

    %Remaining image sets stitched identically to first
    MIJ.run('Grid/Collection stitching', 'type=[Positions from file] order=[Defined by TileConfiguration] directory=[/Volumes/NO NAME/20160518_Pho4_RepetitivePlateWithCGC_Postwash_eGFP_3000ms_3x3Binning/renamed] layout_file=TileConfiguration.registered.txt fusion_method=[Linear Blending] regression_threshold=0.30 max/avg_displacement_threshold=2.50 absolute_displacement_threshold=3.50 computation_parameters=[Save computation time (but use more RAM)] image_output=[Write to disk] output_directory=[/Volumes/NO NAME/20160518_Pho4_RepetitivePlateWithCGC_Postwash_eGFP_3000ms_3x3Binning/renamed]');


%% Gibberish to help with coding

    % "TIFilterBlock1-Label":"4-Cy5"
    % "TIObjectiveTurret-Label":"2-Plan Apo 10x NA 0.45 Dry"
    % "Binning":"3x3"
    % "ROI":"0-0-2048-2048" {1x1}
    % "ROI":"-1--1-1024-1024" {2x2}
    % "ROI":"-1--1-682-682" {3x3}
    % "ROI":"257-257-1536-1536" {custom roi 1x1}
    % "ROI":"128-128-768-768" {custom roi 2x2}
    % "ROI":"85-85-512-512" {custom roi 3x3}
    % "Andor sCMOS Camera-Exposure":"34.0000" {ms}
    % "ZPositionUm": 5222.150077816099 {um}

%}