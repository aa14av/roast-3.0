%% Perform leadField Estimation from ROAST
% Created by Alejandro Albizu on 6/5/2020
% Updated by Alejandro Albizu on 08/10/2020
% 20200810: Changed Disk electrodes to pads via seperate mesh
cd /blue/camctrp/working/Alejandro/toolboxes/roast-3.0/
addpath(genpath(pwd))
%%
clear

% Settings
%------------------------------------------------------
rootDir = '/blue/camctrp/working/aprinda/freesurfer_output';
subs = [9009 9015 9022 9040 9044 9045 9051 9021 9023 9031 9032 9047 9048 9054]; %
simTag = 'FStarget';
% conductivities = struct('white',0.3835,'gray',0.1,'csf',1.8,'bone',0.0109,...
%     'skin',0.43,'air',2.5e-14,'gel',0.3,'electrode',5.9e7); % Indahlastari2016
conductivities = struct('white',0.126,'gray',0.276,'csf',1.65,'bone',0.01,...
    'skin',0.465,'air',2.5e-14,'gel',0.3,'electrode',5.9e7);
%------------------------------------------------------

fid = fopen('./elec72.loc'); C = textscan(fid,'%d %f %f %s'); fclose(fid);
elecName = C{4}; for e=1:length(elecName), elecName{e} = strrep(elecName{e},'.',''); end
if length(conductivities.gel(:))==1
    conductivities.gel = repmat(conductivities.gel,1,length(elecName));
end
if length(conductivities.electrode(:))==1
    conductivities.electrode = repmat(conductivities.electrode,1,length(elecName));
end
indInUsrInput = [1;66;2;48;30;67;31;49;56;11;38;3;20;17;21;4;39;12;57;58;50;40;32;22;68;23;33;41;51;59;13;42;5;24;18;25;6;43;14;52;44;34;26;69;27;35;45;53;60;15;46;7;28;19;29;8;47;16;61;62;54;36;70;37;55;63;9;71;10;64;65;72];
elecName = elecName(indInUsrInput);
conductivities.gel = conductivities.gel(indInUsrInput);
conductivities.electrode = conductivities.electrode(indInUsrInput);

capType = '1010';
elecType = 'pad';
elecSize = [70 50 3];
elecOri = 'lr';

elecPar = struct('capType',capType,'elecType',elecType,...
    'elecSize',elecSize,'elecOri',elecOri);
configTxt = 'leadFieldGeneration';
meshOpt = struct('radbound',5,'angbound',30,'distbound',0.3,'reratio',3,'maxvol',10);

logname = fullfile(rootDir,['RoastFSTarget_' datestr(now,30) '_log.txt']);
logfile(logname,'Beginning hiperoast-FStarget for SB');

tic
for s = 1:length(subs)
    subDir = fullfile(rootDir,['FS6.0_sub-' num2str(subs(s)) '_ses01'],simTag);
    if ~exist(subDir,'dir');mkdir(subDir);end
    if ~exist(fullfile(subDir,'T1.nii'),'file') || ...
            ~exist(fullfile(subDir,'T1_T1orT2_masks.mat'),'file') || ...
            ~exist(fullfile(subDir,['T1_' simTag.mat']),'file') || ...
            ~exist(fullfile(subDir,['T1_' simTag '_header.mat']),'file')
        copyfile(fullfile(fileparts(subDir),'T1.nii'),fullfile(subDir,'T1.nii')); % T1
        copyfile(fullfile(fileparts(subDir),'T1_T1orT2_masks.nii'),fullfile(subDir,'T1_T1orT2_masks.nii')); % allMask
        copyfile(fullfile(fileparts(subDir),'T1_T1orT2_seg8.mat'),fullfile(subDir,'T1_T1orT2_seg8.mat')); % seg8
    end
    subj = fullfile(subDir,'T1.nii');
    [dirname,baseFilename] = fileparts(subj);
    if ~exist(fullfile(subDir,'T1_target_roastResult.mat'),'file')
        [~,indRef] = ismember('Iz',elecName);
        indStimElec = setdiff(1:length(elecName),indRef);
        [isInRoastCore,indInRoastCore] = ismember(elecName,elecName(indStimElec));
        isSolved = zeros(length(indStimElec),1);
        parfor e=1:length(indStimElec)
            if exist([dirname filesep baseFilename '_' simTag '_e' num2str(indStimElec(e)) '.pos'],'file')
                isSolved(e) = 1;
            end
        
            % ELECTRODE PARAMETERS
            try
                [elecPara,~] = elecPreproc(subj,elecName([e end]),elecPar);
                options = struct('configTxt',configTxt,'elecPara',elecPara,'T2',[],'meshOpt',meshOpt,'conductivities',conductivities,'uniqueTag',simTag,'resamp',0,'zeroPad',0,'isNonRAS',0);
            catch
                logfile(logname,['sub-' num2str(subs(s)) ' FAILED at PRE-PLACE...']);
                continue
            end
            
            % ELECTRODE PLACEMENT
            %         if ~exist(fullfile(dirname,[baseFilename '_header.mat']),'file')
            try
                hdrInfo = electrodePlacement(subj,subj,[],elecName([e end]),elec(e).options,simTag);
            catch
                logfile(logname,['sub-' num2str(subs(s)) ' FAILED at PLACE...']);
                continue
            end
            %         else
            %             hdr = load(fullfile(dirname,[baseFilename '_header.mat']),'hdrInfo');
            %             hdrInfo = hdr.hdrInfo;
            %         end
            
            % MESH
            %         if ~exist(fullfile(dirname,[baseFilename '_' simTag '.mat']),'file')
            try
                [node,elem,~] = meshByIso2mesh(subj,subj,[],meshOpt,hdrInfo,simTag);
            catch
                logfile(logname,['sub-' num2str(subs(s)) ' FAILED at MESH...']);
                continue
            end
            %         else
            %             mesh = load(fullfile(dirname,[baseFilename '_' simTag '.mat']),'node','elem');
            %             node = mesh.node;
            %             elem = mesh.elem;
            %         end
            
            if ~exist(fullfile(dirname,[baseFilename '_' simTag '_v.pos']),'file')
                % PREPARE FOR FE SOLVER
                try
                    prepareForGetDP(subj,node,elem,elecName,simTag);
                catch
                    logfile(logname,['sub-' num2str(subs(s)) ' FAILED at PREP...']);
                    continue
                end
                
                % FE SOLVER
                try
                    injectCurrent = ones(length(elecName),1); % 1 mA at each candidate electrode
                    injectCurrent(indRef) = -1;
%                     for i=1:length(indStimElec)
%                         if ~isSolved(i)
                            fprintf('\n======================================================\n');
                            disp(['SOLVING FOR ELECTRODE ' num2str(e) ' OUT OF ' num2str(length(indStimElec)) ' ...']);
                            fprintf('======================================================\n\n');
                            tmpDir = fullfile(dirname,['tmp' num2str(e)]);
                            if ~exist(tmpDir,'dir'); mkdir(tmpDir); end
                            indElecSolve = [indStimElec(e) indRef];
                            copyfile(fullfile(dirname,[baseFilename '_' simTag '_ready.msh']),fullfile(tmpDir,[baseFilename '_' simTag '_ready.msh']))
                            par_solveByGetDP(tmpDir,subj,injectCurrent,conductivities,indElecSolve,simTag,num2str(indStimElec(e)));
%                         else
%                             disp(['ELECTRODE ' num2str(i) ' HAS BEEN SOLVED, SKIPPING...']);
%                         end
%                     end
                catch
                    logfile(logname,['sub-' num2str(subs(s)) ' FAILED at SOLVE...']);
                    continue
                end
            end
            
            try
                hdr=load(fullfile(subDir,'T1_header.mat'),'hdrInfo');
                postGetDP(subj,[],mesh.node,hdr.hdrInfo,simTag,indStimElec,indInRoastCore(isInRoastCore));
                logfile(logname,['sub-' num2str(subs(s)) ' COMPLETE !']);
            catch
                logfile(logname,['sub-' num2str(subs(s)) ' FAILED at POST...']);
                continue
            end
        end
        close all
    else
        logfile(logname,['sub-' num2str(subs(s)) ' COMPLETE !']);
    end
end
toc
disp('All MATLAB Workers Done !')