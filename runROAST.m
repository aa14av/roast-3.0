%% Run ROAST Manual Segmentation
%====================================================
% Created by: Alejandro Albizu for the Center for Cognitive Aging and Memory
% University of Florida
% Email: aa14av@gmail.com
% Created: 02/02/2022
% Updated: 08/08/2022
%====================================================
clear % Clear Workspace

% Settings
%--------------------------------
rootDir = 'directory containing individual subject folder(s)';
recipe = {}; % ROAST recipe (e.g. {'F3',-2,'F4',2})
elecType = {}; % ROAST elecTypes (e.g. {'pad','pad'})
elecSize = {}; % ROAST elecSizes (e.g. {[70 50 3],[70 50 3]})
elecOri = {}; % ROAST elecOris (e.g. {'lr','lr'})
simTag = 'unique tag for ROAST output file(s)';
%--------------------------------

% Locate Subject Folders
subfdr = dir(fullfile(rootDir,'sub*'));
subnames = {subfdr.name}';

tic
missing = ones(length(subnames),1); % Pre-allocate
for s = 1:length(subnames)
    subDir = fullfile(rootDir,subnames{s},'ROAST_Output');
    T1 = 'absolute path to T1 nifti file';
    condFile = 'absolute path to conductivity mat file';
    segFile = 'absolute path to segmentation nifti file';
    if ~exist(subDir,'dir'); mkdir(subDir); end % Create Output Folder 
    
    if ~exist(fullfile(subDir,[baseFilename '_' simTag '_Jbrain.nii']),'file') % Check if ROAST is already completed
        c = load(condFilename),'cond');
        cond = cell2struct(c.cond(:,4),c.cond(:,3)); % convert to struct for ROAST
        cond.gel = 1;                                                                   % HARDCODED GEL CONDUCTIVITY
        cond.electrode = 2.5e7;                                                         % HARDCODED ELEC CONDUCTIVITY
        cond.index = cell2mat(c.cond(:,1)); % Get Unique Tissue Indexes
        cond.brain = cell2mat(c.cond(:,2)); % Boolean Index of Brain vs Non-brain

        % Repeat Gel Conductivity for each electrode
        if length(cond.gel(:))==1
            cond.gel = repmat(cond.gel,1,length(elecSize));                                               
        end
        
        % Repeat Electrode Conductivity for each Electrode
        if length(cond.electrode(:))==1
            cond.electrode = repmat(cond.electrode,1,length(elecSize));                               
        end
        
        % Copy T1 to ROAST Output Directory
        if exist(T1,'file') &&...
            ~exist(fullfile(subDir,'T1.nii'),'file')
            copyfile(T1,fullfile(subDir,'T1.nii')); 
        end 

        % Copy Segmentation to ROAST Output Directory
        if exist(segFile,'file') &&...
            ~exist(fullfile(subDir,'T1_T1orT2_masks.nii'),'file')
            copyfile(segFile,fullfile(subDir,'T1_T1orT2_masks.nii')); 
        end
        
        % Run ROAST with specified settings (no need for T2 with manual seg)
        try
            roast(fullfile(subDir,'T1.nii') ,recipe, ...
                'electype', elecType, ...
                'elecsize', elecSize, ...
                'elecOri', elecOri, ...
                'conductivities',cond, ...
                'T2', [], 'simulationTag', simTag);
             missing(s) = 0; % ROAST Complete
             disp([subnames{s} ' Complete !']); % lmk when finished
             close all; % Close ROAST figures
        catch ME
            delete(fullfile(subDir,'*')); % START OVER
            warning(ME.message); % Print ROAST fail error
        end
        
    else
        missing(s) = 0; % ROAST already complete 
    end
    
end
toc
