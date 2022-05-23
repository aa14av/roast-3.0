%% Perform leadField Estimation from ROAST
% By Alejandro Albizu on 6/5/2020
cd /ufrc/camctrp/working/Alejandro/toolboxes/roast-3.0/
addpath(genpath(pwd))
%%
clear

% Settings
%------------------------------------------------------
rootDir = '/ufrc/camctrp/working/aprinda/freesurfer_output';
subs = [9009 9015 9022 9040 9044 9045 9051 9021 9023 9031 9032 9047 9048 9054]; %
conductivities = struct('white',0.126,'gray',0.276,'csf',1.65,'bone',0.01,...
    'skin',0.465,'air',2.5e-14,'gel',0.3,'electrode',5.9e7); % Indahlastari2016
%------------------------------------------------------

% Create Electrode Pairs
fid = fopen('./elec72.loc'); C = textscan(fid,'%d %f %f %s');
fclose(fid); elecName = C{4}; for i=1:length(elecName), elecName{i} = strrep(elecName{i},'.',''); end

Nelec = length(elecName)-1; % Remove Reference (Iz)
elecIdx = zeros(2,Nelec*((Nelec-1)/2)); butt = 1;
elecPair = triu(repmat(1:Nelec,Nelec,1),1);
for e1 = 1:Nelec
    for e2 = 1:Nelec
        if elecPair(e1,e2) ~= 0
            elecIdx(1,butt) = e1;
            elecIdx(2,butt) = e2;
            butt = butt + 1;
        end
    end
end

logname = fullfile(rootDir,['RoastFSTarget_' datestr(now,30) '_log.txt']);
logfile(logname,'Beginning hiperoast-FStarget for SB ...');

for s = 1:length(subs)
    tic
    parfor e = 1:length(elecIdx)
        % Label current Iteration
        simTag = ['target-' elecName{elecIdx(1,e)} '-' elecName{elecIdx(2,e)}];
        subDir = fullfile(rootDir,['FS6.0_sub-' num2str(subs(s)) '_ses01'],'FStarget',simTag);
        if ~exist(subDir,'dir');mkdir(subDir);end
        
        % Copy Data
        if ~exist(fullfile(subDir,'T1.nii'),'file') || ...
                ~exist(fullfile(subDir,'T1_T1orT2_masks.mat'),'file') || ...
                ~exist(fullfile(subDir,'T1_T1orT2_seg8.mat'),'file')
            copyfile(fullfile(fileparts(subDir),'T1.nii'),fullfile(subDir,'T1.nii')); % T1
            copyfile(fullfile(fileparts(subDir),'T1_T1orT2_masks.nii'),fullfile(subDir,'T1_T1orT2_masks.nii')); % allMask
            copyfile(fullfile(fileparts(subDir),'T1_T1orT2_seg8.mat'),fullfile(subDir,'T1_T1orT2_seg8.mat')); % seg8
        end
        
        % Perform iterative ROAST
        subj = fullfile(subDir,'T1.nii');
        if ~exist(fullfile(subDir,['T1_' simTag '_result.mat']),'file')
            hipertarget(s,subj,elecName{elecIdx(1,e)},elecName{elecIdx(2,e)},simTag,logname);
        else
            logfile(logname,['sub-' num2str(subs(s)) ' COMPLETE !']);
        end
    end
    toc
end
disp('All MATLAB Workers Done !')