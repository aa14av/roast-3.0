%% Place Electrodes on Manual Segmentation
% Created by Alejandro Albizu on 02/02/2022
% Last Updated: 03/02/2023 by AA
cd /blue/camctrp/working/Alejandro/toolboxes/roast-3.0
% cd W:/camctrp/working/Alejandro/toolboxes/roast-3.0
addpath(genpath(pwd))
clear

rootDir = '/blue/camctrp/working/Alejandro/ACT';
% rootDir = 'W:/camctrp/working/Alejandro/ACT'; % HIPERGATOR
% rootDir = 'P:\ACT-head_models\FEM\manual_segmentation';
dims = [256 256 256];
mnames = flip({'wm','gm','eyes','csf','air','blood','cancellous','cortical','skin','fat','muscle'});
idx = flip(1:11); % HARDCODED
recipe = {'F3',-2,'F4',2};

% Mesh Options
rb = 5; % Radial Bound
ab = 30; % Angular Bound
db = 0.3; % Distance Bound
rr = 3; % Radius Edge Ratio
mv = 10; % Max Elem Volume

% adata = readtable('../../ACT/ACT_Data_w.csv'); 

% Run only CT + tDCS
% subs = adata{:,2};
% subs = adata{contains(adata{:,3},'Cognitive Training + tDCS') & adata.levels >= 720 & adata.stim >= 16,2};
% subnames = cellfun(@(x) ['FS6.0_sub-' num2str(x) '_ses-01_T1w'],
% num2cell(subs),'uni',0); % HIPERGATOR
% subnames = cellfun(@(x) ['FS6.0_sub-' num2str(x) '_ses-01_T1w'], num2cell(subs),'uni',0);
subfdr = dir(fullfile(rootDir,'FS6.0*'));
subnames = {subfdr.name}';

tic
missing = ones(length(subnames),1);
timing = zeros(length(subnames),1);
parfor s = 1:length(subnames)
    substart = tic;
    if exist(fullfile(rootDir,subnames{s}),'dir')
        subDir = fullfile(rootDir,subnames{s},'ROAST_11tis_Output_III');
%         subDir = fullfile(rootDir,subnames{s},'ROAST_11tis_Output');
        if ~exist(subDir,'dir'); mkdir(subDir); end

        if ~exist(fullfile(subDir,'T1_tDCSLAB_Jbrain.nii'),'file')
            c = load(fullfile(pwd,'conductivities','cond_11tis.mat'),'cond');
            c.cond.gel = 0.3; c.cond.electrode = 5.9e7;
            if length(c.cond.gel(:))==1
                c.cond.gel = repmat(c.cond.gel,1,2); % HARDCODED for 2 electrodes
            end
            if length(c.cond.electrode(:))==1
                c.cond.electrode = repmat(c.cond.electrode,1,2); % HARDCODED for 2 electrodes
            end
            try
                T1 = fullfile(subDir,'T1.nii');
                copyfile(fullfile(fileparts(subDir),'T1.nii'),T1)
                % PDRIVE
                if exist(fullfile(fileparts(subDir),'11tis','T1_T1orT2_masks.nii'),'file')
                    copyfile(fullfile(fileparts(subDir),'11tis','T1_T1orT2_masks.nii'),fullfile(subDir,'T1_T1orT2_masks.nii'))
                else
                    warning([subnames{s} ': 11-Tissue Missing !!']);
                    continue;
                end
                try
                    roast(s, T1 ,recipe, ...
                        'electype', {'pad','pad'}, ...
                        'elecsize', {[70 50 3],[70 50 3]}, ...
                        'elecOri', {'lr','lr'}, ...
                        'conductivities',c.cond, ...
                        'meshoptions',struct( ...
                        'radbound',rb, ...
                        'angbound',ab, ...
                        'distbound',db, ...
                        'reratio', rr, ...
                        'maxvol',mv), ...
                        'simulationTag', 'tDCSLAB');
                    missing(s) = 0; close all;
                    disp([subnames{s} ' Complete !'])
                    timing(s) = toc(substart);
                catch ME
                    delete(fullfile(subDir,'*')); % START OVER
                    warning(ME.message); close all;
                end
            catch ME
                warning(ME.message)
            end
        else
            missing(s) = 0;
        end
    end
end
toc