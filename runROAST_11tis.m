%% Place Electrodes on Manual Segmentation
% Created by Alejandro Albizu on 02/02/2022
% Last Updated: 06/07/2022 by AA
clear

% rootDir = '/Volumes/woodslab/ACT-head_models/FEM/manual_segmentation';
% condDir = '/Volumes/camctrp/working/Alejandro/toolboxes/tDCSLAB-2.0/';
%aprinda
% rootDir = '/Volumes/projects/woodslab/ACT-head_models/FEM/manual_segmentation';
% condDir = '/Volumes/working/Alejandro/toolboxes/tDCSLAB-2.0/';

rootDir = 'P:\WoodsLab\ACT-head_models\FEM\manual_segmentation\Education_Training';
condDir = 'P:\WoodsLab\ACT-head_models\tDCSLAB-3.0';
hrDir = 'W:\camctrp\working\Alejandro\ACT\headreco';
dims = [256 256 256];
mnames = flip({'wm','gm','eyes','csf','air','blood','cancellous','cortical','skin','fat','muscle'});
idx = flip(1:11); % HARDCODED
recipe = {'F3',-2,'F4',2};

subfdr = dir(fullfile(rootDir,'FS*'));
subnames = {subfdr.name}';

tic
missing = ones(length(subnames),1);
% for s = 140:length(subnames)
for s = 1:length(subnames)
    subDir = fullfile(rootDir,subnames{s},'ROAST_11tis_Output_II');
    if ~exist(subDir,'dir'); mkdir(subDir); end
    
    if ~exist(fullfile(subDir,'T1_tDCSLAB_Jbrain.nii'),'file')
        c = load(fullfile(condDir,'cond_11tis.mat'),'cond');
        cond = cell2struct(c.cond(:,4),c.cond(:,3)); cond.gel = 1; cond.electrode = 2.5e7;
        cond.index = cell2mat(c.cond(:,1)); cond.brain = cell2mat(c.cond(:,2));
        if length(cond.gel(:))==1
            cond.gel = repmat(cond.gel,1,2); % HARDCODED for 2 electrodes        
        end
        if length(cond.electrode(:))==1
            cond.electrode = repmat(cond.electrode,1,2); % HARDCODED for 2 electrodes
        end
        try
            T1 = fullfile(subDir,'T1.nii');
            copyfile(fullfile(fileparts(subDir),'T1.nii'),T1)
            if ~exist(fullfile(fileparts(subDir),'T1_T1orT2_masks.nii'),'file')
                copyfile(fullfile(fileparts(subDir),'T1_T1orT2_masks.nii'),fullfile(subDir,'T1_T1orT2_masks.nii'))
            else
            end
%             copyfile(fullfile(fileparts(subDir),'ROAST_11tis_Output','T1_T1orT2_seg8.mat'),fullfile(subDir,'T1_T1orT2_seg8.mat'))
%             copyfile(fullfile(fileparts(subDir),'ROAST_11tis_Output','T1_header.mat'),fullfile(subDir,'T1_header.mat'))
%             copyfile(fullfile(fileparts(subDir),'ROAST_11tis_Output','T1_tDCSLAB_mask_elec.nii'),fullfile(subDir,'T1_tDCSLAB_mask_elec.nii'))
%             copyfile(fullfile(fileparts(subDir),'ROAST_11tis_Output','T1_tDCSLAB_mask_gel.nii'),fullfile(subDir,'T1_tDCSLAB_mask_gel.nii'))
%             if ~exist(fullfile(subDir,'c1T1_T1orT2.nii'),'file'); fid = fopen(fullfile(subDir,'c1T1_T1orT2.nii'),'w'); fprintf(fid,'temp'); end
            try
                roast(T1 ,recipe, ...
                    'electype', {'pad','pad'}, ...
                    'elecsize', {[70 50 3],[70 50 3]}, ...
                    'elecOri', {'lr','lr'}, ...
                    'conductivities',cond, ...
                    'T2', [], 'simulationTag', 'tDCSLAB'); 
                 missing(s) = 0; disp([subnames{s} ' Complete !'])
                 close all;
%                 if exist(fullfile(subDir,'c1T1_T1orT2.nii'),'file'); delete(fullfile(subDir,'c1T1_T1orT2.nii')); end
            catch ME
                delete(fullfile(subDir,'*')); % START OVER
                warning(ME.message);
            end
        catch ME
            warning(ME.message)
        end
    else
        missing(s) = 0;
    end
end
toc