%% Run ROAST
%----------------------------------------
% Created By Alejandro Albizu
% Center for Cognitive Aging and Memory
% University of Florida
% 06/31/2021
%----------------------------------------
% Last Updated: 06/31/2021 by AA
parfor s = 1:length(subnames)
    while ~exist(fullfile(subDir,subnames{s},'ROAST_output','T1_Jbrain.nii'),'file')
        subj = fullfile(subDir,subnames{s},'ROAST_output','T1.nii');
        if ~exist(fullfile(subDir,subnames{s},'ROAST_output'),'dir')
            mkdir(fullfile(subDir,subnames{s},'ROAST_output'))
            copyfile(fullfile(subDir,subnames{s},'T1.nii'),fullfile(subDir,subnames{s},'ROAST_output','T1.nii'))
            copyfile(fullfile(subDir,subnames{s},'c1T1_T1orT2.nii'),fullfile(subDir,subnames{s},'ROAST_output','c1T1_T1orT2.nii'))
            copyfile(fullfile(subDir,subnames{s},'T1_T1orT2_masks.nii'),fullfile(subDir,subnames{s},'ROAST_output','T1_T1orT2_masks.nii'))
            copyfile(fullfile(subDir,subnames{s},'T1_T1orT2_rmask.mat'),fullfile(subDir,subnames{s},'ROAST_output','T1_T1orT2_rmask.mat'))
            copyfile(fullfile(subDir,subnames{s},'T1_T1orT2_seg8.mat'),fullfile(subDir,subnames{s},'ROAST_output','T1_T1orT2_seg8.mat'))
            copyfile(fullfile(subDir,subnames{s},'T1_customLocations.txt'),fullfile(subDir,subnames{s},'ROAST_output','T1_customLocations'))
        end
        try
            roast(s,subj,{'F3',-2,'F4',2}, ...
                'electype', {'pad','pad'}, ...
                'elecsize', {[70 50 3],[70 50 3]}, ...
                'elecori', {'lr','lr'}, ...
                'simulationTag', 'MVPA');
        catch
           disp(['RESTARTING ' subnames{s} '...'])
           rmdir(fullfile(subDir,subnames{s},'ROAST_output'),'s');
        end
    end
    disp([subnames{s} ' COMPLETE !'])
end