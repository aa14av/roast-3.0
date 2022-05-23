function matlabbatch = makeMLB_UFAB(subDir,baseFilename)        
matlabbatch{1}.spm.spatial.normalise.write.subj.def = {fullfile(subDir,'y_T1.nii')};
matlabbatch{1}.spm.spatial.normalise.write.subj.resample = {fullfile(subDir,[baseFilename '.nii,1'])
    fullfile(subDir,[baseFilename '.nii,2'])
    fullfile(subDir,[baseFilename '.nii,3'])};
matlabbatch{1}.spm.spatial.normalise.write.woptions.bb = [-90.5 -108.5 -90.5
    89.5 107.5 89.5];
matlabbatch{1}.spm.spatial.normalise.write.woptions.vox = [1 1 1];
matlabbatch{1}.spm.spatial.normalise.write.woptions.interp = 4;
matlabbatch{1}.spm.spatial.normalise.write.woptions.prefix = 'wuf';
