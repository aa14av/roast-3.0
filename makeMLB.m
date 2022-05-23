function matlabbatch = makeMLB(dirname,baseFilename)
% for d = 1:3

matlabbatch{1}.spm.spatial.normalise.write.subj.def = {fullfile(fileparts(dirname),'y_T1.nii')};
matlabbatch{1}.spm.spatial.normalise.write.subj.resample = {
                                                            fullfile(dirname,[baseFilename '.nii,1'])
                                                            fullfile(dirname,[baseFilename '.nii,2'])
                                                            fullfile(dirname,[baseFilename '.nii,3'])
                                                            };
matlabbatch{1}.spm.spatial.normalise.write.woptions.bb = [-90.5 -108.5 -90.5
                                                          89.5 107.5 89.5];
matlabbatch{1}.spm.spatial.normalise.write.woptions.vox = [1 1 1];
matlabbatch{1}.spm.spatial.normalise.write.woptions.interp = 4;
matlabbatch{1}.spm.spatial.normalise.write.woptions.prefix = 'wuf';