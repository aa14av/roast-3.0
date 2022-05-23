%% Perform leadField Estimation from ROAST
% Created by Alejandro Albizu on 6/5/2020
% Updated by Alejandro Albizu on 08/15/2020
% 20200810: Changed Disk electrodes to pads via seperate mesh
cd /blue/camctrp/working/Alejandro/toolboxes/roast-3.0/
addpath(genpath(pwd))
%%
clear

% Settings
%------------------------------------------------------
rootDir = '/blue/camctrp/working/aprinda/freesurfer_output';
subs = [9054 1892 101549 300100 300700 301051]; %
simTag = 'FStarget';
dims = [256 256 256];
%------------------------------------------------------

fid = fopen('./elec72.loc'); C = textscan(fid,'%d %f %f %s'); fclose(fid);
elecName = C{4}; for i=1:length(elecName), elecName{i} = strrep(elecName{i},'.',''); end

tic
for s = 1:length(subs)
    sub = subs(s);
    subj = fullfile(rootDir,['FS6.0_sub-' num2str(sub) '_ses01'],simTag,'T1.nii');
    [dirname,baseFilename] = fileparts(subj);
    if exist(dirname,'dir');  rmdir(dirname,'s'); end; mkdir(dirname); % Start Fresh
    copyfile(fullfile(fileparts(dirname),'T1.nii'),fullfile(dirname,'T1.nii'))
    copyfile(fullfile(fileparts(dirname),'T1_T1orT2_masks.nii'),fullfile(dirname,'T1_T1orT2_masks.nii'))
    copyfile(fullfile(fileparts(dirname),'T1_T1orT2_seg8.mat'),fullfile(dirname,'T1_T1orT2_seg8.mat'))
    parfor e = 1:71 % HARDCODED (cuz im lazy)
        simTag = elecName{e};
        if ~exist(fullfile(dirname,[baseFilename '_' simTag '_roastResult.mat']),'file')
    %         try
        
            hiperskipseg(e,[],subj,'elec',elecName([e,end])')
        
%         catch
%             warning(['sub-' num2str(sub) ': ' elecName{e} ' FAILED ...'])
%         end
        else
            disp(['sub-' num2str(sub) ': ' elecName{e} ' SKIPPED ...'])
        end
    end
    A_all = zeros(prod(dims)*3,71); inc = zeros(71,1);
    for e = 1:71
        simTag = elecName{e};
        if exist(fullfile(dirname,[baseFilename '_' simTag '_roastResult.mat']),'file')
            load(fullfile(dirname,[baseFilename '_' simTag '_roastResult.mat']),'ef_all')
            A_all(:,e) = ef_all(:);
            inc(e) = 1;
        end
    end
    save(fullfile(dirname,[baseFilename '_leadField.mat']),'A_all','-v7.3')
end
toc
disp('All MATLAB Workers Done !')