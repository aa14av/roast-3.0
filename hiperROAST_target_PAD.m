%% Perform leadField Estimation from ROAST
% Created by Alejandro Albizu on 06/05/2020
% Updated by Alejandro Albizu on 07/09/2021
% 20200810: Changed Disk electrodes to pads via seperate mesh
cd /blue/camctrp/working/Alejandro/toolboxes/roast-3.0/
addpath(genpath(pwd)); addpath ../NIFTI_20130306/
%%
clear

% Settings
%------------------------------------------------------
rootDir = '/blue/camctrp/working/Alejandro/StimBrain';
simTag = 'FStarget';
%------------------------------------------------------

subs = [9009 9015 9022 9040 9044 9045 9051 9021 9023 9031 9032 9047 9048 9054]';
%subs = [9048]'; %9051 9021 9023 9031 1892 101549 300100 300700 301051]; % 9054 300700
subnames = cellfun(@(x) ['FS6.0_sub-' num2str(x,['%0' num2str(length(num2str(subs(1)))) '.f']) '_ses01'],mat2cell(subs,ones(length(subs),1)),'uni',0);
% subfdr = dir(fullfile(rootDir,'FS6.0_sub-*'));
% subnames = {subfdr.name}';

fid = fopen('./elec72.loc'); C = textscan(fid,'%d %f %f %s'); fclose(fid);
elecName = C{4}; for i=1:length(elecName), elecName{i} = strrep(elecName{i},'.',''); end

tic
for s = 1:length(subnames)
    subname = subnames{s}; if ~exist(fullfile(rootDir,subname,simTag),'dir'); mkdir(fullfile(rootDir,subname,simTag)); end
    subj = fullfile(rootDir,subname,simTag,'T1.nii'); [dirname,baseFilename] = fileparts(subj);
%     if exist(dirname,'dir');  rmdir(dirname,'s'); end; mkdir(dirname); % Start Fresh
    copyfile(fullfile(fileparts(dirname),'T1.nii'),fullfile(dirname,'T1.nii'))
    copyfile(fullfile(fileparts(dirname),'T1_T1orT2_masks.nii'),fullfile(dirname,'T1_T1orT2_masks.nii'))
    copyfile(fullfile(fileparts(dirname),'T1_T1orT2_seg8.mat'),fullfile(dirname,'T1_T1orT2_seg8.mat'))
    parfor e = 1:71 % HARDCODED (cuz im lazy)
        eTag = elecName{e}; m=1;
        while ~exist(fullfile(dirname,[baseFilename '_' eTag '_roastResult.mat']),'file')
            if m > 2; disp([subname ': ITERATION LIMIT REACHED FOR ' elecName{e}]); break; end % QUIT AFTER 5 ITERATIONS
            try; hiperskipseg(e,[],subj,'elec',elecName([e,end])'); close all; catch ERROR; warning([subname ' (' eTag '): ' ERROR.identifier]); end
            if ~exist(fullfile(dirname,[baseFilename '_' eTag '_roastResult.mat']),'file') % START FRESH
                disp(['RESTARTING ' subname ': ' elecName{e}]); m = m + 1;
                delete(fullfile(dirname,[baseFilename '_' eTag '_*.*']));
            end
        end
    end
end
toc; disp('All MATLAB Workers Done !');