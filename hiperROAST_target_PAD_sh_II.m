%% Perform leadField Estimation from ROAST
% Created by Alejandro Albizu on 06/05/2020
% Updated by Alejandro Albizu on 07/30/2021
% 20200810: Changed Disk electrodes to pads via seperate mesh
% 20210730: Converted to Shell Script for HiPerGator

function hiperROAST_target_PAD_sh_II(rootDir,subname,idx,lim)
if ~isa(idx,'double'); idx = str2double(idx); end
tic; simTag = 'FStarget'; fid = fopen('./elec72.loc'); 
C = textscan(fid,'%d %f %f %s'); fclose(fid); elecName = C{4}; 
for i=1:length(elecName), elecName{i} = strrep(elecName{i},'.',''); end
Nelec = length(elecName)-1; % Remove Reference (Iz)

e = Nelec-(ceil(idx/Nelec)*Nelec-idx);
eTag = elecName{e}; m = 1;
disp(['Working on ' subname ': ' eTag]);

if ~exist(fullfile(rootDir,subname,simTag),'dir'); mkdir(fullfile(rootDir,subname,simTag)); end
subj = fullfile(rootDir,subname,simTag,'T1.nii'); [dirname,baseFilename] = fileparts(subj);
copyfile(fullfile(fileparts(dirname),'ROAST_output','T1.nii'),fullfile(dirname,'T1.nii'))
copyfile(fullfile(fileparts(dirname),'ROAST_output','T1_T1orT2_masks.nii'),fullfile(dirname,'T1_T1orT2_masks.nii'))
copyfile(fullfile(fileparts(dirname),'ROAST_output','T1_T1orT2_seg8.mat'),fullfile(dirname,'T1_T1orT2_seg8.mat'))

while ~exist(fullfile(dirname,[baseFilename '_' eTag '_roastResult.mat']),'file')
    if m > lim; disp([subname ': ITERATION LIMIT REACHED FOR ' eTag]); break; end % QUIT AFTER 5 ITERATIONS
    try; hiperskipseg(idx,[],subj,'elec',elecName([e,end])'); end
    if ~exist(fullfile(dirname,[baseFilename '_' eTag '_roastResult.mat']),'file') % START FRESH
        disp(['RESTARTING ' subname ': ' eTag]); m = m + 1;
        delete(fullfile(dirname,[baseFilename '_' eTag '_*.*']));
    end
end
toc