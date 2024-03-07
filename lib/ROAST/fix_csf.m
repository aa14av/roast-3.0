%% CSF Correction Function
%=========================================================================================================
% Fix brain touching skull (this can create extreme FEM values).  
%
% OUTPUT:
% Saves current segmentation mask with _preCorr appended filename and saves
% corrected as segmentation mask that will be used by future steps of ROAST
%
% INPUT:
% T1 full filepath for determining current subject folder
% T2 full filepath for determining baseFilename
% cond 'conductivities' variable from ROAST containing conductivity, 
% tissue name, and tissue index information 
%
% Created by: Alejandro Albizu for the Center for Cognitive Aging and Memory
% University of Florida
% Email: aa14av@gmail.com
% Created: 02/02/2022
% Updated: 03/07/2024
% 20230324: Replaced `cond.index` with `find` to remove dependency
%=========================================================================================================
function fix_csf(T1,T2,cond)

% File Path
[dirname,baseFilenameRasRSPD] = fileparts(T1);
if isempty(dirname), dirname = pwd; end
if isempty(T2)
    baseFilenameRasRSPD = [baseFilenameRasRSPD '_T1orT2'];
else
    baseFilenameRasRSPD = [baseFilenameRasRSPD '_T1andT2'];
end

% Get Tissue Names
mnames = fieldnames(cond);
mnames = mnames(1:end-4); % Remove Non-Tissue

% Throw error if cannot locate gray matter label
if sum(cellfun(@(x) sum(regexpi(x,'gm|gray|grey')),mnames,'uni',1)) == 0
    error('Cannot locate gray matter conductivity in ''conductivities'' variable...')
end

% Throw error if cannot locate bone labels
if sum(cellfun(@(x) sum(regexpi(x,'bone|skull|cancellous|cortical')),mnames,'uni',1)) == 0
    error(sprintf(['Cannot locate bone conductivity in ''conductivities''',...
        ' variable... \nPlease relabel [options: ''bone'', ''skull'',',...
        ' ''cancellous'', ''cortical'']']))
end

% Throw error if cannot locate CSF labels
if sum(cellfun(@(x) sum(regexpi(x,'csf|cerebrospinal')),mnames,'uni',1)) == 0
    error(sprintf(['Cannot locate CSF conductivity in ''conductivities''',...
        ' variable... \nPlease relabel [options: ''csf'',',...
        ' ''cerebrospinal fluid'']']))
end

% Load Segmentations
am = load_untouch_nii(fullfile(dirname,[baseFilenameRasRSPD '_masks.nii']));

masks = am.img;
gm_mask = ismember(masks, find(cellfun(@(x) sum(regexpi(x,'gm|gray|grey')),mnames,'uni',1)~=0)); % Get Gray Matter
wm_mask = ismember(masks, find(cellfun(@(x) sum(regexpi(x,'wm|white')),mnames,'uni',1)~=0)); % Get White Matter
bone_mask = ismember(masks, find(cellfun(@(x) sum(regexpi(x,'bone|skull|cancellous|cortical')),mnames,'uni',1)~=0)); % Get Combined Bone

% Locate GM/bone Intersection
dil_bone = imdilate(bone_mask,ones(3,3,3)); % Dilate bone mask
dil_bone(bone_mask) = 0; % Dilated portion of bone mask ONLY

% Replace GM touching skull with csf
if length(find(cellfun(@(x) sum(regexpi(x,'csf|cerebrospinal')),mnames,'uni',1)~=0))>1
    error(['Loacted ', ...
        num2str(length(find(cellfun(@(x) sum(regexpi(x,'csf|cerebrospinal')),mnames,'uni',1)~=0))), ...
        ' csf conductivities in ''condutivities'' variable, expected 1...'])
end
masks(gm_mask&dil_bone) = find(cellfun(@(x) sum(regexpi(x,'csf|cerebrospinal')),mnames,'uni',1)~=0); % Replace GM touching Bone with CSF
masks(wm_mask&dil_bone) = find(cellfun(@(x) sum(regexpi(x,'csf|cerebrospinal')),mnames,'uni',1)~=0); % Replace WM touching Bone with CSF

nii = am; nii.img = masks;

% Save Corrected Segmentation
try
    save_untouch_nii(nii, fullfile(dirname,[baseFilenameRasRSPD '_masks.nii']));
catch ME
    warning(ME.message)
    save_nii(nii,fullfile(dirname,[baseFilenameRasRSPD '_masks.nii'])); 
end

% Save Old Segmentations
try 
    save_untouch_nii(am,fullfile(dirname,[erase(baseFilenameRasRSPD,'.nii') '_masks_preCorr.nii'])); 
catch ME
    warning(ME.message)
    save_nii(am,fullfile(dirname,[erase(baseFilenameRasRSPD,'.nii') '_masks_preCorr.nii']));
end