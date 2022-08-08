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
% Updated: 08/06/2022
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

mnames = cond(:,3);

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
am = load_nii(fullfile(dirname,baseFilenameRasRSPD));
save_nii(am,fullfile(dirname,[erase(baseFilenameRasRSPD,'.nii') '_preCorr.nii'])); % Save Old Segmentations
masks = am.img;
gm_mask = ismember(masks, [cond{cellfun(@(x) sum(regexpi(x,'gm|gray|grey')),mnames,'uni',1)~=0,1}]); % Get Gray Matter
bone_mask = ismember(masks, [cond{cellfun(@(x) sum(regexpi(x,'bone|skull|cancellous|cortical')),mnames,'uni',1)~=0,1}]); % Get Combined Bone

% Locate GM/bone Intersection
dil_bone = imdilate(bone_mask,ones(3,3,3)); % Dilate bone mask
dil_bone(bone_mask) = 0; % Dilated portion of bone mask ONLY

% Replace GM touching skull with csf
if length([cond{cellfun(@(x) sum(regexpi(x,'csf|cerebrospinal')),mnames,'uni',1)~=0,1}])>1
    error(['Loacted ', ...
        num2str(length([cond{cellfun(@(x) sum(regexpi(x,'csf|cerebrospinal')),mnames,'uni',1)~=0,1}])), ...
        ' csf conductivities in ''condutivities'' variable, expected 1...'])
end
masks(gm_mask&dil_bone) = cond{cellfun(@(x) sum(regexpi(x,'csf|cerebrospinal')),mnames,'uni',1)~=0,1}; % Replace GM touching Bone with CSF

nii = make_nii(masks); 
nii.hdr = am.hdr;
filename = fullfile(dirname,baseFilenameRasRSPD);
save_nii(nii,filename); % Save Corrected Segmentation
