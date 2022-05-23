%% Convert ROAST MESH to Voxels
% Created by Alejandro Albizu on 06/08/2020
% Last Updated: 11/05/2020 by AA
function exhaustiveSearch_v4(sub,e)
if ~isa(e,'double'); e = str2double(e); end
if ~isa(sub,'double'); sub = str2double(sub); end

% Settings
%--------------------------
rootDir = '/blue/camctrp/working/aprinda/freesurfer_output';
psy = [32 26 19 20 20 17 15 2 4 6 10 16 18 7]'; % ACC + RT
dtype = 'DI';
I_max = 4; % Max Intensity of Precision Dose (mA)
simTag = 'FStarget';
dims = [256 256 256 3];
logname = ['exhaustiveSearch_sub-' num2str(sub) '_' dtype '_log.txt'];
%--------------------------

rng default

subj = fullfile(rootDir,['FS6.0_sub-' num2str(sub) '_ses01'],simTag,'T1.nii');
[dirname,baseFilename] = fileparts(subj);

I = 0.1:0.1:I_max; % HARD-CODED

% Define Electrodes
fid = fopen('./elec72.loc'); C = textscan(fid,'%d %f %f %s');
fclose(fid); elecName = C{4}; for y=1:length(elecName), elecName{y} = strrep(elecName{y},'.',''); end
Nelec = length(elecName)-1; % Remove Reference (Iz)
elecIdx = zeros(2,71*(70/2)); butt = 1;
elecPair = triu(repmat(1:Nelec,Nelec,1),1);
for e1 = 1:Nelec
    for e2 = 1:Nelec
        if elecPair(e1,e2) ~= 0
            elecIdx(1,butt) = e1;
            elecIdx(2,butt) = e2;
            butt = butt + 1;
        end
    end
end
clear elecPair; % save RAM

tic
% Load and register electrodes to UFAB
if ~exist(fullfile(dirname,[baseFilename '_' elecName{elecIdx(1,e)} '_roastResult_xform_1mm.nii']),'file') || ...
        ~exist(fullfile(dirname,[baseFilename '_' elecName{elecIdx(2,e)} '_roastResult_xform_1mm.nii']),'file')
    melec = find([~exist(fullfile(dirname,[baseFilename '_' elecName{elecIdx(1,e)} '_roastResult_xform_1mm.nii']),'file'),...
        ~exist(fullfile(dirname,[baseFilename '_' elecName{elecIdx(2,e)} '_roastResult_xform_1mm.nii']),'file')]);
    for i = 1:length(melec)
        if exist(fullfile(dirname,[baseFilename '_' elecName{elecIdx(melec(i),e)} '_roastResult.mat']),'file')
            disp('FLIRTing Electrodes (this may take a while) ...')
            load(fullfile(dirname,[baseFilename '_' elecName{elecIdx(melec(i),e)} '_roastResult.mat']),'ef_all'); a = ef_all(:);
            img = reshape(a,dims); am = load_nii(fullfile(dirname,'T1_T1orT2_masks.nii'));
            img(am.img ~= 1 & am.img ~= 2) = 0;
            fl = make_nii(img); fl.hdr.hist = am.hdr.hist; fl.hdr.dime.dim = [4 dims 1 1 1];
            save_nii(fl,fullfile(dirname,[baseFilename '_' elecName{elecIdx(melec(i),e)} '_roastResult.nii']));
            flirt = sprintf(['ml fsl/6.0.3\n', ...
                'applyxfm4D %s %s %s %s -singlematrix\n', ...
                'gunzip %s'], ...
                ...
                fullfile(dirname,[baseFilename '_' elecName{elecIdx(melec(i),e)} '_roastResult.nii']), ...                    % native space 4D image
                fullfile(rootDir,'maps_587.nii'), ...                                 % UFAB template
                fullfile(dirname,[baseFilename '_' elecName{elecIdx(melec(i),e)} '_roastResult_xform_1mm.nii']), ... % Inverse Affine Matrix
                fullfile(fileparts(dirname),['sub-' num2str(sub) '_T1_omat_1mm.mat']), ...              % Output Filename
                ...
                fullfile(dirname,[baseFilename '_' elecName{elecIdx(melec(i),e)} '_roastResult_xform_1mm.nii.gz']));           % Gunzip Output
            system(flirt);
            clear a fl am; % save RAM
        else
            error([' Could not locate LEAD FIELD for ' elecName{elecIdx(melec(i),e)} ' ...'])
        end
    end
end

% Define Group Labels as above/below Median
label = zeros(length(psy),1);
label(psy >= median(psy)) = 1; % Responders
label(psy < median(psy)) = -1; % Non Responders

load(fullfile(rootDir,'DEFCON',['a_' dtype '_UFAB_ideal.mat']),'activ')
% fprec = (pos_weights ~= 0)';
masks = load_nii(fullfile(fileparts(dirname), ['sub-' num2str(sub) '_allMask_xform_1mm.nii']));
allMask = double(masks.img); clear masks;
Nvox = numel(allMask); nd = length(a)/Nvox;
brain = repmat(allMask(:) == 1 | allMask(:) == 2,3,1);
weights = weights(brain); % shrink weights to size of group brain

% CREATE GMM WEIGHT VECTOR
%=====================================================================================
fprec = zeros(Nvox*nd,1);
fprec(brain(1:Nvox*nd)) = ...
    ( weights' .* (sum(weights ~= 0) / length(weights)) ) + ... % Scale SVM weights
    ( 1 - sum(weights' .* (sum(weights ~= 0) / length(weights))) ) ... % Distribute weight amongst remaining voxels
    / length(weights);
%=====================================================================================             

A = zeros(Nvox*3,2);
for i = 1:size(elecIdx,1)
    nii = load_nii(fullfile(dirname,[baseFilename '_' elecName{elecIdx(i,e)} '_roastResult_xform_1mm.nii']));
    A(:,i) = nii.img(:);
end
clear mag nii

% Fit Gaussian Model to model J-map of Responders (Gaussian mixture model (GMM) is possible with more subjects)
load(fullfile(rootDir,'DEFCON',['alldata_UFAB587_ideal_' dtype '.mat']),'alldata')
J_R = alldata(label == 1,:); g_mean = mean(J_R,1)'; g_std = std(J_R,[],1)'; 
% g_mean(aparc) = g_mean(aparc)+(abs(min(g_mean(aparc)))+1); g_std(aparc) = g_std(aparc)+(abs(min(g_std(aparc)))+1); 
% g_mean = g_mean .* fprec; g_std = g_std .* fprec;
% clear alldata J_R;

% Precision Current Prediction
disp('Computing the optimized electric field ...');
LL = zeros(length(I),2);
for i = 1:length(I)  
    J_sim = (A*[I(i);-I(i)]);
    J_sim(isnan(J_sim)) = 0; J_sim(~brain) = 0; % J_sim = J_sim .* fprec;
    if strcmp(dtype,'I'); J_sim = sqrt(J_sim(1:Nvox).^2 + J_sim(Nvox+1:Nvox*2).^2 + J_sim((Nvox*2)+1:end).^2);end
    J_sim(repmat(allMask(:) == 1,nd,1)) = J_sim(repmat(allMask(:) == 1,nd,1)) .* 0.126; % WM
    J_sim(repmat(allMask(:) == 2,nd,1)) = J_sim(repmat(allMask(:) == 2,nd,1)) .* 0.276; % GM
%     J_sim(brain) = J_sim(brain)+(abs(min(J_sim(brain)))+1); % Standardize

    % log-likelihood
    LL(i,:) = [exp(-sum((((J_sim' - g_mean')).^2 ./ (g_std'.^2 + 1)).*fprec',2)), ... % 
        exp(-sum((((-J_sim' - g_mean')).^2 ./ (g_std'.^2 + 1)).*fprec',2))]; 
end
[maxcol,~] = max(LL);
pol = find(maxcol-max(maxcol) == 0);
[maxLL,ind_LL] = max(LL(:,pol));
Iprec = I(ind_LL);

% Write out results
if ~exist(logname,'file')
    msg = 'Beginning Exhaustive Search ...';
    save(logname,'msg','-ascii')
end
fid = fopen(logname, 'a');
if pol == 1
    fprintf(fid, '%s,%smA,%s,%smA,%s\n',elecName{elecIdx(1,e)},num2str(Iprec),elecName{elecIdx(2,e)},num2str(-Iprec),num2str(maxLL));
else
    fprintf(fid, '%s,%smA,%s,%smA,%s\n',elecName{elecIdx(1,e)},num2str(-Iprec),elecName{elecIdx(2,e)},num2str(Iprec),num2str(maxLL));
end
fclose(fid);
toc
