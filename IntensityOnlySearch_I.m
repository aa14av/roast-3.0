%%
cd /blue/camctrp/working/Alejandro/toolboxes/roast-3.0
addpath ../NIFTI_20130306/
addpath ../plotSpread/; addpath ../plot_topography/plot_topography/
addpath(genpath('/blue/camctrp/working/Alejandro/toolboxes/PRoNTo_dev-2.1.3/machines/libsvm-3.20/')); % LibSVM
%%
clear

rootDir = '/blue/camctrp/working/Alejandro/StimBrain';
outDir = fullfile(rootDir,'DEFCON_II');
% rootDir = 'L:\working\aprinda\freesurfer_output';
subs = [9051 9021 9023 9031 9032 9047 9054]; %1892 101549 300100 300700 301051]; % 9054 300700
psy = [32 26 19 20 20 17 15 2 4 6 10 16 18 7]'; % ACC + RT
I_max = 4; % Max Intensity of Precision Dose (mA)
dtype = 'DI';
p = 0.6666; % Percent contribution of SVM weights
uthr = 0.15; % Voxel Overlap Outlier Threshold
recipe = {'F3',-2,'F4',2}; % Original Montage
slice = 120; dims = [182 218 182];
simTag = 'FStarget';
%--------------------------

tic

% Define Group Labels as above/below Median
label = zeros(length(psy),1);
label(psy >= median(psy)) = 1; % Responders
label(psy < median(psy)) = -1; % Non Responders

I = 0.1:0.1:I_max; % HARD-CODED

% Define Electrodes
fid = fopen('./elec72.loc'); C = textscan(fid,'%d %f %f %s');
fclose(fid); elecName = C{4}; for y=1:length(elecName), elecName{y} = strrep(elecName{y},'.',''); end
Nelec = length(elecName)-1; % Remove Reference (Iz)

% Number of Dimensions
switch dtype; case 'I'; nd = 1; case 'DI'; nd = 3; end 

% Fit Gaussian Model to model J-map of Responders (Gaussian mixture model (GMM) is possible with more subjects)
load(fullfile(outDir,'interpdata_MNI_DI_VIII.mat'),'alldata')
J_RinterpY = alldata(label == 1,:); J_RinterpY(J_RinterpY == 0) = NaN; clear alldata;
r_mean = mean(J_RinterpY,1,'omitnan')'; r_std = std(J_RinterpY,[],1,'omitnan')';

% Get Non-Interp Data
load(fullfile(outDir,'alldata_MNI_DI_II.mat'),'alldata')
J_RinterpN = alldata(label == 1,:); J_RinterpN(J_RinterpN == 0) = NaN; clear alldata;

% Load SVM Weights
load(fullfile(outDir,'a_I_MNI.mat'),'pos_act','Mdl'); % Changed from v1 20211202
POSJ = mean(J_RinterpN(:,pos_act ~= 0),2,'omitnan');

% Union Mask
load(fullfile(outDir,'allAM_MNI_III.mat'),'allAM'); % Tissue Segmentations
umask = repmat(sum(allAM == 1 | allAM == 2)' > size(allAM,1)*uthr,nd,1);
clear allAM; % Save RAM

prs = zeros(length(subs),1); frs = zeros(length(subs),1); I_opt = zeros(length(subs),2);
headvol = zeros(length(subs),1); bratio = zeros(length(subs),1);
fLL = zeros(length(subs),1); LL = zeros(length(subs),1); acdist = zeros(length(subs),1);
for s = 1:length(subs)
    sub = subs(s); logname = ['exhaustiveSearch_sub-' num2str(sub) '_' dtype '_log.txt'];
    subj = fullfile(rootDir,['FS6.0_sub-' num2str(sub) '_ses01'],simTag,'T1.nii');
    [dirname,baseFilename] = fileparts(subj);
    
    % Load Brain Mask
    masks = niftiread(fullfile(fileparts(dirname),'ROAST_output','wufT1_T1orT2_masks.nii'));
    allMask = zeros(size(masks));
    bin = imfill(imerode(masks ~= 0 & ~isnan(masks),strel('cube',4)),'holes');
    allMask(bin) = masks(bin); clear masks; Nvox = numel(allMask(:)); 
    brain = repmat(allMask(:) == 1 | allMask(:) == 2,nd,1); % activ = abs(activ);
    
    % Create Weight Vector
    w = repmat(pos_act,1,nd)./3 * (p/100); w=w'; % Repeat along # dims
    w(brain & w == 0) = (1/sum(brain & w == 0))*(1-p); % Distribute remaining 50% to whole brain
    
    e1 = load_untouch_nii(fullfile(dirname,['wuf' baseFilename '_' recipe{1} '_roastResult.nii'])); % HARDCODED
    B(:,1) = e1.img(:);
    e2 = load_untouch_nii(fullfile(dirname,['wuf' baseFilename '_' recipe{3} '_roastResult.nii'])); % HARDCODED
    B(:,2) = e2.img(:); B = B*[recipe{2};recipe{4}]; B(~brain) = NaN;
    B(isnan(B)) = NaN; B(repmat(allMask(:) == 1,nd,1)) = B(repmat(allMask(:) == 1,nd,1)) .* 0.126;
    B(repmat(allMask(:) == 2,nd,1)) = B(repmat(allMask(:) == 2,nd,1)) .* 0.276; fJbrain = B(:);
    
    % Compute Intensity from LL
    %=================================================
%     N = zeros(length(I),1); %toc
%     for ii = 1:length(I)
%         x = gpuArray(fJbrain.*I(ii));
%         N(ii) = sum(normpdf(x,r_mean,r_std).*w,'omitnan');
% %         N(ii) = exp(-sum( ( (x - r_mean).^2  ./ (1 + r_std.^2) ) .* w,'omitnan'));
%         clear x;
%     end; [LL(s), indLL] = max(N); I_opt(s,:) = [I(indLL),I(indLL)] .* I_d';
    %=================================================

    % Compute Intensity from paper
    %=================================================
    I_opt(s,:) = round([2 -2] .* (1 + ( median(POSJ) - median(fJbrain(pos_act ~= 0),'omitnan') ) / median(fJbrain(pos_act ~= 0),'omitnan')),1); % IPrec from equation
    %=================================================
end; toc