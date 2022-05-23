%% Convert ROAST MESH to Voxels
% Created by Alejandro Albizu on 06/08/2020
% Last Updated: 08/31/2020 by AA
function exhaustiveSearch(sub,e)
if ~isa(e,'double'); e = str2double(e); end
if ~isa(sub,'double'); sub = str2double(sub); end

% Settings
%--------------------------
rootDir = '/blue/camctrp/working/aprinda/freesurfer_output';
% rootDir = 'L:\working\aprinda\freesurfer_output';
I_max = 4; % Max Intensity of Precision Dose (mA)
simTag = 'FStarget';
logname = ['exhaustiveSearch_sub-' num2str(sub) '_' datestr(now,30) '_log.txt'];
%--------------------------

subj = fullfile(rootDir,['FS6.0_sub-' num2str(sub) '_ses01'],simTag,'T1.nii');
[dirname,baseFilename] = fileparts(subj);

mA = 1:0.1:I_max; % HARD-CODED

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
if exist(fullfile(dirname,[baseFilename '_leadField.mat']),'file')
    load(fullfile(dirname,[baseFilename '_leadField.mat']),'A_all')
    A = A_all(:,[elecIdx(1,e) elecIdx(2,e)]); A(isnan(A)) = 0; % remove NaN
    clear A_all; % save RAM
    
    % disp('Inverse Normalization: Warping Weight map to Native Subject Space ...')
    % % Inverse Normalize: Weight map
%     matlabbatch{1}.spm.util.defs.comp{1}.inv.comp{1}.def = {fullfile(rootDir,['FS6.0_sub-' num2str(sub) '_ses01'],'y_T1.nii')};
%     matlabbatch{1}.spm.util.defs.comp{1}.inv.space = {fullfile(rootDir,['FS6.0_sub-' num2str(sub) '_ses01'],'T1.nii')};
%     matlabbatch{1}.spm.util.defs.out{1}.pull.fnames = {fullfile(rootDir,'DEFCON','pos_weights_MNI.nii')};
%     matlabbatch{1}.spm.util.defs.out{1}.pull.savedir.saveusr = {fullfile(rootDir,['FS6.0_sub-' num2str(sub) '_ses01'],'defcon')};
%     matlabbatch{1}.spm.util.defs.out{1}.pull.interp = 1;
%     matlabbatch{1}.spm.util.defs.out{1}.pull.mask = 1;
%     matlabbatch{1}.spm.util.defs.out{1}.pull.fwhm = [0 0 0];
%     matlabbatch{1}.spm.util.defs.out{1}.pull.prefix = 'iwuf';
    
    % % Inverse Normalize: Responder Average
    %    matlabbatch{2}.spm.util.defs.comp{1}.inv.comp{1}.def = {fullfile(rootDir,['FS6.0_sub-' num2str(sub) '_ses01'],'y_T1.nii')};
    %    matlabbatch{2}.spm.util.defs.comp{1}.inv.space = {fullfile(rootDir,['FS6.0_sub-' num2str(sub) '_ses01'],'T1.nii')};
    %    matlabbatch{2}.spm.util.defs.out{1}.pull.fnames = {fullfile(rootDir,'DEFCON','RESP_yEFavg.nii')};
    %    matlabbatch{2}.spm.util.defs.out{1}.pull.savedir.saveusr = {fullfile(rootDir,['FS6.0_sub-' num2str(sub) '_ses01'],'defcon')};
    %    matlabbatch{2}.spm.util.defs.out{1}.pull.interp = 1;
    %    matlabbatch{2}.spm.util.defs.out{1}.pull.mask = 1;
    %    matlabbatch{2}.spm.util.defs.out{1}.pull.fwhm = [0 0 0];
    %    matlabbatch{2}.spm.util.defs.out{1}.pull.prefix = 'iwuf';
    
    % % Inverse Normalize: Responder Standard Deviation
    % matlabbatch{3}.spm.util.defs.comp{1}.inv.comp{1}.def = {fullfile(rootDir,['FS6.0_sub-' num2str(sub) '_ses01'],'y_T1.nii')};
    % matlabbatch{3}.spm.util.defs.comp{1}.inv.space = {fullfile(rootDir,['FS6.0_sub-' num2str(sub) '_ses01'],'T1.nii')};
    % matlabbatch{3}.spm.util.defs.out{1}.pull.fnames = {fullfile(rootDir,'DEFCON','RESP_EFstd.nii')};
    % matlabbatch{3}.spm.util.defs.out{1}.pull.savedir.saveusr = {fullfile(rootDir,['FS6.0_sub-' num2str(sub) '_ses01'],'defcon')};
    % matlabbatch{3}.spm.util.defs.out{1}.pull.interp = 1;
    % matlabbatch{3}.spm.util.defs.out{1}.pull.mask = 1;
    % matlabbatch{3}.spm.util.defs.out{1}.pull.fwhm = [0 0 0];
    % matlabbatch{3}.spm.util.defs.out{1}.pull.prefix = 'iwuf';
    % if ~exist(fullfile(fileparts(dirname),'defcon','iwufpos_weights_MNI.nii'),'file') ...
    %         || ~exist(fullfile(fileparts(dirname),'defcon','RESP_EFavg.nii'),'file') ...
    %         || ~exist(fullfile(fileparts(dirname),'defcon','RESP_EFstd.nii'),'file')
    %     spm_jobman('run',matlabbatch)
    % end
    
    masks = load_nii(fullfile(dirname, [baseFilename '_T1orT2_masks.nii']));
    allMask = double(masks.img); clear masks;
    Nvox = numel(allMask);
    
    wnii = load_nii(fullfile(fileparts(dirname),'defcon','iwufpos_weights_MNI.nii'));
    native_w = repmat(wnii.img(:),3,1)'; native_w(isnan(native_w)) = 0;
    avg = load_nii(fullfile(fileparts(dirname),'defcon','iwufRESP_EFavg.nii'));
    gmean = avg.img(:); gmean(isnan(gmean)) = 0;
    gmean = gmean./sqrt(sum(avg.img(:).^2,2)); % unit vector
    sd = load_nii(fullfile(fileparts(dirname),'defcon','iwufRESP_EFstd.nii'));
    gstd = sd.img(:); clear wnii avg; % save RAM
    gstd = gstd./sqrt(sum(sd.img(:).^2,2)); % unit vector
    
    % Create Explicit Masks
    brain = repmat(allMask(:) == 1 | allMask(:) == 2,3,1)';
%     wmask = repmat(allMask(:) == 2,3,1)'; 
    wmask =  brain & native_w ~= 0;
    
    % Convert Voxel Coord to Mesh Nodes
    disp('Computing the optimized electric field (this may take a while) ...');
    for i = 1:length(mA)
        I_opt = [-mA(i);mA(i)];
        fJbrain = A * I_opt; fJbrain(brain == 0) = 0;
        fJbrain = fJbrain./repmat(sqrt(sum([fJbrain(1:Nvox) fJbrain(Nvox+1:Nvox*2) fJbrain((Nvox*2)+1:end)].^2,2)),3,1);% unit vector
        %         fJbrain = sqrt(sum([fJbrain(1:Nvox) fJbrain(Nvox+1:Nvox*2) fJbrain((Nvox*2)+1:end)].^2,2)); % ef_mag
        
        if ~exist(logname,'file')
            msg = 'Beginning Exhaustive Search ...';
            save(logname,'msg','-ascii')
        end
        LL = mean(exp(-sum((fJbrain(wmask) - gmean(wmask)).^2 ./ (gstd(wmask).^2 + 1),2)),'omitnan'); % log-likelihood
        %         dp = fJbrain(wmask)' * gmean(wmask); % dot-product
        %         udp = (fJbrain(wmask)' * gmean(wmask))/(norm(fJbrain(wmask),1)+norm(gmean(wmask),1)% unit dot-product
%         [~,~,mse] = lscov(fJbrain(brain),gmean(brain),(native_w(brain)/max(native_w(brain)))); % weighted least squares
        
        % Save Results
        fid = fopen(logname, 'a');
        fprintf(fid, '%s %smA %s %smA: %s\n',elecName{elecIdx(1,e)},num2str(-mA(i)),elecName{elecIdx(2,e)},num2str(mA(i)),num2str(LL));
        fclose(fid);
                
        % Swap polarity for inverse montage
        ifJbrain = -fJbrain;
        LL = mean(exp(-sum((ifJbrain(wmask) - gmean(wmask)).^2 ./ (gstd(wmask).^2 + 1),2)),'omitnan'); % log-likelihood
%         [~,~,mse] = lscov(ifJbrain(brain),gmean(brain),(native_w(brain)/max(native_w(brain))));
        
        % Save Results
        fid = fopen(logname, 'a');
        fprintf(fid, '%s %smA %s %smA: %s\n',elecName{elecIdx(1,e)},num2str(mA(i)),elecName{elecIdx(2,e)},num2str(-mA(i)),num2str(LL));
        fclose(fid);
        %         disp(['Electrode ' elecName{elecIdx(1,e)} ' ' num2str(mA(i)) 'mA ' elecName{elecIdx(2,e)} ' ' num2str(-mA(i)) 'mA computed: ' num2str(LL) ' log-likelihood'])
    end
    close all;
end
toc
