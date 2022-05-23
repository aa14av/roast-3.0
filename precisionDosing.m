%% Convert ROAST MESH to Voxels
% Created by Alejandro Albizu on 06/08/2020
% Last Updated: 08/04/2020 by AA
cd /blue/camctrp/working/Alejandro/toolboxes/roast-3.0
addpath(genpath(pwd))
addpath ../../Scripts/
addpath ../spm12/

spm fmri
close all

subs = [9009 9015 9022 9040 9044 9045 9051 9021 9023 9031 9032 9047 9048 9054];
psy = [32 26 19 20 20 17 15 2 4 6 10 16 18 7]'; % ACC + RT

label = zeros(length(psy),1);
label(psy >= median(psy)) = 1; % Responders
label(psy < median(psy)) = -1; % Non Responders
%%
clear
% if ~isa(e,'double'); e = str2double(e); end
% if ~isa(i,'double');i = str2double(i); end
% if ~isa(sub,'double'); sub = str2double(sub); end

% Settings
%--------------------------
rootDir = '/blue/camctrp/working/aprinda/freesurfer_output';
% SVMDir = fullfile(rootDir,'DEFCON');
I_max = 4; % Max Intensity of Precision Dose (mA)
minElecs = 2; % Min # of Electrodes in Precision Dose
maxElecs = 2; % Max # of Electrodes in Precision Dose
% norm_dims = [182 218 182];
sub = 9023;
simTag = 'FStarget';
% elecArray = {{'Fp'},{'AF'},{'F'},{'FC','FT'},{'C','T'},{'CP','TP'},{'P'},{'PO'},{'O'}};
%--------------------------

% ADD hiperROAST_target_gui function %

% load(fullfile(SVMDir,'alldata_I_MNI.mat'),'alldata');
% load(fullfile(SVMDir,'w_I.mat'),'pos_weights')
% respJ = zeros(sum(label == 1), prod(norm_dims));
% g_mean = zeros(1,prod(norm_dims)); g_std = zeros(1,prod(norm_dims));
% respJ(:,pos_weights ~= 0) = alldata(label == 1,pos_weights ~= 0);
% g_mean(:,pos_weights ~= 0) = mean(respJ(:,pos_weights ~= 0),1); % mean of Gaussian model of key features
% g_std(:,pos_weights ~= 0) = std(respJ(:,pos_weights ~= 0),1); % standard deviation of key features
% LL_R = exp(-sum((respJ - g_mean).^2 ./ (g_std.^2 + 1),2)); % likelihood
subj = fullfile(rootDir,['FS6.0_sub-' num2str(sub) '_ses01'],simTag,'T1.nii');
[dirname,baseFilename] = fileparts(subj);

% mA = 1:0.1:max_I; % HARD-CODED
%
% % Define Electrodes
fid = fopen('./elec72.loc'); C = textscan(fid,'%d %f %f %s');
fclose(fid); elecName = C{4}; for y=1:length(elecName), elecName{y} = strrep(elecName{y},'.',''); end
% Nelec = length(elecName)-1; % Remove Reference (Iz)
% elecIdx = zeros(2,71*(70/2)); butt = 1;
% elecPair = triu(repmat(1:Nelec,Nelec,1),1);
% for e1 = 1:Nelec
%     for e2 = 1:Nelec
%         if elecPair(e1,e2) ~= 0
%             elecIdx(1,butt) = e1;
%             elecIdx(2,butt) = e2;
%             butt = butt + 1;
%         end
%     end
% end
tic
if exist(fullfile(dirname,[baseFilename '_' simTag '_roastResult.mat']),'file')
%     load(fullfile(dirname,[baseFilename '_' simTag '_roastResult.mat']),'A_all');
    load(fullfile(dirname,[baseFilename '_leadField.mat']),'A_all')
%     load(fullfile(dirname,[baseFilename '_' simTag '.mat']),'node','elem');
%     load(fullfile(dirname,[baseFilename '_header.mat']),'hdrInfo');
    
    disp('Inverse Normalization: Warping Weight map to Native Subject Space ...')
    % Inverse Normalize: Weight map
    matlabbatch{1}.spm.util.defs.comp{1}.inv.comp{1}.def = {['/blue/camctrp/working/aprinda/freesurfer_output/FS6.0_sub-' num2str(sub) '_ses01/y_T1.nii']};
    matlabbatch{1}.spm.util.defs.comp{1}.inv.space = {['/blue/camctrp/working/aprinda/freesurfer_output/FS6.0_sub-' num2str(sub) '_ses01/T1.nii']};
    matlabbatch{1}.spm.util.defs.out{1}.pull.fnames = {'/blue/camctrp/working/aprinda/freesurfer_output/DEFCON/pos_weights_MNI.nii'};
    matlabbatch{1}.spm.util.defs.out{1}.pull.savedir.saveusr = {['/blue/camctrp/working/aprinda/freesurfer_output/FS6.0_sub-' num2str(sub) '_ses01/defcon']};
    matlabbatch{1}.spm.util.defs.out{1}.pull.interp = 1;
    matlabbatch{1}.spm.util.defs.out{1}.pull.mask = 1;
    matlabbatch{1}.spm.util.defs.out{1}.pull.fwhm = [0 0 0];
    matlabbatch{1}.spm.util.defs.out{1}.pull.prefix = 'iwuf';
    
    % Inverse Normalize: Responder Average
    matlabbatch{2}.spm.util.defs.comp{1}.inv.comp{1}.def = {['/blue/camctrp/working/aprinda/freesurfer_output/FS6.0_sub-' num2str(sub) '_ses01/y_T1.nii']};
    matlabbatch{2}.spm.util.defs.comp{1}.inv.space = {['/blue/camctrp/working/aprinda/freesurfer_output/FS6.0_sub-' num2str(sub) '_ses01/T1.nii']};
    matlabbatch{2}.spm.util.defs.out{1}.pull.fnames = {'/blue/camctrp/working/aprinda/freesurfer_output/DEFCON/RESP_EFavg.nii'};
    matlabbatch{2}.spm.util.defs.out{1}.pull.savedir.saveusr = {['/blue/camctrp/working/aprinda/freesurfer_output/FS6.0_sub-' num2str(sub) '_ses01/defcon']};
    matlabbatch{2}.spm.util.defs.out{1}.pull.interp = 1;
    matlabbatch{2}.spm.util.defs.out{1}.pull.mask = 1;
    matlabbatch{2}.spm.util.defs.out{1}.pull.fwhm = [0 0 0];
    matlabbatch{2}.spm.util.defs.out{1}.pull.prefix = 'iwuf';
    
    % Inverse Normalize: Responder Standard Deviation
    matlabbatch{3}.spm.util.defs.comp{1}.inv.comp{1}.def = {['/blue/camctrp/working/aprinda/freesurfer_output/FS6.0_sub-' num2str(sub) '_ses01/y_T1.nii']};
    matlabbatch{3}.spm.util.defs.comp{1}.inv.space = {['/blue/camctrp/working/aprinda/freesurfer_output/FS6.0_sub-' num2str(sub) '_ses01/T1.nii']};
    matlabbatch{3}.spm.util.defs.out{1}.pull.fnames = {'/blue/camctrp/working/aprinda/freesurfer_output/DEFCON/RESP_EFstd.nii'};
    matlabbatch{3}.spm.util.defs.out{1}.pull.savedir.saveusr = {['/blue/camctrp/working/aprinda/freesurfer_output/FS6.0_sub-' num2str(sub) '_ses01/defcon']};
    matlabbatch{3}.spm.util.defs.out{1}.pull.interp = 1;
    matlabbatch{3}.spm.util.defs.out{1}.pull.mask = 1;
    matlabbatch{3}.spm.util.defs.out{1}.pull.fwhm = [0 0 0];
    matlabbatch{3}.spm.util.defs.out{1}.pull.prefix = 'iwuf';
    spm_jobman('run',matlabbatch)
    
    masks = load_nii(fullfile(dirname, [baseFilename '_T1orT2_masks.nii']));
    allMask = double(masks.img);
    
    wnii = load_nii(fullfile(fileparts(dirname),'defcon','iwufpos_weights_MNI.nii'));
    native_w = wnii.img(:);
    avg = load_nii(fullfile(fileparts(dirname),'defcon','iwufRESP_EFavg.nii'));
    gmean = avg.img(:);
    sd = load_nii(fullfile(fileparts(dirname),'defcon','iwufRESP_EFstd.nii'));
    gstd = sd.img(:);
    
    [ndx,ndy,ndz] = ndgrid(1:hdrInfo.dim(1),1:hdrInfo.dim(2),1:hdrInfo.dim(3));
    coord(:,1) = ndx(:); coord(:,2) = ndy(:); coord(:,3) = ndz(:);
    wmask = native_w > median(native_w((allMask(:) == 1 | allMask(:) == 2) & native_w > 0)); % top 50%
    targetCoord = coord(wmask,:); % Targets
    numTargs = length(targetCoord);
    
    indBrain = elem((elem(:,5)==1 | elem(:,5)==2),1:4); indBrain = unique(indBrain(:));
    A = A_all(indBrain,:,:);
    
    % convert pseudo-world coordinates back to voxel coordinates for targeting,
    % as targeting code works in the voxel space
    nodeV = zeros(size(node,1),3);
    for i=1:3, nodeV(:,i) = node(:,i)/hdrInfo.pixdim(i); end
    locs = nodeV(indBrain,1:3);
    isNaNinA = isnan(sum(sum(A,3),2)); % make sure no NaN is in matrix A or in locs
    if any(isNaNinA), error('Recompute Lead Fields: leadFields contain NaN ...'); end
    Nlocs = size(locs,1);
    Nelec = size(A,3);
    A = reshape(A,Nlocs*3,Nelec);
    %     for i=1:3; locs(:,i) = rescale(locs(:,i),0,norm_dims(i)); end
    
    % Convert Voxel Coord to Mesh Nodes
    disp('Computing the optimized electric field (this may take a while) ...');
%     dist = pdist2(locs,targetCoord);
    Cf = zeros(Nelec,numTargs); % targetRadius = max(min(dist));
    parfor n = 1:numTargs
        dist = pdist2(locs,targetCoord(n,:));
        mindist(n) = min(dist);
    end
    targetRadius = max(mindist);
    disp(['Target Radius: ' num2str(targetRadius) 'mm'])
    
    nodeSearch = cell(numTargs,1);
    w = ones(3*Nlocs,1); d = ones(3*Nlocs,1);
    targW = native_w(wmask); targD = gmean(wmask); 
    for n = 1:numTargs
        dist = pdist2(locs,targetCoord(n,:));
        nodeSearch(n) = {find(dist <= targetRadius)};
%         w(nodeSearch{n}) = targW(n)+1;
%         d(nodeSearch{n}) = targD(n);
%         w(nodeSearch{n}+Nlocs) = targW(n)+1;
%         d(nodeSearch{n}+Nlocs) = targD(n);
%         w(nodeSearch{n}+(2*Nlocs)) = targW(n)+1;
%         d(nodeSearch{n}+(2*Nlocs)) = targD(n);
        %         nodeSearch(n) = rangesearch(locs,targetCoord(n,:),targetRadius);
        %         disp(['Target Coordinates: ' num2str(targetRadius) 'mm Sphere generated at voxel [' num2str(targetCoord(n,1)) ',' num2str(targetCoord(n,2)) ',' num2str(targetCoord(n,3)) ']'])
        %         disp(['Target Coordinates: Sphere contains ' num2str(length(nodeSearch{n})) ' nodes'])
        %     target_nodes{n} = find(distances_to_target(:,n)<targetRadius);
        if isempty(nodeSearch{n})
            error('No nodes found near target. Please increase the value of ''targetRadius''.');
        end
        Cx = mean(A(nodeSearch{n},:),1);
        Cy = mean(A(nodeSearch{n}+Nlocs,:),1);
        Cz = mean(A(nodeSearch{n}+2*Nlocs,:),1);
        C = [Cx;Cy;Cz];
        f = [1;1;1];
        Cf(:,n) = C'*f;
    end
%     padCf = Disk2Pad(Cf,elecName,elecArray,Nelec,numTargs)
    
%     % GRADIENT DESCENT WITH MOMENTUM
%     for i =1 : maxiter
%         % Derivative(Finding the residual)
%         r = (1/n)*(padCf*(padCf' * theta' - gmean(wmask)));
%         % Updating Theta
%         vt = gamma * vt + alpha * r';
%         theta = theta - vt;
%         % Mean-squared Error
%         mse = (1/(2*n)) * sum((padCf' * theta'-gmean(wmask)).^2);
%         % Relative change in Theta using squared norm
%         err_theta = norm(theta - t)/norm(theta);
%         % Display iteration info
%         disp(['Iteration ' num2str(i) ' : Mean-squared error = ' num2str(mse)]);
%         % Check for convergence
%         if err_theta <= tol
%             fprintf('\nChange in theta less tha specified tolerance\n')
%             break
%         end
%         t = theta;
%     end
    [U,S,V] = svd(repmat(sqrt(w),1,size(A,2)).*A, 0);
    [I_hat,~] = wls_dp(S*V',U'*(sqrt(w).*d),I_max,1);
    [I_hat,~] = maxLL_cvx(Cf',I_max,gmean(wmask)',gstd(wmask)',true);
%     mA = 1:0.1:4; elecIdx = find(I_hat);
%     allef = repmat(Cf(I_hat > 0,:),31,1) .* mA' + repmat(Cf(I_hat < 0,:),31,1) .* -mA';
%     for i = 1:length(mA); alldp(i) = allef(i,:) * gmean(wmask); end
%     [~,bestIdx] = max(alldp);
    fJbrain = (Cf * I_hat)';   
%     I_hat = I_hat.*median(1+((gmean(wmask)'-fJbrain)./fJbrain));
    disp(['Precision Dose: ' elecName{elecIdx(1)} ' ' num2str(I_hat(elecIdx(1))) 'mA + ' elecName{elecIdx(2)} ' ' num2str(I_hat(elecIdx(2))) 'mA'])
    udp = (fJbrain./norm(fJbrain,2))*(gmean(wmask)./norm(gmean(wmask),2)); dp = (fJbrain)*(gmean(wmask));
    disp(['Precision Dose: ' num2str(udp) ' (' num2str(dp) ') unit dot product for ' elecName{elecIdx(1)} ' + ' elecName{elecIdx(2)} ])
        
    % Save Precision Dose as NIFTI
    Jbrain = convertMesh2Vox(dirname,elecName,hdrInfo,allMask,nodeV,A_all,I_hat);
    %     percLL = exp(-sum((Jbrain(wmask) - gmean(wmask)).^2 ./ (gstd(wmask).^2 + 1),1)); % likelihood

    jn = make_nii(Jbrain); jn.hdr.hist = masks.hdr.hist;
    save_nii(jn,fullfile(dirname,['PRECISE_' elecName{elecIdx(1)} '_' elecName{elecIdx(2)} '_' num2str(norm(I_hat,1)/2) 'mA_EFbrain.nii']))
end
toc

%     masks = load_nii(fullfile(dirname, [baseFilename '_T1orT2_masks.nii']));
%     allMask = double(masks.img); % LL = zeros(Nelec,Nelec,length(mA)); %LL_C = zeros(Nelec,Nelec,length(mA));
%     jb = load_nii(fullfile(fileparts(dirname),'defcon','T1_Jbrain.nii'));
%     %         ideal = load_nii(fullfile(fileparts(dirname),'defcon','wufT1_Jbrain.nii'));
%     %         fJideal = double(ideal.img(:))'; % nanIdx = isnan(fJideal);
%
%     if e == 1 && i == 1; bestIdx = [0 0 0]; if exist(fullfile(dirname,'PRECISE.mat'),'file'); delete(fullfile(dirname,'PRECISE.mat')); end; end
%
%     % Reverse Calculate Jmap from leadField
%     I_opt = zeros(Nelec,1); I_opt([elecIdx(1,e) elecIdx(2,e)]) = [mA(i) -mA(i)];
%     Jbrain = convertMesh2Vox(dirname,elecName,hdrInfo,allMask,nodeV,A_all,I_opt);
%     jn = make_nii(Jbrain); jn.hdr = jb.hdr;
%     save_nii(jn,fullfile(dirname,['T1_' elecName{elecIdx(1,e)} '_' elecName{elecIdx(2,e)} '_' num2str(mA(i)) 'mA_Jbrain.nii']))
%
%     % SPM normalise to MNI with UFAB-587 TPM
%     matlabbatch = makeMLB(dirname,['T1_' elecName{elecIdx(1,e)} '_' elecName{elecIdx(2,e)} '_' num2str(mA(i)) 'mA_Jbrain']);
%     save(fullfile(dirname,['matlabbatch_' elecName{elecIdx(1,e)} '_' elecName{elecIdx(2,e)} '_' num2str(mA(i)) 'mA.mat']),'matlabbatch')
%     [~,res]=system(['/blue/camctrp/working/Alejandro/toolboxes/standalone/run_spm12.sh /apps/matlab/mcr/2017a/v92 batch ' fullfile(dirname,['matlabbatch_' elecName{elecIdx(1,e)} '_' elecName{elecIdx(2,e)} '_' num2str(mA(i)) 'mA.mat'])]);
%     disp(res)
%     % FLIRT new Jmap to UFAB-587
%     % affdr = dir(fullfile(fileparts(dirname),'*T1_omat*.mat'));
%     % if exist(fullfile(dirname,['T1_' elecName{elecIdx(1,e)} '_' elecName{elecIdx(2,e)} '_' num2str(mA) 'mA_Jbrain_xform_1mm.nii']),'file')
%     %     delete(fullfile(dirname,['T1_' elecName{elecIdx(1,e)} '_' elecName{elecIdx(2,e)} '_' num2str(mA) 'mA_Jbrain_xform_1mm.nii']))
%     % end
%     % system(sprintf('ml fsl/6.0.1\nflirt -interp nearestneighbour -in %s -ref /blue/camctrp/working/aprinda/freesurfer_output/maps_587.nii -applyxfm -init %s -out %s\ngunzip %s',fullfile(dirname,['T1_' elecName{elecIdx(1,e)} '_' elecName{elecIdx(2,e)} '_' num2str(mA) 'mA_Jbrain.nii']),fullfile(affdr.folder,affdr.name),fullfile(dirname,['T1_' elecName{elecIdx(1,e)} '_' elecName{elecIdx(2,e)} '_' num2str(mA) 'mA_Jbrain_xform_1mm.nii']),fullfile(dirname,['T1_' elecName{elecIdx(1,e)} '_' elecName{elecIdx(2,e)} '_' num2str(mA) 'mA_Jbrain_xform_1mm.nii.gz'])));
%
%     % Load Warped Jmap
%     jw = load_nii(fullfile(dirname,['wufT1_' elecName{elecIdx(1,e)} '_' elecName{elecIdx(2,e)} '_' num2str(mA(i)) 'mA_Jbrain.nii']));
%     Jwarp = double(jw.img); fJbrain = Jwarp(:)'; fJbrain(g_mean == 0) = 0;
%
%     % Compute Guassian Model
%     if exist(fullfile(dirname,'PRECISE.mat'),'file'); load(fullfile(dirname,'PRECISE.mat'),'bestIdx','bestLL'); end
%
%     if exp(-sum((fJbrain - g_mean).^2 ./ (g_std.^2 + 1),2)) > bestLL
%         bestIdx = [elecIdx(1,e) elecIdx(2,e) i];
%         bestLL = exp(-sum((fJbrain - g_mean).^2 ./ (g_std.^2 + 1),2)); % likelihood                                disp(['Computing: Electrode (' num2str(elecIdx(1,e)) ') ' elecName{elecIdx(1,e)} ' with Electrode (' num2str(elecIdx(2,e)) ') ' elecName{elecIdx(2,e)} ' @ ' num2str(mA(i)) 'mA'])
%         disp(['NEW BEST: Electrode (' num2str(elecIdx(1,e)) ') ' elecName{elecIdx(1,e)} ' + (' num2str(elecIdx(2,e)) ') ' elecName{elecIdx(2,e)} ' @ ' num2str(mA(i)) 'mA'])
%
%         % Save Results
%         save(fullfile(dirname,'PRECISE_.mat'),'bestIdx','bestLL')
%     end
% end
%
% if e == length(elecIdx) && i == length(mA); disp(['PRECISION DOSE: Electrode (' num2str(bestIdx(1)) ') ' elecName{bestIdx(1)} ' + (' num2str(bestIdx(2)) ') ' elecName{bestIdx(2)} ' @ ' num2str(mA(bestIdx(3))) 'mA']); end
