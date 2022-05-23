%% Convert ROAST MESH to Voxels
% Created by Alejandro Albizu on 06/08/2020
% Last Updated: 09/22/2020 by AA
function exhaustiveSearch_v2(sub,e)
if ~isa(e,'double'); e = str2double(e); end
if ~isa(sub,'double'); sub = str2double(sub); end

% Settings
%--------------------------
rootDir = '/blue/camctrp/working/aprinda/freesurfer_output';
psy = [32 26 19 20 20 17 15 2 4 6 10 16 18 7]'; % ACC + RT
dtype = 'DI';
I_max = 4; % Max Intensity of Precision Dose (mA)
% pslt = 0.5;
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

load(fullfile(rootDir,'DEFCON',['w_' dtype '_UFAB_ideal.mat']),'weights')
% fprec = (pos_weights ~= 0)';
masks = load_nii(fullfile(fileparts(dirname), ['sub-' num2str(sub) '_allMask_xform_1mm.nii']));
allMask = double(masks.img); clear masks;
Nvox = numel(allMask); nd = length(weights)/Nvox;
brain = repmat(allMask(:) == 1 | allMask(:) == 2,3,1);
weights = abs(weights(brain));
fprec = zeros(Nvox*nd,1);
% fprec = (weights'.*(sum(weights ~= 0)/length(weights)))+((1-sum((weights'.*(sum(weights ~= 0)/length(weights)))))/length(weights));
fprec(brain(1:Nvox*nd)) = (weights'.*(sum(weights ~= 0)/length(weights)))+((1-sum((weights'.*(sum(weights ~= 0)/length(weights)))))/length(weights));
% fprec = (weights'/100*pslt);
% fprec(weights ~= 0) = fprec(weights ~= 0)+((1-pslt)/sum(weights ~=0));

A = zeros(Nvox*3,2);
for i = 1:size(elecIdx,1)
    nii = load_nii(fullfile(dirname,[baseFilename '_' elecName{elecIdx(i,e)} '_roastResult_xform_1mm.nii']));
    A(:,i) = nii.img(:);
end
clear mag nii

% Fit Gaussian Model to model J-map of Responders (Gaussian mixture model (GMM) is possible with more subjects)
load(fullfile(rootDir,'DEFCON',['alldata_UFAB587_ideal_' dtype '.mat']),'alldata')
J_R = alldata(label == 1,:); g_mean = mean(J_R,1)'; g_std = std(J_R,[],1)';
% g_mean = g_mean .* fprec; g_std = g_std .* fprec;
% clear alldata J_R;

% Precision Current Prediction
disp('Computing the optimized electric field ...');
LL = zeros(length(I),2);
for i = 1:length(I)  
    J_sim = (A*[I(i);-I(i)]);
    J_sim(isnan(J_sim)) = 0; J_sim(~brain) = 0; % J_sim = J_sim .* fprec;
    if strcmp(dtype,'I'); J_sim = sqrt(J_sim(1:Nvox).^2 + J_sim(Nvox+1:Nvox*2).^2 + J_sim((Nvox*2)+1:end).^2);end
    J_sim(repmat(allMask(:) == 1,nd,1)) = J_sim(repmat(allMask(:) == 1,nd,1)) .* 0.126;
    J_sim(repmat(allMask(:) == 2,nd,1)) = J_sim(repmat(allMask(:) == 2,nd,1)) .* 0.276;
    
    % log-likelihood
    LL(i,:) = [exp(-sum((((J_sim' - g_mean')).^2 ./ (g_std'.^2 + 1)).*fprec',2)), ... % 
        exp(-sum((((-J_sim' - g_mean')).^2 ./ (g_std'.^2 + 1)).*fprec',2))]; 
end
[maxcol,~] = max(LL);
pol = find(maxcol-max(maxcol) == 0);
[maxLL,ind_LL] = max(LL(:,pol));
Iprec = I(ind_LL);

% PLOTTING (COMMENT BEFORE COMPILING) !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
%==========================================================================
% Standardize Vars
% mu = 100; sig = 15;
% J_prec = (A*[I(ind_LL);-I(ind_LL)]);
% J_prec(repmat(allMask(:) == 1,nd,1)) = J_prec(repmat(allMask(:) == 1,nd,1)) .* 0.126;
% J_prec(repmat(allMask(:) == 2,nd,1)) = J_prec(repmat(allMask(:) == 2,nd,1)) .* 0.276;
% prec = NaN([1 Nvox*nd]);
% prec(brain) = zscore(J_prec(brain))*sig+mu; 
% rmu = NaN([1 Nvox*nd]);
% rmu(brain) = zscore(g_mean(brain))*sig+mu;
% 
% % Compute
% bidx = find(brain);
% sim_norm = NaN(Nvox,1);
% sim_norm(brain(1:Nvox)) = zscore(sqrt(sum(reshape(J_prec(brain).^2,[sum(brain)/3 3]),2)))*sig+mu;
% rsp = reshape(sim_norm,[121 145 118]);
% % nonrJ = alldata(label == -1,:);
% % rsi = reshape(nonrJ(2,:),[121 145 118 nd]);
% r_norm = NaN(Nvox,1);
% r_norm(brain(1:Nvox)) = zscore(sqrt(sum(reshape(g_mean(brain).^2,[sum(brain)/3 3]),2)))*sig+mu;
% rsg = reshape(r_norm, [121 145 118]);
% diff = rsp-rsg;
% range = [mu-sig*3 mu+sig*3];
% % subplot(1,2,1); imagesc(rsp(:,:,70,di),range); axis off; axis tight; title('Simulated'); subplot(1,2,2); imagesc(rsi(:,:,70,di),range); axis off; axis tight; title('ROAST');
% 
% % Plot
% figure; cm = colormap('jet'); % cm(1,:) = [0 0 0];
% s1=subplot(2,2,1); imagesc(rsp(:,:,70),range); axis off; axis tight; title(['x: ' elecName{elecIdx(1,e)} '+' elecName{elecIdx(2,e)} ' @ ' num2str(I(ind_LL)) ' mA']); colormap(s1,cm); colorbar
% s2=subplot(2,2,2); imagesc(rsg(:,:,70),range); axis off; axis tight; title('Responder Mean (\mu)'); colormap(s2,cm); colorbar 
% if pol ~= 1; prec = -prec; end
% subplot(2,2,3); scatter(zscore(prec(bidx(weights ~= 0)))*sig+mu,zscore(rmu(bidx(weights ~= 0)))*sig+mu,'filled'); lsline; xlim(range); ylim(range); xlabel('J_{x} of weighted voxels'); ylabel('J_{\mu} of weighted voxels');
% title(['r = ' num2str(corr(zscore(prec(bidx(weights ~= 0))')+sig+mu,zscore(rmu(bidx(weights ~= 0))')*sig+mu)) ' + ' num2str(round(maxLL,3)*100) '% likelihood']);
% s4=subplot(2,2,4); imagesc(diff(:,:,70),[-sig*2 sig*2]); axis off; axis tight; title(['x - ' char(181)]); colormap(s4,'gray'); colorbar
%==========================================================================

% Write out results
if ~exist(logname,'file')
    msg = 'Beginning Exhaustive Search ...';
    save(logname,'msg','-ascii')
end
fid = fopen(logname, 'a');
if pol == 1
    fprintf(fid, '%s,%smA,%s,%smA,%s\n',elecName{elecIdx(1,e)},num2str(-Iprec),elecName{elecIdx(2,e)},num2str(Iprec),num2str(maxLL));
else
    fprintf(fid, '%s,%smA,%s,%smA,%s\n',elecName{elecIdx(1,e)},num2str(Iprec),elecName{elecIdx(2,e)},num2str(-Iprec),num2str(maxLL));
end
fclose(fid);
toc
