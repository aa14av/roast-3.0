%% Convert ROAST MESH to Voxels
% Created by Alejandro Albizu on 06/08/2020
% Last Updated: 09/15/2020 by AA
function exhaustiveSearch(sub,e)
if ~isa(e,'double'); e = str2double(e); end
if ~isa(sub,'double'); sub = str2double(sub); end

% Settings
%--------------------------
rootDir = '/blue/camctrp/working/aprinda/freesurfer_output';
% rootDir = 'L:\working\aprinda\freesurfer_output';
I_max = 4; % Max Intensity of Precision Dose (mA)
simTag = 'FStarget';
logname = ['exhaustiveSearch_sub-' num2str(sub) '_log.txt'];
%--------------------------

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
if exist(fullfile(dirname,[baseFilename '_' elecName{elecIdx(1,e)} '_roastResult.mat']),'file') && ...
        exist(fullfile(dirname,[baseFilename '_' elecName{elecIdx(2,e)} '_roastResult.mat']),'file')
    load(fullfile(dirname,[baseFilename '_' elecName{elecIdx(1,e)} '_roastResult.mat']),'ef_all')
    A = ef_all(:);
    load(fullfile(dirname,[baseFilename '_' elecName{elecIdx(2,e)} '_roastResult.mat']),'ef_all')
    A(:,2) = ef_all(:); clear ef_all; % save RAM

    masks = load_nii(fullfile(dirname, [baseFilename '_T1orT2_masks.nii']));
    allMask = double(masks.img); clear masks;
    Nvox = numel(allMask);
    
    wnii = load_nii(fullfile(fileparts(dirname),'iUFmnipos_weights_UFmni_ideal.nii'));
    native_w = wnii.img(:)'; native_w(isnan(native_w)) = 0;
    avg = load_nii(fullfile(fileparts(dirname),'iUFmniRESP_EFavg_UFmni_defcon.nii'));
    gmean = avg.img(:); 
%     gmean = gmean./repmat(sqrt(sum([gmean(1:Nvox) gmean(Nvox+1:Nvox*2) gmean((Nvox*2)+1:end)].^2,2)),3,1); % unit vector
    gmean(isnan(gmean)) = 0; 
    sd = load_nii(fullfile(fileparts(dirname),'iUFmniRESP_EFstd_UFmni_defcon.nii'));
    gstd = sd.img(:); 
%     gstd = gstd./repmat(sqrt(sum([gstd(1:Nvox) gstd(Nvox+1:Nvox*2) gstd((Nvox*2)+1:end)].^2,2)),3,1); % unit vector
    gstd(isnan(gstd)) = 0;
    clear wnii avg sd; % save RAM
    
    % Create Explicit Masks
    brain = repmat(allMask(:) == 1 | allMask(:) == 2,3,1)';
    gmean(~brain) = 0; gstd(~brain) = 0;
%     wmask = repmat(allMask(:) == 2,3,1)'; 
    wmask =  brain & native_w > 0;
    
    % Convert Voxel Coord to Mesh Nodes
    disp('Computing the optimized electric field (this may take a while) ...');
    fJbrain = zeros(length(I),sum(wmask));
    for i = 1:length(I)
        I_opt = [-I(i);I(i)];
        fJbrain(i,:) = A(wmask,:) * I_opt; % fJbrain(i,~brain) = 0;
    end
    fJbrain(isnan(fJbrain)) = 0;
    ifJbrain = [fJbrain; -fJbrain];
    LL= exp(-sum((ifJbrain - gmean(wmask)').^2 ./ (gstd(wmask)'.^2 + 1),2)); % log-likelihood
    [maxLL,ind_LL] = max(LL,[],1);
    mA = [I; I]; Iprec = mA(ind_LL);
    
    
    %         fJbrain = fJbrain./repmat(sqrt(sum([fJbrain(1:Nvox) fJbrain(Nvox+1:Nvox*2) fJbrain((Nvox*2)+1:end)].^2,2)),3,1);% unit vector
    %         fJbrain = sqrt(sum([fJbrain(1:Nvox) fJbrain(Nvox+1:Nvox*2) fJbrain((Nvox*2)+1:end)].^2,2)); % ef_mag
        
        if ~exist(logname,'file')
            msg = 'Beginning Exhaustive Search ...';
            save(logname,'msg','-ascii')
        end
        LL = exp(-sum((fJbrain(wmask) - gmean(wmask)).^2 ./ (gstd(wmask).^2 + 1),'omitnan')); % log-likelihood
%         [~,~,mse] = lscov(fJbrain(brain),gmean(brain),(native_w(brain)/max(native_w(brain)))); % weighted least squares
        
        % Save Results
        fid = fopen(logname, 'a');
        fprintf(fid, '%s %smA %s %smA: %s\n',elecName{elecIdx(1,e)},num2str(-I(i)),elecName{elecIdx(2,e)},num2str(I(i)),num2str(LL));
        fclose(fid);
                
        % Swap polarity for inverse montage
        
        LL = exp(-sum((ifJbrain(wmask) - gmean(wmask)).^2 ./ (gstd(wmask).^2 + 1))); % log-likelihood
%         [~,~,mse] = lscov(ifJbrain(brain),gmean(brain),(native_w(brain)/max(native_w(brain))));
        
        % Save Results
        fid = fopen(logname, 'a');
        fprintf(fid, '%s %smA %s %smA: %s\n',elecName{elecIdx(1,e)},num2str(I(i)),elecName{elecIdx(2,e)},num2str(-I(i)),num2str(LL));
        fclose(fid);
        %         disp(['Electrode ' elecName{elecIdx(1,e)} ' ' num2str(mA(i)) 'mA ' elecName{elecIdx(2,e)} ' ' num2str(-mA(i)) 'mA computed: ' num2str(LL) ' log-likelihood'])
    end
    close all;
end
toc
