function Jbrain = convertMesh2Vox(subDir,elecName,hdrInfo,allMask,nodeV,A_all,I_opt)
% if ~exist(fullfile(subDir,['T1_' elecName{I_opt > 0} '_leadField.mat']),'file')
    [xi,yi,zi] = ndgrid(1:hdrInfo.dim(1),1:hdrInfo.dim(2),1:hdrInfo.dim(3));
    ef_all = zeros([hdrInfo.dim(1:3) 3]); Jroast = zeros(hdrInfo.dim(1:3));
    isNaNinA = isnan(sum(sum(A_all,3),2)); % handle NaN properly
    xopt = zeros(sum(~isNaNinA),4);
    xopt(:,1) = find(~isNaNinA);
    for i=1:size(A_all,2), xopt(:,i+1) = squeeze(A_all(~isNaNinA,i,:))*I_opt; end
    
    F = TriScatteredInterp(nodeV(~isNaNinA,1:3), xopt(:,2));
    ef_all(:,:,:,1) = F(xi,yi,zi);
    F = TriScatteredInterp(nodeV(~isNaNinA,1:3), xopt(:,3));
    ef_all(:,:,:,2) = F(xi,yi,zi);
    F = TriScatteredInterp(nodeV(~isNaNinA,1:3), xopt(:,4));
    ef_all(:,:,:,3) = F(xi,yi,zi);
    ef_mag = sqrt(sum(ef_all.^2,4));
    
    Jroast(allMask == 1) = ef_mag(allMask == 1); % * 0.126; % WM
    Jroast(allMask == 2) = ef_mag(allMask == 2); % * 0.276; % GM
    Jroast(allMask == 3) = ef_mag(allMask == 3); % * 1.65; % CSF
    Jroast(allMask == 4) = ef_mag(allMask == 4); % * 0.01; % BONE
    Jroast(allMask == 5) = ef_mag(allMask == 5); % * 0.465; % SKIN
    Jroast(allMask == 6) = ef_mag(allMask == 6); % * 2.5e-14; % AIR
    brain = (allMask==1 | allMask==2);
    Jbrain = zeros(hdrInfo.dim(1:3)); Jbrain(brain) = Jroast(brain);
    save(fullfile(subDir,['T1_' elecName{I_opt > 0} '_leadField.mat']),'Jroast','Jbrain')
% else
%     load(fullfile(subDir,['T1_' elecName{I_opt > 0} '_leadField.mat']),'Jbrain')
% end