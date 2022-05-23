function [electrode_coord,midline_opt,center,rscalp]= getIdeal(FB,scalp_surface,dims,landmarks,capInfo,indNeed)
% [electrode_coord,center]= fitCap2individual(scalp,scalp_surface,landmarks,P2,capInfo,indNeed,isBiosemi,isEGI)
%
% Place the electrodes with pre-defined coordinates in the standard EEG
% system (e.g., 10/05, BioSemi, or EGI system).
%
% Landmarks follow the order of: nasion, inion, right, left, front neck,
% and back neck.
%
%
% (c) Yu (Andy) Huang, Parra Lab at CCNY
% yhuang16@citymail.cuny.edu
% April 2018
% August 2019 added EGI layout
% November 2020 Rewritten for tDCSLAB by Alejandro Albizu

nasion = landmarks(1,:);
inion = landmarks(2,:);
right = landmarks(3,:);
left = landmarks(4,:);

% front_neck = landmarks(5,:);
% back_neck = landmarks(6,:);

disp('measuring head size...')
L = norm(inion-nasion); % Distance between nasion and inion
line_center = (inion+nasion)/2; % Midpoint between nasion and inion

disp('wearing the cap...')
elec = capInfo{1};
centralElec = {'Oz';'POz';'Pz';'CPz';'Cz';'FCz';'Fz';'AFz';'Fpz'};
[~,indCentralElec] = ismember(centralElec,elec);

indFit = cat(1,indCentralElec,indNeed); % only fit those elec specified by users (to save time)
elec_template = cell2mat(capInfo(2:4));
elec_template = elec_template(indFit,:);
elec_template = elec_template./repmat([0.5 0.5 0.5],length(indFit),1);

theta = 23;
alpha = ((360-10*theta)/2)*(pi/180);
h = (L/2)*(1/tan(alpha));
% For the calculation of the center of electrode coordinates
% This applies to 1010 and BioSemi systems, but only roughly applies to EGI system

s = right-left; s = s/norm(s);
c = nasion-inion; c = c/norm(c);
a = cross(s,c); a = a/norm(a); % vectors to be used in affine transform
disp('adjust the cap for optimized position...this will take a while...')
factor = 1:-0.05:0.5; % Adjusting factor
F = zeros(length(factor),1);

CENTER = zeros(length(factor),3); midline = zeros(length(centralElec)+2,3,length(factor));
ELEC_COORD = zeros(size(elec_template,1),3,length(factor));
for n = 1:length(factor)
    fprintf('Iteration No. %d...\n', n);
    center = line_center + h*factor(n)*a; % Adjust the center
    CENTER(n,:) = center; % buffer
    scale = round(max(dims)/2);
    shift = center';
    
    affine = scale * [s' c' a' shift/scale;0 0 0 1/scale]; % Affine transform matrix
    
    elec_adjusted = [elec_template';ones(1,size(elec_template,1))];
    elec_adjusted(3,:) = elec_adjusted(3,:)*factor(n);
    % Adjust the z-coordinate correspondingly
    elec_transformed = affine * elec_adjusted;
    elec_transformed = elec_transformed(1:3,:)';
    % Affine transform the EasyCap coordinate to an approximate position for each electrode outside of the scalp surface
    idx = zeros(size(elec_transformed,1),1);
    [cosineAngle,indOnScalpSurf] = project2ClosestSurfacePoints(elec_transformed,scalp_surface,center);
    for i = 1:length(idx)
        %         testPts = scalp_surface(indOnScalpSurf(cosineAngle(:,i) > max(cosineAngle(:,i))*0.99993,i),:);
        testPts = scalp_surface(indOnScalpSurf(cosineAngle(:,i) > prctile(cosineAngle(:,i),99.99),i),:);
        [~,indFarthestOnTestPts] = map2Points(center,testPts,'farthest');
        idx(i) = indOnScalpSurf(indFarthestOnTestPts,i);
        % Find the only point on the outer surface of the scalp for each electrode,
        % i.e., the exact coordinates for each electrode on the scalp surface
    end
    
    elec_interp = scalp_surface(idx,:); % exact coordinates for each electrode
    ELEC_COORD(:,:,n) = elec_interp; % buffer
    
    rscalp = round(scalp_surface,1);
    
    %     center_points = [inion;elec_interp(indCentralElec,:);nasion];
    center_points = [inion;elec_interp(1:length(indCentralElec),:);nasion];
    midline(:,:,n) = center_points;
    % coordinates for electrodes on central sagittal line
    center_fit = rscalp(rscalp(:,1) == round(line_center(1),1) ... % HARDCODED
        & rscalp(:,3) >= line_center(3),:); % [centralSag*ones(length(yi),1) yi zi]; % coordinates for each point on central sagittal line
    yi = center_fit(:,2); zi = center_fit(:,3);
    
    indxsave = 1;
    distance = zeros(size(center_points,1)-1,1);
    for ii = 2:size(center_points,1)
        electemp = repmat(center_points(ii,:),[size(center_fit,1) 1]);
        [~,indx] = min(sqrt(sum((electemp - center_fit).^2,2)));
        distance(ii-1) = sum(sqrt(diff(yi(indxsave:indx)).^2+diff(zi(indxsave:indx)).^2));
        % Calculate the distance between every adjacent pair of electrodes on central sagittal line
        indxsave = indx;
    end
    F(n) = sum(abs(distance-FB/10)); % the total error compared to ideal 10-10 system, this error needs to be minimized
    % Optimize the location of the center, to make the distance between each adjacent electrode on central sagittal line equal to
    % 1/10 of the total distance from nasion to inion
end

[~,index] = min(F); midline_opt = midline(:,:,index);
% electrode_coord = ELEC_COORD(:,:,index); % exact coordinate for each electrode projected on the scalp surface
electrode_coord = ELEC_COORD(length(indCentralElec)+1:end,:,index); % exact coordinate for each electrode projected on the scalp surface
center = CENTER(index,:); % center of electrode coordinates