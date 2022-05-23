function [electrode_coord,center,elecName]= placeActualLandmarks(sub_info,filenames,rs,rootDir)
% (c) Yu (Andy) Huang, Parra Lab at CCNY
% yhuang16@citymail.cuny.edu
% April 2018

%----------------------------------------
% Adapted By Alejandro Albizu
% Center for Cognitive Aging and Memory
% University of Florida
% 11/09/2020
%----------------------------------------
% Last Updated: 05/10/2021 by AA

% Setup Variables
s = sub_info(1); d = sub_info(2); n = sub_info(3); 
T1 = rs.output{1}; T2 = rs.output{2};
[dirname,baseFilename] = fileparts(T1);
if isempty(T2); baseFilenameRasRSPD = [baseFilename '_T1orT2'];
else; baseFilenameRasRSPD = [baseFilename '_T1andT2'];
end

% Load Data
load(fullfile(rootDir,'cap1005FullWithExtra.mat'),'capInfo');
if exist(fullfile(rootDir,'Data_Output',strcat('sub-',num2str(s)),strcat('sub',num2str(s),'_day',num2str(d),'_scan',num2str(n),'_dist.mat')),'file')
    load(fullfile(rootDir,'Data_Output',strcat('sub-',num2str(s)),strcat('sub',num2str(s),'_day',num2str(d),'_scan',num2str(n),'_dist.mat')),'head3D','midline','actual','elecName');
else; error('Please run ''Placement Accuracy'' module first');
end

% Rotate to RAS Orientation
disp('converting head to RAS orientation...')
cap = [capInfo{2:4}]; [cperm,~,cisFlipInner,~] = how2getRAS( ...
    [cap(ismember(capInfo{1},{'Nz'}),:); cap(ismember(capInfo{1},{'Iz'}),:); ...
    cap(ismember(capInfo{1},{'RPA'}),:); cap(ismember(capInfo{1},{'LPA'}),:); ...
    [cap(ismember(capInfo{1},{'Nz'}),1) cap(ismember(capInfo{1},{'Nz'}),2) min(cap(:,3))]]);
rhead = head3D.vertices; rcap = changeOrientationPointCloud(cap,cperm,cisFlipInner,[1 1 1]);

% Extract Electrode Surface 
disp('extracting electrode surface...')
[elec,~] = getElecOri(rootDir,filenames,sub_info,rs);
rpad = [elec{1};elec{2}]; % changeOrientationPointCloud([elec{1};elec{2}],perm,isFlipInner,dims);

% Center Vector at Head Fpz Centroid
rctr = [midline(10,1) mean(rhead(:,2)) midline(10,3)]; % Fpz Centroid
ulm = (actual - rctr); upad = rpad - rctr; %./sum(abs(rlm - rctr),2); 

% Create Cap Sphere
rrad = norm(rcap(ismember(capInfo{1},'Fpz'),:) - rcap(ismember(capInfo{1},'Oz'),:),2)/2;
[sx,sy,sz]=sphere; rcctr = table2array(regionprops3(cat(3,sx,sy,sz),'centroid')); % [rcap(ismember(capInfo{1},'Fpz'),1) mean(rcap(:,2)) rcap(ismember(capInfo{1},'Fpz'),3)];

% Scale Unit Vector to Cap
disp('normalizing extracted electrode locations...')
P = ulm - rcctr; Q = (rrad./sqrt(sum(P.^2,2))).*P; slm = Q; % + rcctr;
spad = (rrad./sqrt(sum((upad - rcctr).^2,2))).*(upad - rcctr); 

% Add Actual Electrodes to Cap
for e = 1:length(elecName)
    capInfo{2}(ismember(capInfo{1},elecName{e})) = slm(e,1);
    capInfo{3}(ismember(capInfo{1},elecName{e})) = slm(e,2);
    capInfo{4}(ismember(capInfo{1},elecName{e})) = slm(e,3);
end

% Load the scalp mask; template is used for saving the results with the same header info as the input
template = load_untouch_nii([dirname filesep baseFilenameRasRSPD '_masks.nii']);
pixdim = template.hdr.dime.pixdim(2:4);
dim = size(template.img); scalp = template.img==5; 
v2w = [template.hdr.hist.srow_x;template.hdr.hist.srow_y;template.hdr.hist.srow_z;0 0 0 1];
hdrInfo = struct('pixdim',pixdim,'dim',dim,'v2w',v2w);
if ~exist([dirname filesep baseFilename '_header.mat'],'file')
    save([dirname filesep baseFilename '_header.mat'],'hdrInfo');
end
[nx,ny,nz]=ndgrid(1:dim(1),1:dim(2),1:dim(3));

% Extract Scalp Boundary
disp('getting scalp surface...')
pad_coor = [nx(:) ny(:) nz(:)]; scoord = pad_coor(scalp(:),:); k = boundary(scoord,1);
scalp_surface = scoord(k,:); landmarks = getLandmarks(T1,T2);
nasion = landmarks(1,:); inion = landmarks(2,:); right = landmarks(3,:); left = landmarks(4,:);

% Measure Head
disp('measuring head size...')
L = norm(inion-nasion); % Distance between nasion and inion
line_center = (inion+nasion)/2; % Midpoint between nasion and inion

% Make sure the central sagittal slice can be closed completely in order to detect the edge correctly
centralSag = round(line_center(1));
img_c = squeeze(scalp(centralSag,:,:)); % The central sagittal slice
indc = find(sum(img_c)>0,1,'first');
indr1 = find(img_c(:,indc)>0,1,'first');
indr2 = find(img_c(:,indc)>0,1,'last');
img_c(indr1:indr2,indc) = 255;

% Get the edge of the central sagittal slice
se = 0; isFilled = 0;
[ytemp,ztemp] = ind2sub(size(img_c),find(img_c));
centroid = round(mean([ytemp ztemp]));
while ~all(isFilled)
    se = se+8;
    im_test = imfill(imclose(img_c,ones(se,se)),'holes');
    isFilled = im_test(centroid(1),centroid(2));
end
bw_c = edge(imopen(im_test,ones(3,3)));
[r_c,c_c] = find(bw_c==1);

% Preparation for the calculation of the distance between nasion and inion
% along scalp surface using Natural Cubic Spline
indxinion = find(inion(3)==c_c,1,'first');
indxnasion = find(nasion(3)==c_c,1,'last');
[~,I] = max(c_c);
temp_right_up = find((c_c>=c_c(indxinion))&(r_c<(r_c(I))));
[~,I_up] = sort(c_c(temp_right_up));
temp_right_up = temp_right_up(I_up);
temp_right_down = find((c_c>=c_c(indxnasion))&(r_c>=(r_c(I))));
[~,I_down] = sort(c_c(temp_right_down),'descend');
temp_right_down = temp_right_down(I_down);
index = [temp_right_up; temp_right_down];

% Approximation of 2-D Data by Natural Cubic Spline
[bx,by,finalbreaks]=ncs2dapprox(r_c(index),c_c(index));
% http://www.mathworks.co.jp/matlabcentral/fileexchange/7617

t = finalbreaks'; pp1= spline(t,[bx,by]');
range = linspace(1,finalbreaks(end),finalbreaks(end));
yizi = ppval(pp1,range); yi=yizi(1,:)'; zi=yizi(2,:)';

% Calculate the distance between nasion and inion along the scalp 
if exist(fullfile(rootDir,'Measurements',strcat('sub',num2str(s)),strcat('sub',num2str(s),'_Measurements.mat')),'file')
    load(fullfile(rootDir,'Measurements',strcat('sub',num2str(s)),strcat('sub',num2str(s),'_Measurements.mat')),'Measurements');
    distance_all = Measurements.FB; disp('loading head measurements...')
else; disp('measuring head..'); distance_all = sum(sqrt(diff(yi).^2+diff(zi).^2));
end

centralElec = {'Oz';'POz';'Pz';'CPz';'Cz';'FCz';'Fz';'AFz';'Fpz'};
[~,indCentralElec] = ismember(centralElec,capInfo{1});

indNeed = find(ismember(capInfo{1},elecName));
indFit = cat(1,indCentralElec,indNeed); % only fit those elec specified by users (to save time)
elec_template = cell2mat(capInfo(2:4)); elec_template = elec_template(indFit,:);
elec_template = elec_template./repmat([0.5 0.5 0.5],length(indFit),1);
pads = spad./repmat([0.5 0.5 0.5],length(spad),1);

theta = 23; alpha = ((360-10*theta)/2)*(pi/180); h = (L/2)*(1/tan(alpha));
% For the calculation of the center of electrode coordinates
% This applies to 1010 and BioSemi systems, but only roughly applies to EGI system

s = right-left; s = s/norm(s);
c = nasion-inion; c = c/norm(c);
a = cross(s,c); a = a/norm(a); % vectors to be used in affine transform

disp('optimizing electrode fit... this will take a while...')
factor = 1:-0.05:0.5; % Adjusting factor
F = zeros(length(factor),1);

CENTER = zeros(length(factor),3);
ELEC_COORD = zeros(size(elec_template,1),3,length(factor));
for n = 1:length(factor)
    center = line_center + h*factor(n)*a; % Adjust the center
    CENTER(n,:) = center; % buffer
    scale = round(max(size(scalp))/2);
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
        testPts = scalp_surface(indOnScalpSurf(cosineAngle(:,i) > prctile(cosineAngle(:,i),99.99),i),:);
        [~,indFarthestOnTestPts] = map2Points(center,testPts,'farthest');
        idx(i) = indOnScalpSurf(indFarthestOnTestPts,i);
        % Find the only point on the outer surface of the scalp for each electrode,
        % i.e., the exact coordinates for each electrode on the scalp surface
    end
    elec_interp = scalp_surface(idx,:); % exact coordinates for each electrode
    ELEC_COORD(:,:,n) = elec_interp; % buffer
    
    pads_adjusted = [pads';ones(1,size(pads,1))];
    pads_adjusted(3,:) = pads_adjusted(3,:)*factor(n);
    pads_transformed = affine * pads_adjusted;
    pads_transformed = pads_transformed(1:3,:)';
    
    vecP = pads_transformed - repmat(center,size(pads_transformed,1),1);
    % vector connecting center to the point that is to be projected onto a surface
    normVecP = repmat(sqrt(sum(vecP.^2,2)),1,size(vecP,2));
    vecP = vecP./normVecP;
    vecP = single(vecP);
    
    vecS = scalp_surface - repmat(center,size(scalp_surface,1),1);
    % vectors connecting center to each point on the surface
    normVecS = repmat(sqrt(sum(vecS.^2,2)),1,size(vecS,2));
    vecS = vecS./normVecS;
    vecS = single(vecS);
    [~,idx2] = max(vecS*vecP',[],1);
    pad_interp(:,:,n) = scalp_surface(idx2,:);
    
    center_points = [inion;elec_interp(1:length(indCentralElec),:);nasion];
    % coordinates for electrodes on central sagittal line
    center_fit = [centralSag*ones(length(yi),1) yi zi]; % coordinates for each point on central sagittal line
    
    indxsave = 1;
    distance = zeros(size(center_points,1)-1,1);
    for ii = 2:size(center_points,1)
        electemp = repmat(center_points(ii,:),[size(center_fit,1) 1]);
        [~,indx] = min(sqrt(sum((electemp - center_fit).^2,2)));
        distance(ii-1) = sum(sqrt(diff(yi(indxsave:indx)).^2+diff(zi(indxsave:indx)).^2));
        % Calculate the distance between every adjacent pair of electrodes on central sagittal line
        indxsave = indx;
    end
    F(n) = sum(abs(distance-distance_all/10)); % the total error compared to ideal 10-10 system, this error needs to be minimized
    % Optimize the location of the center, to make the distance between each adjacent electrode on central sagittal line equal to
    % 1/10 of the total distance from nasion to inion
end
[~,index] = min(F);
electrode_coord = ELEC_COORD(length(indCentralElec)+1:end,:,index); % exact coordinate for each electrode projected on the scalp surface
[~,ro] = sort(capInfo{1}(ismember(capInfo{1},elecName))); elecName = flip(sort(elecName));
electrode_coord = electrode_coord(flip(ro),:); center = CENTER(index,:); % center of electrode coordinates

% Model Pad HARDCODED FOR 2 ELECTRODES
pad_coord = pad_interp(:,:,index); e1_mask = zeros(dim); e2_mask = zeros(dim); 
e1_mask(sub2ind(dim,pad_coord(1:length(elec{1}),1),pad_coord(1:length(elec{1}),2),pad_coord(1:length(elec{1}),3))) = 1; % Elec 1
e2_mask(sub2ind(dim,pad_coord(length(elec{1})+1:end,1),pad_coord(length(elec{1})+1:end,2),pad_coord(length(elec{1})+1:end,3))) = 2; % Elec 2
for i = 1:size(e1_mask,3); e1_mask(:,:,i) = bwconvhull(e1_mask(:,:,i),'union'); end
for i = 1:size(e2_mask,3); e2_mask(:,:,i) = bwconvhull(e2_mask(:,:,i),'union'); end
 
% Clear Scalp
[scalp_clean,scalpFilled] = cleanScalp(scalp,scalp_surface);
temp1 = scalpFilled(:,:,[1 end]); temp2 = scalpFilled(:,[1 end],:); temp3 = scalpFilled([1 end],:,:);
if any([temp1(:);temp2(:);temp3(:)])
    warning('Scalp touches image boundary. Electrodes may go out of image boundary. ROAST can continue but results may not be accurate. It is recommended that you expand the input MRI by specifying the ''zeroPadding'' option.');
end    
scalpCleanSurf = mask2EdgePointCloud(scalp_clean,'erode',ones(3,3,3)); res= mean(pixdim);

disp('placing electrodes...')
padH = zeros(length(rs.usedElec),1); 
for i=1:length(rs.usedElec)
    elecDims = str2num(rs.output{3}{i,3});
    if strcmpi(rs.output{3}{i,2},'pad')
        padH(i) = elecDims(3)/res;
    end
end
% pad_mask = zeros(dim); pad_mask(e1_mask ~= 0 | e2_mask ~= 0) = 1;
% pad_mask = imdilate(pad_mask,ones(elecDims(3)*2)); dpad = coord(pad_mask == 1,:); 

indPad = find(padH>0);
if ~isempty(indPad)
    ind2allPH = zeros(length(indPad),1);
    [allPH,~,ind2allPH(indPad)] = unique(padH(indPad));
    gel_layer = cell(length(allPH),1); elec_layer = cell(length(allPH),1);
    for i=1:length(allPH)
        strel = ones(round(allPH(i)),round(allPH(i)),round(allPH(i))); % this is not perfect yet
        [gel_layer{i},scalpDilated] = mask2EdgePointCloud(scalpFilled,'dilate',strel);
        elec_layer{i} = mask2EdgePointCloud(scalpDilated,'dilate',strel);
        % Get the layer of electrode/gel to intersect with placed pad
    end
end

figName = 'Electrode placement in Simulation: tDCSLAB';
figure('Name',[figName '. Move your mouse to rotate.'],'NumberTitle','off');
set(gcf,'color','w');
plot3(scalpCleanSurf(:,1),scalpCleanSurf(:,2),scalpCleanSurf(:,3),'y.');
hold on;

rotation = zeros(size(rs.usedElec,1),1); butt = 1;
elec_C = cell(size(rs.usedElec,1),1); gel_C = cell(size(rs.usedElec,1),1);
% buffer for coordinates of each electrode and gel point
for i = 1:length(rs.usedElec)   
    elecRange = pad_coord(butt:(length(elec{i})+butt)-1,:);
    [U,D] = eig(cov(elecRange)); [~,ind] = min(diag(D));
    nv = U(:,ind)'; normal = nv/norm(nv); % Local normal for each electrode
    
    lenTry=1;
    testPointIn = round(electrode_coord(strcmp(rs.usedElec{i},elecName),:) - lenTry*normal);
    testPointOut = round(electrode_coord(strcmp(rs.usedElec{i},elecName),:) + lenTry*normal);
    while all(min([testPointIn;testPointOut])>0) && all(max([testPointIn;testPointOut])<=size(scalpFilled)) && ...
            ~xor(scalpFilled(testPointIn(1),testPointIn(2),testPointIn(3)),scalpFilled(testPointOut(1),testPointOut(2),testPointOut(3)))
        lenTry = lenTry+1;
        testPointIn = round(electrode_coord(strcmp(rs.usedElec{i},elecName),:) - lenTry*normal);
        testPointOut = round(electrode_coord(strcmp(rs.usedElec{i},elecName),:) + lenTry*normal);
    end
    if all(min([testPointIn;testPointOut])>0) && all(max([testPointIn;testPointOut])<=size(scalpFilled)) && ...
            scalpFilled(testPointIn(1),testPointIn(2),testPointIn(3))==0
        normal = -normal; 
    end % make sure the normal is pointing out
    
    switch lower(rs.output{3}{i,2}) 
        case 'pad'
            fprintf('placing pad at %s (%d out of %d) in the 1020 layout...\n',rs.usedElec{i},i,size(rs.usedElec,1));

            pad_length = elecDims(1)/res;
            pad_width = elecDims(2)/res;
            
            dimTry = mean([pad_length pad_width]);
            % bigger electrode needs bigger dimTry (needs to establish a better relation)
            
            % Determine Electrode Orientation
            %             [V,~] = eig(cov(pad_coord(butt:(length(elec{i})+butt)-1,:))); butt = butt + length(elec{i});
            %             plot3(scalpCleanSurf(:,1),scalpCleanSurf(:,2),scalpCleanSurf(:,3),'y.'); hold on;
            %             plot3(elecRange(:,1),elecRange(:,2),elecRange(:,3),'b.'); hold on
%             quiver3(electrode_coord(strcmp(rs.usedElec{i},elecName),1),electrode_coord(strcmp(rs.usedElec{i},elecName),2),electrode_coord(strcmp(rs.usedElec{i},elecName),3),normal(:,1)',normal(:,2)',normal(:,3)',0,'LineWidth',5);
            [~,~,U] = svd(elecRange); % EIGENVECTORS TO RIGHT SINGULAR VALUES
            quiver3(repmat(electrode_coord(strcmp(rs.usedElec{i},elecName),1),3,1),repmat(electrode_coord(strcmp(rs.usedElec{i},elecName),2),3,1),repmat(electrode_coord(strcmp(rs.usedElec{i},elecName),3),3,1),U(1,:)'*10,U(2,:)'*10,U(3,:)'*10,0,'LineWidth',5);
            if ischar(rs.output{3}{i,4})
                switch lower(rs.output{3}{i,4})
                    case 'lr'; iOri = [1 0 0]; % [~,Sind] = min(abs(U(1,:)));
                    case 'ap'; iOri = [0 1 0]; % [~,Sind] = min(abs(U(2,:)));
                    case 'si'; iOri = [0 0 1]; % [~,Sind] = min(abs(U(3,:)));
                end
            end
            
            % Generate Electrode
            iOriS = cross(normal,iOri); iOriS = iOriS/norm(iOriS);
            iOriL = cross(normal,iOriS); iOriL = iOriL/norm(iOriL);
            
%             padOriLong = zeros(3); padOriShort = zeros(3);
%             for o = 1:3
%                 padOriLong(o,:) = cross(normal,V(:,o)'); padOriLong(o,:) = padOriLong(o,:)/norm(padOriLong(o,:)); 
%                 padOriShort(o,:) = cross(normal,padOriLong(o,:)); padOriShort(o,:) = padOriShort(o,:)/norm(padOriShort(o,:));
%             end
%             Lind = knnsearch(U,iOriL,'k',1); longAxis = U(:,Lind);
            Sind = knnsearch(U,iOriS,'k',1); shortAxis = U(:,Sind)';

%             shortAxis = U(:,Sind)'; % for SVD
            longAxis = cross(normal,shortAxis); longAxis = longAxis/norm(longAxis); % make long orthoganal to short
            
            cosT = max(min(dot(longAxis,iOriL)/(norm(longAxis)*norm(iOriL)),1),-1);
            rotation(i) = real(acosd(cosT)); den = 2; % 2 points per pixel
            corner = electrode_coord(strcmp(rs.usedElec{i},elecName),:) - longAxis*pad_length/2 - shortAxis*pad_width/2 - normal*dimTry/2;
            
            numSamp_shortAxis = round(pad_width*den+1);
            numSamp_normalAxis = round(dimTry*den+1);
            
            pad_coor = cell(numSamp_shortAxis*numSamp_normalAxis,1); repr=1;
            for j = 1:numSamp_shortAxis
                for k = 1:numSamp_normalAxis
                    start = corner + 1/den*shortAxis*(j-1) + 1/den*normal*(k-1);
                    pad_coor{repr} = drawLine(start,longAxis,pad_length,den);
                    repr = repr+1;
                end
            end
            pad_coor = cell2mat(pad_coor);
            pad_coor = unique(round(pad_coor),'rows'); % clean-up of the coordinates
            
            % Separate Gel from Electrode
            gel_coor = intersect(pad_coor,gel_layer{ind2allPH(i)},'rows');
            elec_coor = intersect(pad_coor,elec_layer{ind2allPH(i)},'rows');

            plot3(elec_coor(:,1),elec_coor(:,2),elec_coor(:,3),'.b');
            plot3(gel_coor(:,1),gel_coor(:,2),gel_coor(:,3),'.m');

            gel_C{i} = gel_coor; elec_C{i} = elec_coor; % buffer for coordinates of each electrode and gel point
            butt = butt + length(elec{i});
        case 'disc'
            fprintf('placing disc at %s (%d out of %d) in the 1020 layout...\n',rs.usedElec{i},i,size(rs.usedElec,1));
            
            disc_radius = elecDims(1)/res;
            disc_height = elecDims(2)/res;
            
            dimTry = disc_radius;
            
            gel_out = electrode_coord(strcmp(rs.usedElec{i},elecName),:) +  2*disc_height*normal;
            electrode = gel_out + disc_height*normal;
            gel_in = gel_out - dimTry*normal; % coordinates of the boundaries of gel and electrode
            
            den = 2; % 2 points per pixel
            gel_coor = drawCylinder(0,disc_radius,gel_in,gel_out,den);
            elec_coor = drawCylinder(0,disc_radius,gel_out,electrode,den);
            % Use cylinders to model electrodes and gel, and calculate the coordinates of the points that make up the cylinder
            
            gel_coor = unique(round(gel_coor),'rows');
            elec_coor = unique(round(elec_coor),'rows'); % clean-up of the coordinates

            plot3(elec_coor(:,1),elec_coor(:,2),elec_coor(:,3),'.b');
            plot3(gel_coor(:,1),gel_coor(:,2),gel_coor(:,3),'.m');
            
            gel_C{i} = gel_coor; elec_C{i} = elec_coor; % buffer for coordinates of each electrode and gel point
            
        case 'ring'
            fprintf('placing ring at %s (%d out of %d) in the 1020 layout...\n',rs.usedElec{i},i,size(rs.usedElec,1));
            
            ring_radiusIn = elecDims(1)/res;
            ring_radiusOut = elecDims(2)/res;
            ring_height = elecDims(3)/res;
            
            dimTry = mean([ring_radiusOut ring_radiusIn]);
            
            gel_out = electrode_coord(strcmp(rs.usedElec{i},elecName),:) +  2*ring_height*normal;
            electrode = gel_out + ring_height*normal;
            gel_in = gel_out - dimTry*normal; % coordinates of the boundaries of gel and electrode
            
            den = 2; % 2 points per pixel
            gel_coor = drawCylinder(ring_radiusIn,ring_radiusOut,gel_in,gel_out,den);
            elec_coor = drawCylinder(ring_radiusIn,ring_radiusOut,gel_out,electrode,den);
            % Use cylinders to model electrodes and gel, and calculate the coordinates of the points that make up the cylinder
            
            gel_coor = unique(round(gel_coor),'rows');
            elec_coor = unique(round(elec_coor),'rows'); % clean-up of the coordinates
            
            plot3(elec_coor(:,1),elec_coor(:,2),elec_coor(:,3),'.b');
            plot3(gel_coor(:,1),gel_coor(:,2),gel_coor(:,3),'.m');
            
            gel_C{i} = gel_coor; elec_C{i} = elec_coor; % buffer for coordinates of each electrode and gel point
    end
end
xlabel('x');ylabel('y');zlabel('z'); axis equal; title('Reconstructed Electrode Surface');
hold off; rotate3d on; % Place electrodes and visualize the results

volume_elec_C = generateElecMask(elec_C,size(scalp),rs.usedElec,1);
volume_gel_C = generateElecMask(gel_C,size(scalp),rs.usedElec,0);

disp('final clean-up...')
volume_elec = volume_elec_C>0;
volume_gel = volume_gel_C>0;
volume_gel = xor(volume_gel,volume_gel & volume_elec); % remove the gel that overlaps with the electrode
for i=1:6
    volume_tissue = template.img==i;
    volume_gel = xor(volume_gel,volume_gel & volume_tissue);
end % remove the gel that goes into other tissue masks

disp('saving the results...')
template.fileprefix = [dirname filesep baseFilename '_tDCSLAB_mask_elec'];
template.hdr.hist.descrip = 'electrode mask';
template.img = uint8(volume_elec_C.*volume_elec);
save_untouch_nii(template,[dirname filesep baseFilename '_tDCSLAB_mask_elec.nii']);
template.fileprefix = [dirname filesep baseFilename '_tDCSLAB_mask_gel'];
template.hdr.hist.descrip = 'gel mask';
template.img = uint8(volume_gel_C.*volume_gel);
save_untouch_nii(template,[dirname filesep baseFilename '_tDCSLAB_mask_gel.nii']);