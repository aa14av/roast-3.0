clear

% Settings
%--------------------------
rootDir = '/blue/camctrp/working/Alejandro/StimBrain';
outDir = fullfile(rootDir,'DEFCON_II');
logDir = fullfile(pwd,'ES_DOSE');

subs = [9051 9021 9023 9031 9032 9047 9054]; %1892 101549 300100 300700 301051]; % 9054 300700
psy = uint8([32 26 19 20 20 17 15 2 4 6 10 16 18 7])'; % ACC + RT
dtype = 'DI'; % Data Type
I_max = 4; % Max Intensity of Precision Dose (mA)
uthr = 0.15; % Voxel Overlap Outlier Threshold
recipe = {'F3',-2,'F4',2}; % Original Montage
dims = [182 218 182];
range = [-0.1 0.1]; 
title_font = 12; sub_font = 10;
Nsub = 7; 

color(:,:,1) = [.635 .078 .184]; 
color(:,:,2) = [0 .447 .741]; 
color(:,:,3) = [0.466 0.674 0.188];% Red/Blue/Green

simTag = 'FStarget';
uniTag = '20220110';
%--------------------------
tic; % Start Timing

rng default; % For Reproducibility

I = 1:0.1:I_max; % HARD-CODED

switch dtype; case 'I'; nd = 1; case 'DI'; nd = 3; end % # of Dims

% Define Every Possible Electrode Combination
fid = fopen('./elec72.loc'); T = textscan(fid,'%d %f %f %s');
fclose(fid); elecName = strrep(T{4},'.','');
Nelec = length(elecName)-1; % Remove Reference (Iz)
elecIdx = uint8(nchoosek(1:Nelec,2)');
clear Nelec T; % save RAM

% Define Group Labels as above/below Median
label = int8(zeros(Nsub,1));
label(psy >= median(psy)) = 1; % Responders
label(psy < median(psy)) = -1; % Non Responders

% Load Data
load(fullfile(outDir,'interpdata_MNI_DI_VIII.mat'),'alldata');
% load(fullfile(outDir,'alldata_MNI_DI_II.mat'),'alldata')
load(fullfile(outDir,'allAM_MNI_III.mat'),'allAM'); % Tissue Segmentations
umask = repmat(sum(allAM == 1 | allAM == 2)' > Nsub*uthr,nd,1);
clear allAM; % Save RAM

% Get Responder Mean / Std
J_R = alldata(label == 1,:); J_R(J_R == 0) = NaN;
J_R(isnan(J_R)) = NaN; MU = double(mean(J_R,1,'omitnan'))';
J_NR = alldata(label == -1,:); J_NR(J_NR == 0) = NaN;
J_NR(isnan(J_NR)) = NaN; clear alldata;

% Load Individual Mask and Erode
for s = 1:length(subs)
    sub = subs(s); logname = ['exhaustiveSearch_sub-' num2str(sub) '_' dtype '_log.txt'];
    dirname = fullfile(rootDir,['FS6.0_sub-' num2str(sub) '_ses01'],simTag); % T1 dir
    masks = niftiread(fullfile(fileparts(dirname),'ROAST_output','wufT1_T1orT2_masks.nii'));
    allMask = zeros(size(masks));  bin = imfill(imerode(masks ~= 0 & ~isnan(masks),strel('cube',4)),'holes');
    allMask(bin) = masks(bin); Nvox = numel(allMask);
    brain = repmat(allMask(:) == 1 | allMask(:) == 2,nd,1); clear bin masks;
    
    data = readtable(fullfile(logDir,uniTag,logname),'ReadVariableNames',false);
    data.Var5(isnan(data.Var5)) = 0; [~,sidx] = sort(data.Var5,'descend'); 
    sprec = data(sidx,:); elecs = [sprec{1,1} sprec{1,3}];
    
    I_d = [str2double(sprec{1,2}{1}(1:regexp(sprec{1,2}{1},'mA')-1));
         str2double(sprec{1,4}{1}(1:regexp(sprec{1,4}{1},'mA')-1))];
    
    % Load Lead Fields
    A = zeros(Nvox*nd,2);
    for ee = 1:size(elecIdx,1)
        nii = niftiread(fullfile(dirname,['wufT1_' elecs{ee} '_roastResult.nii']));
        A(:,ee) = nii(:);
    end; A(~brain,:) = NaN; clear ee mag nii;
    
    % Fit Gaussian Model to model J-map of Responders (Precision Current Prediction)
    A = (A*I_d)'; % Compute EF from Lead Fields <---------------------------- HARDCODED
    A(A == 0) = NaN; A(~umask) = NaN; % Remove Non-Brain Voxels
    A(repmat(allMask(:) == 1,1,nd)) = A(repmat(allMask(:) == 1,1,nd)) .* 0.126; % Convert to J for WM
    A(repmat(allMask(:) == 2,1,nd)) = A(repmat(allMask(:) == 2,1,nd)) .* 0.276; % Convert to J gor GM
    idata = reshape(A,[Nvox nd]); clear allAM; % X = sign(X).*log(abs(X));
    
    % Impute Missing Data
    for d = 1:nd; idata(:,d) = impute_nii_old2(idata(:,d)',umask(1:Nvox)',dims); end; A = idata(:); clear idata;
    J_P(s,:) = A; clear A;
    
    e1 = load_untouch_nii(fullfile(dirname,['wufT1_' recipe{1} '_roastResult.nii'])); % HARDCODED
    B(:,1) = e1.img(:);
    e2 = load_untouch_nii(fullfile(dirname,['wufT1_' recipe{3} '_roastResult.nii'])); % HARDCODED
    B(:,2) = e2.img(:); B = B*[recipe{2};recipe{4}]; B(~brain) = NaN;
    B(isnan(B)) = NaN; B(repmat(allMask(:) == 1,nd,1)) = B(repmat(allMask(:) == 1,nd,1)) .* 0.126;
    B(repmat(allMask(:) == 2,nd,1)) = B(repmat(allMask(:) == 2,nd,1)) .* 0.276;
    J_NR(s,:) = B; clear B;
end
    
% SVM Weight Mask
load(fullfile(outDir,'a_DI_MNI.mat'),'pos_act');
wbin = pos_act ~= 0; wmag = wbin(1:Nvox); % weight mask

% Responder Mean
J_R(J_R == 0) = NaN; J_Rmag = sqrt(sum(reshape(J_R(1:Nsub,:).^2,[Nsub Nvox 3]),3));

% Fixed Dose
J_NR(J_NR == 0) = NaN; J_NRmag = sqrt(sum(reshape(J_NR(1:Nsub,:).^2,[Nsub Nvox 3]),3));

% Precision Dose
J_P(J_P == 0) = NaN; J_Pmag = sqrt(sum(reshape(J_P(1:Nsub,:).^2,[Nsub Nvox 3]),3));

J_stack = [J_P;J_NR;J_R]; %[J_P;J_NR;J_R]; 
J_stack(isnan(J_stack)) = 0; % PCA cannot handle NaN
for s = 1:Nsub
    J_all(s,:) = [mean(J_NR(s,wbin),'omitnan'),...
                  mean(J_P(s,wbin),'omitnan'),...
                  mean(J_R(s,wbin),'omitnan')]; 
    J_mag(s,:) = [mean(J_NRmag(s,wmag),'omitnan'),...
                   mean(J_Pmag(s,wmag),'omitnan'),...
                   mean(J_Rmag(s,wmag),'omitnan')];
end 
r_norm = mean(J_Rmag,1,'omitnan'); n_mag = mean(J_NRmag,1,'omitnan'); p_mag = mean(J_Pmag,1,'omitnan');
r_mean = mean(J_R,1,'omitnan'); n_mean = mean(J_NR,1,'omitnan'); p_mean = mean(J_P,1,'omitnan');
diff1 = -median(reshape(n_mean,[Nvox 3])-reshape(r_mean,[Nvox 3]),2,'omitnan'); %diff1(temp.img == 0) = NaN; % Remove Zeros for Plotting
diff2 = -median(reshape(p_mean,[Nvox 3])-reshape(r_mean,[Nvox 3]),2,'omitnan'); %diff2(temp.img == 0) = NaN; % Remove Zeros for Plotting
% diff1 = n_norm - r_norm; diff2 = p_norm - r_norm;
rsn = reshape(diff1, dims); rsp = reshape(diff2, dims);

% Create Figure
f1 = figure('units','normalized','position',[0 0 1 1]); cm = redbluecmap; % rngy=1:5:size(xi,2); rngz=slice; cm(6,:) = [NaN NaN NaN];

% Principal Compnent Analysis (takes a while)
[~, svs] = pca(J_stack(:,wbin));

% Clustering
clab = [repelem({'Optimized Doses'},size(J_P,1),1);repelem({'Non-Responder Doses'},size(J_NR,1),1);repelem({'Responder Doses'},size(J_R,1),1)];
s2=subplot(3,3,3); s2pos = s2.Position; gs1 = scatterhist(svs(:,1),svs(:,2),'Group',clab,'kernel','on','Location','Southwest','Direction','out','Color', ...
    [color(:,:,2);color(:,:,1);color(:,:,3)],'parent',f1,'LineStyle','-','LineWidth',2,'legend','off'); alpha(gs1(1),0); hold on;
u = uipanel(f1,'title','Gaussian Mixture Dose Clustering','Position',[s2pos(1:2)-0.065 s2pos(3:4)+0.065], ...
    'bordertype','none','BackgroundColor',[1 1 1]); set(gs1,'Parent',u); gax = get(gs1(1),'Children'); 
for i = 1:length(unique(clab)); gax(i).MarkerFaceColor = 'none';  gax(i).MarkerEdgeColor = 'none'; end 
GMM = fitgmdist([svs(:,1),svs(:,2)],2,'CovarianceType','diagonal','SharedCovariance',true);
gmPDF = @(x,y) arrayfun(@(x0,y0) pdf(GMM,[x0 y0]),x,y); fc = fcontour(gmPDF,[gs1(1).XLim gs1(1).YLim]); hold on;
gs2 = scatter(svs(1:7,1),svs(1:7,2),100,repmat(color(:,:,2),7,1),'filled','MarkerEdgeColor','k'); alpha(gs2,0.7);
gs3 = scatter(svs(8:14,1),svs(8:14,2),100,repmat(color(:,:,1),7,1),'filled','MarkerEdgeColor','k'); alpha(gs3,0.7);
gs4 = scatter(svs(15:21,1),svs(15:21,2),100,repmat(color(:,:,3),7,1),'filled','MarkerEdgeColor','k'); alpha(gs4,0.7);
% annotation(u,'textbox',[.42 .5 .1 .1],'String',sprintf('Average Optimization\nResponse Likelihood = %s%%',num2str(round(mean(LL),4)*100)),'FaceAlpha',0,'LineStyle','none','FontWeight','bold','FontSize',title_font)
xlabel(gs1(1),'PC1','fontname','arial','fontsize',sub_font,'fontweight','bold');  
ylabel(gs1(1),'PC2','fontname','arial','fontsize',sub_font,'fontweight','bold');
legend([gs2 gs3 gs4],{'Optimized Doses','Non-Responder Doses','Responder Doses'},'Location','Southwest')

% Bar Plot
s10 = subplot(3,3,[1 2 4 5],'Parent',f1);
for d = 1:2
    b2 = bar(d,mean(J_mag(:,d)),'barwidth',0.4,'FaceColor','none','LineWidth',3,'FaceAlpha',0,'EdgeAlpha',.6);
    hold on; b2.EdgeColor = color(:,:,d);
end; hold on;
dp1 = distributionPlot(J_mag(:,1:2),'distwidth',0.2,'histOpt',0,'globalNorm',2,'color', ...
    {color(:,:,1),color(:,:,2)},'showMM',0,'histOri','right','xValues',[1.5 2.5]); alpha(dp1{3},0.4)
plotSpread(J_mag(:,1:2), ...
    'distributionColors','k','distributionMarkers',{'o', 's'},'showMM',0, ...
    'xNames',{'Fixed Doses','Optimized Doses'},'yLabel','Average Current Density (Am^{-2})');
ylim([0 max(J_mag(:,3)+std(J_mag(:,3)))])
hold on; set(gca,'FontSize',title_font,'Fontname','arial','FontWeight','bold'); title('Dosing Variability');
gmu = mean(J_mag(:,3)); gsd = std(J_mag(:,3)); %gmu = mean(mean(J_Rmag(:,all(J_Rmag~=0)),2,'omitnan')); gsd = std(mean(J_Rmag(:,all(J_Rmag~=0)),2,'omitnan'));
l1= line([-10 10],[gmu gmu],'Color','k','LineStyle',':','LineWidth',2); hold on;
cip = gmu + (1.96*(gsd)/sqrt(Nsub)); cin = gmu - (1.96*(gsd)/sqrt(Nsub));
pat2=patch([-10 -10 10 10],[cin cip cip cin],'k');
pat2.FaceAlpha = 0.1; pat2.EdgeAlpha = 0;
plot([1 2],[cip + 0.0005 cip + 0.0005],'-k','LineWidth',2) % Effect Size
scatter([1.4 2.4],[mean(J_mag(:,1)) mean(J_mag(:,2))],100,'k','filled') % mean
plot([1.4 1.4],[mean(J_mag(:,1))-std(J_mag(:,1)) mean(J_mag(:,1))+std(J_mag(:,1))],'-k','LineWidth',2) % 95%CI
plot([2.4 2.4],[mean(J_mag(:,2))-std(J_mag(:,2)) mean(J_mag(:,2))+std(J_mag(:,2))],'-k','LineWidth',2) % 95%CI
stats = mes(J_mag(:,2),J_mag(:,1),'hedgesg','isDep',1,'nBoot',1000);
annotation('textbox',[.31 .8 .1 .1],'String',sprintf('Hedges'' g = %.02f',round(stats.hedgesg,2)),'FaceAlpha', ... 
    0,'LineStyle','none','FontWeight','bold','FontSize',title_font)
legend([l1 pat2],{'Responder Mean (\mu)','Responder 95% CI'})

% Optimized Dose
s1=subplot(3,3,8); imagesc(rsp(:,:,slice),[-0.025 0.025]); axis off; axis tight; title('Optimized - Responder Mean','fontname','arial','fontsize',title_font,'fontweight','bold'); 
colormap(s1,cm); alpha(s1,0.5); cb1 = colorbar;
ylabel(cb1,'Difference (Am^{-2})','fontname','arial','fontsize',sub_font,'fontweight','bold'); hold on;
% q1 = quiver(yi(rngx,rngy,rngz),xi(rngx,rngy,rngz),...
%         rspad(rngx,rngy,rngz,2),rspad(rngx,rngy,rngz,1),2,'color','k'); hold off;

% Fixed Dose
s3=subplot(3,3,7); imagesc(rsn(:,:,slice),[-0.025 0.025]); axis off; axis tight; title('Fixed - Responder Mean','fontname','arial','fontsize',title_font,'fontweight','bold'); 
colormap(s3,cm); alpha(s3,0.5); cb2 = colorbar; ylabel(cb2,'Difference (Am^{-2})','fontname','arial','fontsize',sub_font,'fontweight','bold'); hold on;
% q2 = quiver(yi(rngx,rngy,rngz),xi(rngx,rngy,rngz),...
%         rsgd(rngx,rngy,rngz,2),rsgd(rngx,rngy,rngz,1),2,'color','k'); hold off;
    
% Scatter
s9=subplot(3,3,9); scatter(p_mean(aparc),r_mean(aparc),8,'filled'); xlim(range); ylim(range); ls1 = lsline; ls1.LineWidth = 2;
xlabel('Optimized Mean','fontname','arial','fontsize',sub_font,'fontweight','bold'); ylabel('Responder Mean','fontname','arial','fontsize',sub_font,'fontweight','bold');
par = corr(p_mean(aparc)',r_mean(aparc)','rows','complete'); annotation('textbox',[.7 .22 .1 .1],'String',sprintf('R^2 = %.3f',par^2),'FaceAlpha',0,'LineStyle','none','FontWeight','bold','FontSize',sub_font)
title('Optimized vs Responder Mean','fontname','arial','fontsize',title_font,'fontweight','bold'); 

% Histogram
s4=subplot(3,3,6);

ps = abs(p_mag(wmag)); rs = abs(n_mag(wmag)); gs = abs(r_norm(wmag));
d1 = ps; d1d = prctile(d1,95); h1 = histfit(d1,400,'kernel'); alpha(h1,0.01); hold(s4,'on'); %aparc' & p_mean ~= 0 & ~isnan(p_mean)
d2 = rs; d2d = prctile(d2,95); h2 = histfit(d2,400,'kernel'); alpha(h2,0.01); %aparc' & r_mean ~= 0 & ~isnan(r_mean)
d3 = gs; d3d = prctile(d3,95); h3 = histfit(d3,400,'kernel'); alpha(h3,0.01); %aparc' & g_mean ~= 0 & ~isnan(g_mean)

h3(1).FaceColor = color(:,:,3); h3(2).Color = h3(1).FaceColor;
[~,lid3] = min(pdist2(h3(2).XData',d3d));
sc3 = scatter(d3d,h3(2).YData(lid3),100,'MarkerFaceColor',color(:,:,3),'MarkerEdgeColor','k');

h1(1).FaceColor = color(:,:,2); h1(2).Color = h1(1).FaceColor;
[~,lid1] = min(pdist2(h1(2).XData',d1d));
sc1 = scatter(d1d,h1(2).YData(lid1),100,'MarkerFaceColor',color(:,:,2),'MarkerEdgeColor','k');

h2(1).FaceColor = color(:,:,1); h2(2).Color = h2(1).FaceColor;
[~,lid2] = min(pdist2(h2(2).XData',d2d));
sc2 = scatter(d2d,h2(2).YData(lid2),100,'MarkerFaceColor',color(:,:,1),'MarkerEdgeColor','k');

% for s = 1:size(J_R,1)
%     d3 = J_stack(s+14,:);
%     h4 = histfit(d3(~isnan(d3)),400,'kernel'); alpha(h4,0); %h3(1).FaceColor = color(:,:,3); h3(2).Color = h3(1).FaceColor;
%     h4(1).FaceColor = color(:,:,3); h4(2).Color = h4(1).FaceColor;
% end

yt = get(s4,'YTick'); ylim(s4,[0 max(yt)]); set(s4,'YTick',yt,'YTickLabels',round(yt/numel(p_mag(aparc(1:Nvox))),4)*100);
xlim([0 0.1]); xlabel('Current Density (Am^{-2})','fontname','arial','fontsize',sub_font,'fontweight','bold'); ylabel('Percent of Voxels','fontname','arial','fontsize',sub_font,'fontweight','bold')
title(s4,'Electric Field Distribution','fontname','arial','fontsize',title_font,'fontweight','bold');
legend([h1(2),h2(2),h3(2),sc1],{'Optimized Mean','Non-Responder Mean','Responder Mean'},'fontname','arial','fontsize',sub_font,'fontweight','bold')

% Difference
% s4=subplot(2,2,4); imagesc(diff(:,:,slice),[-0.1 0.1]); axis off; axis tight; title(['x - ' char(251)]); colormap(s4,'gray'); colorbar
% title(s4,'Current Intensity Difference','fontname','arial','fontsize',12,'fontweight','bold');
% saveas(gcf,fullfile(rootDir,'ALLDAYS','PRECISE_avg_summary_II.svg'))
set(f1,'Units','inches');
screenposition = get(gcf,'Position');
set(f1,...
    'PaperPosition',[0 0 screenposition(3:4)],...
    'PaperSize',screenposition(3:4));
print(fullfile(outDir,'PRECISE_avg_summary_II'),'-dpdf','-painters') 
print(fullfile(outDir,'PRECISE_avg_summary_II'),'-dsvg','-painters') 
print(fullfile(outDir,'PRECISE_avg_summary_II'),'-depsc','-painters')
print('-r600',fullfile(outDir,'PRECISE_avg_summary_II'),'-dpng','-painters')
% close all
