%% Skip Seg
function hiperskipseg(s,day,subj,varargin)

% Settings
%=======================================
recipe = {'F3', -2, 'F4', 2};
elecType = {'pad','pad'};
elecSize = {[70 50 3],[70 50 3]};
elecOri = 'lr';
capType = '1020';
paddingAmt = 0;
doResamp = 0;
simTag = 'MVPA';
%=======================================

indArg = 1;
while indArg <= length(varargin)
    switch lower(varargin{indArg})
        case 'elec'
            elec = varargin{indArg+1};
            indArg = indArg+2;
        case 'gel'
            gel = varargin{indArg+1};
            indArg = indArg+2;
        case 'unitag'
            simTag = varargin{indArg+1};
            indArg = indArg+2;
        otherwise
            error('Supported options are: ''elec'', ''gel'', or ''uniTag''.');
    end
end
if exist('elec','var') && ~exist('gel','var')
    if strcmp(elec,'custom')
        recipe = {'custom1',-2,'custom2',2};
    else
        recipe = {elec{1},1,elec{2},-1};
        simTag = elec{1};
    end
end
[dirname,baseFilename] = fileparts(subj);
logname = fullfile(dirname,['hiperskipseg_' simTag '_log.txt']);
elecName = (recipe(1:2:end-1))';
elecPara = struct('capType',capType,'elecType',elecType,...
    'elecSize',elecSize,'elecOri',elecOri);
meshOpt = struct('radbound',5,'angbound',30,'distbound',0.4,'reratio',3,'maxvol',10);
% conductivities = struct('white',0.3835,'gray',0.1,'csf',1.8,'bone',0.0109,...
%     'skin',0.43,'air',2.5e-14,'gel',0.3,'electrode',5.9e7); % Indahlastari2016
conductivities = struct('white',0.126, 'gray',0.276,'csf',1.65,'bone',0.01, ...
    'skin',0.465,'air',2.5e-14,'gel',0.3,'electrode',5.9e7); % ROAST Def
% conductivities = struct('white',0.22, 'gray',0.47,'csf',1.71,'bone',0.02, ...
%     'skin',0.41,'air',2.5e-14,'gel',0.3,'electrode',5.9e7); % McCann et al 2019

if length(conductivities.gel(:))==1
    conductivities.gel = repmat(conductivities.gel,1,length(elecName));
end
if length(conductivities.electrode(:))==1
    conductivities.electrode = repmat(conductivities.electrode,1,length(elecName));
end
try
    [elecPara,~] = tdcslab_elecPreproc(subj,elecName,elecPara,day);
    disp 'PRE-PLACE COMPLETE !'; logfile(logname,'PRE-PLACE COMPLETE !');
catch
    logfile(logname,'PRE-PLACE FAILED ...');
    disp 'PRE-PLACE FAILED ...'; return
end
injectCurrent = (cell2mat(recipe(2:2:end)))';
configTxt = [];
for i=1:length(elecName)
    configTxt = [configTxt elecName{i} ' (' num2str(injectCurrent(i)) ' mA), '];
end
configTxt = configTxt(1:end-2);
options = struct('configTxt',configTxt,'elecPara',elecPara,'T2',[],'meshOpt',meshOpt,'conductivities',conductivities,'uniqueTag',simTag,'resamp',doResamp,'zeroPad',paddingAmt);
% load(fullfile(dirname,[baseFilename '_T1orT2_seg8.mat']),'Affine','image','lkp','ll','mg','mn','MT','tpm','Twarp','vr','wp')
% image.fname = subj;
% image.private.dat.fname = subj;
% % image.dim = hdr.dime.dim(2:4);
% save(fullfile(dirname,[baseFilename '_T1orT2_seg8.mat']),'Affine','image','lkp','ll','mg','mn','MT','tpm','Twarp','vr','wp')
% clear('Affine','image','lkp','ll','mg','mn','MT','tpm','Twarp','vr','wp')
if ~exist(fullfile(dirname,[baseFilename '_' simTag '_mask_elec.nii']),'file')
    %     hdrInfo = tdcslab_ElecPlace(day,subj,subj,[],elecName,options,simTag);
    try
        hdrInfo = electrodePlacement(subj,subj,[],elecName,options,simTag);
        disp 'PLACE COMPLETE !'; logfile(logname,'PLACEMENT COMPLETE !');
    catch
        logfile(logname,'PLACE FAILED ...');
        disp 'PLACE FAILED ...'; return
    end
else
    load(fullfile(dirname,[baseFilename '_header.mat']),'hdrInfo');
end
if ~exist(fullfile(dirname,[baseFilename '_' simTag '.mat']),'file')
    % hdr = load_untouch_header_only(subj);
    % v2w = diag(ones(1,4));
    % v2w(:,4) = 1;
    % hdrInfo = struct('pixdim',hdr.dime.pixdim(2:4),'dim',hdr.dime.dim(2:4),'v2w',v2w);
    % [node,elem,hdrInfo] = meshMYelec(s,subj,elec,gel,meshOpt,uniqueTag);
    try
        [node,elem,~] = meshByIso2mesh(s,subj,subj,[],meshOpt,hdrInfo,simTag);
        disp 'MESH COMPLETE !'; logfile(logname,'MESH COMPLETE !');
    catch
        logfile(logname,'MESH FAILED ...');
        disp 'MESH FAILED ...'; return
    end
else
    load(fullfile(dirname,[baseFilename '_' simTag '.mat']),'node','elem')
end
% if ~exist(fullfile(dirname,[baseFilename '_' simTag '_roastResult.mat']),'file')
    try prepareForGetDP(subj,node,elem,elecName,simTag); logfile(logname,'PREPARE COMPLETE !'); catch; logfile(logname,'PREPARE FAILED ...'); end
    try solveByGetDP(subj,injectCurrent,conductivities,1:length(elecName),simTag,[]); logfile(logname,'SOLVE COMPLETE !'); catch; logfile(logname,'SOLVE FAILED ...'); end
    try postGetDP(subj,subj,node,hdrInfo,simTag); logfile(logname,'BOOM, ROASTED !!!'); catch; logfile(logname,'POST FAILED ...'); end
% else
%     logfile(logname,'BOOM, ROASTED !!!');
% end
