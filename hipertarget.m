%% Skip Seg
function hipertarget(s,subj,anode,cathode,simTag,logname)

% Settings
%=======================================
recipe = {cathode, -2, anode, 2};
elecType = {'pad','pad'};
elecSize = {[70 50 3],[70 50 3]};
elecOri = 'lr';
capType = '1020';
paddingAmt = 0;
doResamp = 0;
%=======================================
[dirname,baseFilename] = fileparts(subj);
elecName = (recipe(1:2:end-1))';
elecPara = struct('capType',capType,'elecType',elecType,...
    'elecSize',elecSize,'elecOri',elecOri);
meshOpt = struct('radbound',5,'angbound',30,'distbound',0.4,'reratio',3,'maxvol',10);
conductivities = struct('white',0.3835,'gray',0.1,'csf',1.8,'bone',0.0109,...
    'skin',0.43,'air',2.5e-14,'gel',0.3,'electrode',5.9e7); % Indahlastari2016
if length(conductivities.gel(:))==1
    conductivities.gel = repmat(conductivities.gel,1,length(elecName));
end
if length(conductivities.electrode(:))==1
    conductivities.electrode = repmat(conductivities.electrode,1,length(elecName));
end
try
    [elecPara,~] = elecPreproc(subj,elecName,elecPara);
catch
    logfile(logname,[subj ' ' anode '+' cathode ' FAILED at elec PRE ...']);
end
injectCurrent = (cell2mat(recipe(2:2:end)))';
configTxt = [];
for i=1:length(elecName)
    configTxt = [configTxt elecName{i} ' (' num2str(injectCurrent(i)) ' mA), '];
end
configTxt = configTxt(1:end-2);
options = struct('configTxt',configTxt,'elecPara',elecPara,'T2',[],'meshOpt',meshOpt,'conductivities',conductivities,'uniqueTag',simTag,'resamp',doResamp,'zeroPad',paddingAmt);
load(fullfile(dirname,[baseFilename '_T1orT2_seg8.mat']),'Affine','image','lkp','ll','mg','mn','MT','tpm','Twarp','vr','wp')
image.fname = subj;
image.private.dat.fname = subj;
% image.dim = hdr.dime.dim(2:4);
save(fullfile(dirname,[baseFilename '_T1orT2_seg8.mat']),'Affine','image','lkp','ll','mg','mn','MT','tpm','Twarp','vr','wp')
clear('Affine','image','lkp','ll','mg','mn','MT','tpm','Twarp','vr','wp')
if ~exist(fullfile(dirname,[baseFilename '_header.mat']),'file')
    try
        hdrInfo = electrodePlacement(subj,subj,[],elecName,options,simTag);
    catch
        logfile(logname,[subj ' ' anode '+' cathode ' FAILED at elec PLACE ...']);
        return
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
        [node,elem,~] = par_meshByIso2mesh(s,subj,subj,[],meshOpt,hdrInfo,simTag);
    catch
        logfile(logname,[subj ' ' anode '+' cathode ' FAILED at MESH ...']);  
        return
    end
else
    load(fullfile(dirname,[baseFilename '_' simTag '.mat']),'node','elem')
end
try
    prepareForGetDP(subj,node,elem,elecName,simTag);
catch
    logfile(logname,[subj ' ' anode '+' cathode ' FAILED at PREPARE ...']);
    return
end
try
    solveByGetDP(subj,injectCurrent,conductivities,1:length(elecName),simTag,[]);
catch
    logfile(logname,[subj ' ' anode '+' cathode ' FAILED at SOLVE ...']);
    return
end
try
    postGetDP(subj,subj,node,hdrInfo,simTag);
catch
    logfile(logname,[subj ' ' anode '+' cathode ' FAILED at elec POST ...']);
    return
end