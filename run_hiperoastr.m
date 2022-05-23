rootDir = pwd;
outDir = '/blue/camctrp/working/aprinda/PRECISE/SB';
subfdr = dir(fullfile(outDir,'sub*'));
subnames = {subfdr.name}';
parfor s = 1:length(subnames)
    T1 = fullfile(outDir,subnames{s},[subnames{s} '_ses-01_T1w.nii']);
    T2 = fullfile(outDir,subnames{s},[subnames{s} '_ses-01_FLAIR.nii']);
%     try
    if exist(T2,'file')
        hiperoastr(s,T1,'F3-F4',T2)
    else
        hiperoastr(s,T1,'F3-F4')
    end
        disp([subnames{s} ' Complete !'])
%     catch
%         disp([subnames{s} ' FAILED !!!!!!!!!!!!!! BROTHERRRRRR!'])
%     end
end
disp('HELLLLL YELAHJJJJJ BORTHERRRRRRRR!!!')