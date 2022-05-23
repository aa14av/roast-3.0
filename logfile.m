function fid = logfile(logname,msg)
if ~exist(logname,'file')
    save(logname,'msg','-ascii')
end
fid = fopen(logname, 'a');
fprintf(fid, '%s: %s\n', datestr(now, 0), msg);
fclose(fid);