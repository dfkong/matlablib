% RESTORE_IDL_DEMO
%
% finds the IDL save files included in the restore_idl distribution and
% restores each of them in turn

me = evalc('which restore_idl_demo');
[path, file, type]=fileparts(me);
savefiles=dir(fullfile(path,'idl_save','*.sav'));
for i=1:numel(savefiles)
    fname=savefiles(i).name;
    disp(['press any key to restore ' fname]);
    pause;
    outargs=restore_idl(fullfile(path,'idl_save',fname),'verbose');
    outargs;
    disp('Contents of outargs:');
    fields=fieldnames(outargs)
    for n=1:numel(fields)
        disp(fields{n});
        thevar=outargs.(fields{n});
        whos thevar
        if ischar(thevar) || numel(thevar) <= 10,
            thevar
        end
    end
end