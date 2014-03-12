function [names,old_names,new_names] = tcga_filerenamer(dr,varargin)
% dr - directory containing files to be renamed
% varargin - if renaming tcga ovarian samples set varargin = 'ovarian'
%       note: this may generalize for all affymetrix samples
% Rename files for TCGA datasets
% REQUIRES FileRename function from Jan Simon: download here:
% http://www.mathworks.com/matlabcentral/fileexchange/29569-filerename 
% Also need to download compiled c-mex file, or compile yourself; compiled
% mex file at: http://www.n-simon.de/mex/Matlab76_64/FileRename.mexw64
% Modeled off Oleg Komarov's solution from MATLAB community site: 
% http://www.mathworks.com/matlabcentral/answers/1760-how-to-rename-a-bunch-of-files-in-a-folder

names = dir(dr); names = names(3:end); % first 2 files not samples - find better way to do this, but seems consistent
old_names = {names(~[names.isdir]).name}; % original name has junk  that confuses matlab

if strcmp(varargin,'ovarian') == 1
    %get rid of junk at end of filename - ovarian samples (affymetrix)
    for i = 1:numel(names)
        x = find(names(i).name == '.'); names(i).name = [names(i).name(1:x) 'txt'];        
    end
else
    %get rid of junk at end of filename - works for brca and colon (agilent)
    for i = 1:numel(names)
        x = find(names(i).name == '.'); names(i).name = names(i).name(1:x+3);
    end
end

new_names = {names(~[names.isdir]).name}; % ends at .txt, only includes filename, not path
if strcmp(new_names(1),old_names(1)) == 0
    for i = 1:numel(names)
        oldname = [dr old_names{i}]; 
        newname = [dr new_names{i}];
        FileRename(oldname,newname)
    end
else
    display('Files already renamed, no action taken');
end