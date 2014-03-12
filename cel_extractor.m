function cel_extractor(dr,output_dr)
% dr is directory containing .cel.gz files
% function extracts to output_dr 
assert(length(dir(output_dr)) <= 2 || isempty(dir(output_dr)) == 1,...
    'Output directory already contains files');
        
files = dir(dr); files = files(3:end); % first 2 files not samples - find better way to do this, but seems consistent
full_files = cell(1,numel(files));
for i = 1:numel(files)
    full_files{i} = [dr files(i).name];
end
% unzip files to new directory
gunzip(full_files,output_dr);

