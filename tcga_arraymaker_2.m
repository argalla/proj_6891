function [ex_mat,barcodes,genelist,emptycells] = tcga_arraymaker_2(dr,new_names,numgenes)
% numgenes 17814 for brca and colon, 12042 for ovarian
% new_names and dr from tcga_filerenamer 
% Import files as tdf and form array, genes x sample
% ver2 outputs gene expression data as matrix instead of cell array; this will
% use less memory (with cell array can take up 3 GB or so)
% imputes missing values with row mean

num_samples = numel(new_names);
barcodes = cell(1,num_samples);
ex = cell(numgenes,num_samples);

for i = 1:num_samples
    s = tdfread([dr new_names{i}]);
    fnames = fieldnames(s);
    barcodes(i) = fnames(2);
    sample = s.(fnames{2});
    for ii = 2:numgenes+1
        ex(ii-1,i) = textscan(sample(ii,:),'%f'); %improve by not going through full column?
    end
end

% Remove empty cells from array - just take mean of gene expression across
% samples
emptycells = cellfun('isempty',ex); 
[x,y] = find(emptycells);

for i = 1:length(x)
    ex(x(i),y(i)) = {mean([ex{x(i),1:end}])};
end
       
ex_mat = cell2mat(ex);
genelist = s.(fnames{1});