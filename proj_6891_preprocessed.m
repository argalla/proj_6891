%% 1. Clear workspace
clear all; close all

%% 2. Import preprocessed data (provided by Wei - Yi, lead author of paper)
% problem - takes 10 minutes just to input data!
tic
% ov_gse = tdfread('ov.gse9891.19190x285.jetset.ncbi.txt');
% ov_tcga = tdfread('ov.tcga.11963x582.ncbi.txt');
 co_gse = tdfread('coad.gse14333.19190x290.jetset.ncbi.txt');
% co_tcga = tdfread('coad.tcga.17475x154.ncbi.txt');
% br_gse = tdfread('brca.gse2034.12160x286.jetset.ncbi.txt');
% br_tcga = tdfread('brca.tcga.17475x536.ncbi.txt');
toc

%% 3. Run attractor finding algorithm

% index all samples
samples = fieldnames(co_gse);
num_samples = numel(samples);
genes = co_gse.ExpID;
num_genes = length(genes);

% form matrix of all expression values, gene x sample
ex = zeros(num_genes,num_samples-1);
for i = 2:num_samples
   ex(:,i-1) = co_gse.(samples{i}); 
end

%% find attractor for each metagene
%notes - after running for ~14 hours got through 6697 genes as seed
%7.7 seconds per gene after optimizing
%94 pct of time is just doing mi calculations - hard to make that faster

%topgenes = cell(num_genes,1);
act_samples = num_samples - 1;


%seed = 8478;    %LCP2, ov_gse
%seed = 16830;   %THBS2, ov_gse
%seed = 1579;    %BUB1, ov_gse
%seed = 1724;    %CD53, ov_tcga
%seed = 2223;    %COL5A1, ov_tcga
%seed = 5437;    %KIF4A, ov_tcga
%seed = 17775;   %TYROBP, co_gse
seed = 5554;    %FBN1, co_gse
%seed = 1580;    %BUB1B, co_gse
%seed = 8079;    %LAPTM5, co_tcga
%seed = 690;     %ANTXR1, co_tcga
%seed = 1548;    %BUB1, co_tcga
%seed = 8585;    %PTPRC, br_gse
%seed = 2312;    %COL5A2, br_gse
%seed = 1975;    %CENPA, br_gse
%seed = 2636;    %CD5, br_tcga
%seed = 6021;    %GLT8D2, br_tcga
%seed = 15706;   %TPX2, br_tcga
% run algorithm
iter = 10;
alpha = 5;


% calculate mutual information between seed and all genes
MI = zeros(num_genes,1);
for n = 1:num_genes
    MI(n,1) = mi(ex(seed,:)',ex(n,:)')/max(...
        mi(ex(seed,:)',ex(seed,:)'),...
        mi(ex(n,:)',ex(n,:)'));
end
w = zeros(num_genes,2);
w(:,1) = MI(:,1).^alpha;
% create metagene
for i = 1:iter
    mg_pre = zeros(size(ex));
    for n = 1:act_samples
        mg_pre(:,n) = w(:,1).*ex(:,n);
    end

    mg = sum(mg_pre);
    %calculate mutual information between metagene and all genes
    for n = 1:num_genes
        MI(n,1) = mi(mg',ex(n,:)')/max(...
            mi(mg',mg'),...
            mi(ex(n,:)',ex(n,:)'));
    end
    w(:,2) = MI.^alpha;
    if max(w(:,2)) ~= w(seed,2) || max(w(:,2)) == 0
        break
    else w(:,1) = w(:,2);
    end

end
[mg_sort,idx] = sort(w(:,1),'descend');
attractor_idx = idx(1:50);
attractor = genes(idx(1:50),:); attractor_cell = cellstr(attractor);

%write to excel sheet to compare to supp. data
xlswrite('Attractors.xlsx',attractor_cell,'seed=FBN1','A1');





