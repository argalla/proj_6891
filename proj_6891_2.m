%% 1. Clear Workspace
clear all; close all

%% 2. Rename all files for TCGA data
% REQUIRES FileRename function from Jan Simon: download here:
% http://www.mathworks.com/matlabcentral/fileexchange/29569-filerename 
% Also need to download compiled c-mex file, or compile yourself; compiled
% mex file at: http://www.n-simon.de/mex/Matlab76_64/FileRename.mexw64

% Modeled off Oleg Komarov's solution from MATLAB community site: 
% http://www.mathworks.com/matlabcentral/answers/1760-how-to-rename-a-bunch-of-files-in-a-folder

% file directory with all TCGA breast samples (604)
% safe copy of directory on proj_6891 folder
dr_brca = 'C:\Users\Arthur\Documents\MATLAB\breast_all_tcga\';
% file directory with TCGA colon samples (181)
dr_colon = 'C:\Users\Arthur\Documents\MATLAB\colon_all_tcga\';
%file directory with TCGA ovarian samples (594)
dr_ovarian = 'C:\Users\Arthur\Documents\MATLAB\ovarian_tcga\';
numgenes = 17814; %number of genes assayed for colon and brca
numgenes_ovarian = 12042; %number of genes assayed for ovarian

[names_brca,old_names_brca,new_names_brca] = tcga_filerenamer(dr_brca);
[names_colon,old_names_colon,new_names_colon] = tcga_filerenamer(dr_colon);
[names_ovarian,old_names_ovarian,new_names_ovarian] = tcga_filerenamer(dr_ovarian,'ovarian');

%% 3. Import tcga files as tdf and form array (genes x samples), perform data preprocessing
% note: This can take a long time, tcga raw data need to be converted from
% strings.
% note: can use tcga_arraymaker or tcga_arraymaker_2. 2 is better for
% memory purposes as it converts the cell array to a matrix. Think I should
% keep using 2 going forward. PROBLEM - some empty entries.
% for colon empty value at (23,161)
% temporary solution: just fill empty cells with average of given gene expression
% across all samples

[ex_colon,barcodes_colon,genelist_colon,emptycells_colon] ...
    = tcga_arraymaker_2(dr_colon,new_names_colon,numgenes);
[ex_brca,barcodes_brca,genelist_brca,emptycells_brca] ...
    = tcga_arraymaker_2(dr_brca,new_names_brca,numgenes);
[ex_ovarian,barcodes_ovarian,genelist_ovarian, emptycells_ovarian] ...
    = tcga_arraymaker_2(dr_ovarian,new_names_ovarian,numgenes_ovarian);

%% 4. Import expression data from GEO
% expression data, also includes metadata such as clinical outcomes
% note: this requires an internet connection
GEO_breast = getgeodata('GSE2034'); GEO_colon = getgeodata('GSE14333');
GEO_ovarian = getgeodata('GSE9891'); num_geo = 3;
% column names indicate sample (GSM ID)
samples_breast = colnames(GEO_breast.Data); samples_colon = colnames(GEO_colon.Data);
samples_ovarian = colnames(GEO_ovarian.Data);
% row names: ID ref
rows_breast = rownames(GEO_breast.Data); rows_colon = rownames(GEO_colon.Data);
rows_ovarian = rownames(GEO_ovarian.Data);

%% 5. Get platform data - information from Affymetrix about probes
% need this to map rows to genes
% note: this requires an internet connection
platform = getgeodata('GPL96');
probeIDs = platform.Data(:,strcmp(platform.ColumnNames,'ID'));
platform_genes = platform.Data(:,strcmp(platform.ColumnNames,'Gene Symbol'));

%% 6. Extract raw .cel files (affymetrix) from .gz
dr_raw_breast = 'E:\Data\proj_6891\GSE2034\'; 
dr_breast_unzipped = 'E:\Data\proj_6891\GSE2034_unzipped\';
dr_raw_colon = 'E:\Data\proj_6891\GSE14333\';
dr_colon_unzipped = 'E:\Data\proj_6891\GSE14333_unzipped\';
dr_raw_ovarian = 'E:\Data\proj_6891\GSE9891\';
dr_ovarian_unzipped = 'E:\Data\proj_6891\GSE9891_unzipped\';

breast_cel = dir(dr_breast_unzipped);
colon_cel = dir(dr_colon_unzipped);
ovarian_cel = dir(dr_ovarian_unzipped);

% Only do this once; otherwise an error is asserted as the directories are
% already filled with cel files
%cel_extractor(dr_raw_breast,dr_breast_unzipped);
%cel_extractor(dr_raw_colon,dr_colon_unzipped);
%cel_extractor(dr_raw_ovarian,dr_ovarian_unzipped);

% Includes 3 GEO data sets and ovarian TCGA
% Note: still need to get .cel files for ovarian TCGA
CDF_dr = 'E:\Data\proj_6891\CDF Files\'; %directory with CDF files - necessary for RMA algorithm

% CDF info: breast - U133A
%           colon  - U133Plus 2.0
%          ovarian - U133Plus 2.0
U133A_cdf = 'HG-U133A.cdf';
U133Plus_cdf = 'HG-U133_Plus_2.cdf';

% Make cell arrays of all cel files, 1 array each for breast, cancer, colon
breast_cel_names = {breast_cel(3:end).name};
colon_cel_names = {colon_cel(3:end).name};
ovarian_cel_names = {ovarian_cel(3:end).name};

%% 7. RMA implementation
breast_rma = affyrma(breast_cel_names,U133A_cdf,'CelPath',dr_breast_unzipped, ...
    'CDFPath',CDF_dr);
save('E:\Data\proj_6891\breast_rma.mat','breast_rma');
clearvars 'breast_rma'

colon_rma = affyrma(colon_cel_names,U133Plus_cdf,'CelPath',dr_colon_unzipped,...
    'CDFPath',CDF_dr);
save('E:\Data\proj_6891\colon_rma.mat','colon_rma');
clearvars 'colon_rma'

ovarian_rma = affyrma(ovarian_cel_names,U133Plus_cdf,'CelPath',dr_ovarian_unzipped,...
    'CDFPath',CDF_dr);
save('E:\Data\proj_6891\ovarian_rma.mat','ovarian_rma');
clearvars 'ovarian_rma'









