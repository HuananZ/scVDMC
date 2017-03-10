clc;
clear;
close all;

%%
% This file contains single cell RNA-seq expression data (log2(FPKM) values)
% for all 80 lung epithelial cells at E18.5 together with the putative cell type 
% of each cell in a .txt file.
% last two are bulk values


% T = readtable('nature13173-s4.txt', 'delimiter', '\t');
% genelist_full = T.Properties.VariableNames(5:end)';
% Cellname = table2cell(T(1:80, 1));
% Celltype = table2cell(T(1:80, 4));
% data = table2array(T(1:80, 5:end));

fid = fopen('nature13173-s4.txt');
headerline = textscan(fid,'%s', 23275);
fclose(fid);
genelist_full = headerline{1}(5:23275);

fid = fopen('nature13173-s4.txt');
textformat = ['%s%s%s%s', repmat('%f',1,23271)];
C = textscan(fid, textformat, 'delimiter', '\t', 'HeaderLines',1);
fclose(fid);
Cellname = C{1}(1:80);
Celltype = C{4}(1:80);
data = cell2mat(C(5:end));
data = data(1:80,:);


% replicate information
CellRep = zeros(80,1);
for i = 1:80    
    tmp = regexp(Cellname{i}, '_', 'split');
    CellRep(i) = str2double(tmp{2});
end;

Celltype_list = unique(Celltype);


%% keep high mean features, keep rare but high value features
genelist = genelist_full((sum(data > 1) > 2));
cleandata = data(:, (sum(data > 1) > 2));

save Octave_lungscdata.mat data cleandata genelist Cellname Celltype Celltype_list CellRep;