Multitask Clustering of scRNAseq Data]{Multitask clustering of single-cell RNA-seq data identifies cell subpopulations and markers in recessive dystrophic epidermolysis bullosa


nature13173-s4.txt

This file contains single cell RNA-seq expression data (log2(FPKM) values) for all 80 lung epithelial cells at E18.5 together with the putative cell type of each cell in a .txt file. The last two are bulk values.

Matlab_ProcessData.m
Octave_ProcessData.m

Process data with Matlab or Ocatve.


Matlab_run_test.m
Octave_run_test.m

Run scVDMC algorithm with Matlab or Ocatve.

% cleandata: cell by genes matrix, after filtering 
% data: cell by genes matrix, before filtering 
% genelist: filtered gene list
% Cellname: each single cell with name
% Celltype: each single cell with cell type
% Celltype_list: cell type list
% CellRep: each single cell with replicate information
# scVDMC
