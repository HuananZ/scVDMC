clc;
clear;
close all;

load Matlab_lungscdata;

% cleandata: cell by genes matrix, after filtering 
% data: cell by genes matrix, before filtering 
% genelist: filtered gene list
% Cellname: each single cell with name
% Celltype: each single cell with cell type
% Celltype_list: cell type list
% CellRep: each single cell with replicate information

%% separate data by replicate information
firsttime = 0;
d = 3;
max_iter = 100;

n = zeros(d,1);
sp_data = cell(d,1);
sp_label = cell(d,1);
for dd = 1:d
    sp_data{dd} = cleandata(CellRep == dd, :)';
    [m, n(dd)] = size(sp_data{dd});
    sp_label{dd} = Celltype(CellRep == dd);
end;


%% remove ciliated since only shown in data3
% then all data have 4 labels
ix = find(strcmp(sp_label{3}, 'ciliated'));
sp_data{3}(:,ix) = [];
sp_label{3}(ix) = [];
n(3) = n(3) - length(ix);

k = 4;
% sort sample by cell type
for dd = 1:d
    [sp_label{dd}, IX] = sort(sp_label{dd});
    sp_data{dd} = sp_data{dd}(:,IX);
end;

%% get unuseable features
sp_data_orig = sp_data;

nanix = [];
for dd = 1:d
    mD = mean(sp_data{dd}, 2);
    vD = std(sp_data{dd}, [], 2);
    sp_data{dd} = (sp_data{dd} - repmat(mD, 1, n(dd))) ./ repmat(vD, 1, n(dd));
    nanix = union(nanix, find(vD == 0));
end;

for dd = 1:d
    sp_data{dd}(nanix,:) = [];
    sp_data_orig{dd}(nanix,:) = [];
end;

m_clean = m - length(nanix);
genelist_clean = genelist;
genelist_clean(nanix, :) = [];

%% get pool, unnormalized
% pool data
pool_data = cell2mat(sp_data_orig');
pool_label = [];
for dd = 1:d
    pool_label = [pool_label; sp_label{dd}];
end;

%% original data, check the center and means
V_True_sp = cell(d, 1);
thelabels = unique(pool_label);
for dd = 1:d
    V_True_sp{dd} = zeros(n(dd), k);
    for i = 1:k
        V_True_sp{dd}(strcmp(thelabels{i}, sp_label{dd}), i) = 1;
    end;
end;

%% pool_mean
IX = kmeans(pool_data', k, 'Distance', 'correlation', 'Replicates', 20);

V_PKMS = zeros(sum(n), k);
for i = 1:k
    V_PKMS(IX == i, i) = 1;
end;

%% using kmeans result as our start point
wpool = [0.1 1];
lambdapool = 10:10:200;

err_VDMC = zeros(d, length(wpool), length(lambdapool));

U_VDMC = cell(d, length(wpool), length(lambdapool));
V_VDMC = cell(d, length(wpool), length(lambdapool));

VDMC_sorted_gene = cell(length(wpool), length(lambdapool));
VDMC_sorted_IX = cell(length(wpool), length(lambdapool));

V_ini = cell(d,1);
U_ini = cell(d,1);

iix = [0; cumsum(n)];

% use pool kmeans as start point
V_ini{1} = V_PKMS(iix(1)+1:iix(2),:);
V_ini{2} = V_PKMS(iix(2)+1:iix(3),:);
V_ini{3} = V_PKMS(iix(3)+1:iix(4),:);

for dd = 1:d
    for kk = 1:k
        U_ini{dd}(:, kk)  = mean(sp_data{dd}(:, V_ini{dd}(:,kk) == 1), 2);
    end;
end;


% use house method
a = perms(1:k);
for wi = 1:length(wpool)
    w = wpool(wi);
    for li = 1:length(lambdapool)
        lambda = lambdapool(li);
        [U_VDMC_once, V_VDMC_once, B_TF_once, sortB_TF_once, Obj_TF_once]...
            = VDMC(sp_data, d, k, w, lambda, U_ini, V_ini, max_iter);
        
        for dd = 1:d            
            err = zeros(size(a,1), 1);
            for i = 1:size(a,1)
                err(i) = length(find(sum(V_True_sp{dd} ~= V_VDMC_once{dd}(:, a(i,:)), 2)));
            end;
            [err_VDMC(dd,wi,li), ix] = min(err);
            V_VDMC{dd,wi,li} = V_VDMC_once{dd}(:, a(ix(1),:));
            U_VDMC{dd,wi,li} = U_VDMC_once{dd}(:, a(ix(1),:));            
        end;
        
        VDMC_sorted_gene{wi, li} = genelist_clean(sortB_TF_once);
        VDMC_sorted_IX{wi, li} = sortB_TF_once;
    end;
end;


save('Lung_VDMC.mat',...
    'wpool', 'lambdapool',...
    'VDMC_sorted_gene', 'VDMC_sorted_IX', ...
    'err_VDMC', 'U_VDMC', 'V_VDMC');


%% disp
FigHandle = figure('Position', [100, 100, 1000, 600]);
plot(squeeze(sum(err_VDMC, 1))', 'linewidth', 3, 'LineStyle', '-.');
set(gca, 'xtick', 1:length(lambdapool), 'XTickLabel', lambdapool);
title('Lung Single Cell');
xlabel('\lambda');
ylabel('err');