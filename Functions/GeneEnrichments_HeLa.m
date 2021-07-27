function [enrichment_table] = GeneEnrichments(GeneList,CellLine)
% Select only CRC cell lines
load Hela_EG.mat

% for all possible genes in Recon
load('recon_genes.mat')

GeneList_all = unique(recon_genes.SYMBOL); % 1729 genes in Recon2
M = numel(GeneList_all);
enrichment_table = [];

for i = 1:size(GeneList,2)
    N = numel(GeneList{i});
    Non_EG_names = setdiff(GeneList_all,GeneList{i}); %Predicted non essential genes
    Non_EG = numel(Non_EG_names);

    % We want to check if our predicted essential genes (GeneList) the cell
    % line under study (CellLine)
    idx = (find(table2array(data(:,2))==1));
    cancer_genes = table2array(data(idx,1));
    all_genes = lower(cancer_genes);
    Known_Non_EG_names = setdiff(GeneList_all,cancer_genes); %Genes known to be non essential
    Known_Non_EG = numel(Known_Non_EG_names);

    K=[];
    x=[];
    n = [];
    K = [K;sum(ismember(lower(GeneList_all), all_genes))];
    x = [x;sum(ismember(lower(GeneList{i}), all_genes))];
    n = [n;sum(ismember(lower(Non_EG_names), lower(Known_Non_EG_names)))]; %number of predicted non_EG that are known as non_EG
    enrichment = table({'HeLa'},num2cell(1-hygecdf(x-1,M,K,N)),num2cell(M),num2cell(K),num2cell(N),num2cell(x),num2cell(Known_Non_EG),num2cell(Non_EG),num2cell(n),...
            'VariableNames',{'Cell_Line','enrichment','Metabolic_genes','Metabolic_cancer_genes','Predicted_EG','Known_EG','non_EG_in_model','Predicted','Known_non_EG'});
    enrichment_table = vertcat(enrichment_table,enrichment);
end
end
