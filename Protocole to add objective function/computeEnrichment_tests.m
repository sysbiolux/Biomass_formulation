function [enrichment_table] = True_False_positives_template_code(predictedEssential,CellLine, essential_data,colnames_essential, rownames_essential, dico, all_unique_genes_in_model)
% Select only CRC cell lines
if ischar(CellLine)
  CellLine=cellstr(CellLine); 
  disp('Warning CellLine should be a cell array.')
  disp('CellLine was converted to a cell array.')
end
    
if numel(CellLine)~=numel(predictedEssential)
    disp('CellLine should be a cell array and match the number of fields in Essential')
end 
col = find(sum(ismember(dico,rownames_essential)) == max(sum(ismember(dico,rownames_essential)))); % find matching column, usually the one with the highest number of matches
col2 = find(sum(ismember(dico,all_unique_genes_in_model)) == max(sum(ismember(dico,all_unique_genes_in_model)))); % find matching column, usually the one with the highest number of matches
if col==0
    disp('the dico does not match the data')
elseif numel(col)>1
    col=col(1);
end
if col2==0
    disp('the dico does not match the model')
elseif numel(col2)>1
    col2=col2(1);
end

[~,idico,irownames] = intersect(dico(:,col),rownames_essential); % get indices
essential_data = essential_data(irownames,:); %get essential_data for matching rownames
rownames_essential=lower(dico(idico, col2));

% for all possible genes in Recon
M = numel(all_unique_genes_in_model);
enrichment_table = [];
for z = 1:size(CellLine,1)
    if any((ismember(colnames_essential,CellLine(z)))==1)
        cell_line = find(ismember(colnames_essential,CellLine(z)));

        % for GeneList
        % GeneList = GeneList_all(randi([1 1175],1,50)); %for testing

        N = numel(unique(predictedEssential(z).geneList)); %predicted EG
        Non_EG_names = setdiff(all_unique_genes_in_model,predictedEssential(z).geneList); %Predicted non essential genes
        Non_EG = numel(Non_EG_names);
        
        % We want to check if our predicted essential genes (GeneList) the cell
        % line under study (CellLine)
        [r,~]=find(essential_data(:,cell_line)==1);
        cancer_genes=unique(rownames_essential(r));
        Known_Non_EG_names = setdiff(all_unique_genes_in_model,cancer_genes); %Genes known to be non essential
        Known_Non_EG = numel(Known_Non_EG_names);
        
        K=[];
        x=[]; 
        n = [];
        K = [K;sum(ismember(lower(all_unique_genes_in_model), cancer_genes))]; %number of metabolic genes present in the CRISPR-Cas 9 data
        x = [x;sum(ismember(unique(predictedEssential(z).geneList), cancer_genes))]; %number of predicted EG that are known EG
        n = [n;sum(ismember(lower(Non_EG_names), lower(Known_Non_EG_names)))]; %number of predicted non_EG that are known as non_EG

        pval = table(colnames_essential(cell_line),num2cell(1-hygecdf(x-1,M,K,N)),num2cell(M),num2cell(K),num2cell(N),num2cell(x),num2cell(Known_Non_EG),num2cell(Non_EG),num2cell(n),...
            'VariableNames',{'Cell_Line','p-value','Metabolic_genes','Metabolic_cancer_genes','Predicted_EG','Known_EG','non_EG_in_model','Predicted_non_EG','Known_non_EG'});
        enrichment_table = vertcat(enrichment_table,pval);
     
    end
    
    
end
Sensitivity=[];
Specificity=[];
Precision=[];
for i=1:size(enrichment_table,1)
     sensitivity=cell2mat(enrichment_table.Known_EG(i))/cell2mat(enrichment_table.Metabolic_cancer_genes(i));
    specificity=cell2mat(enrichment_table.Known_non_EG(i))/cell2mat(enrichment_table.non_EG_in_model(i));
    precision=cell2mat(enrichment_table.Known_EG(i))/cell2mat(enrichment_table.Predicted_EG(i));
    
    Sensitivity = [Sensitivity,sensitivity];
    Specificity = [Specificity,specificity];
    Precision = [Precision,precision];
end
names=enrichment_table.Properties.VariableNames;
enrichment_table=table2array(enrichment_table);
enrichment_table(:,end+1)=num2cell(Sensitivity');
enrichment_table(:,end+1)=num2cell(Specificity');
enrichment_table(:,end+1)=num2cell(Precision');
enrichment_table=cell2table(enrichment_table);
names(end+1:end+3)={'Sensitivity', 'Specificity', 'Precision'};
enrichment_table.Properties.VariableNames= names;



end