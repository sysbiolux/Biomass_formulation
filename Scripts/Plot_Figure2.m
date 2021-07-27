altcolor= [255 255 255;255 204 204; 255 153 153; 255 102 102; 255 51 51;...
    255 0 0; 204 0 0; 152 0 0; 102 0 0;  51 0 0; 0 0 0]/255; %some red color
set(0,'defaultTextInterpreter','none')
set(groot, 'DefaultAxesTickLabelInterpreter', 'none')
set(groot, 'defaultLegendInterpreter','none');

%% Sensitivity and specificity analysis iHsa
% Sensitivity, specificity and precision analyses are run to assess the
% performance of the biomass prediction. The predicted essential genes are
% compared against CRISPR-Cas 9 screenings (DepMap Achilles 19Q2). 
% Sensitivity = TP/P
% Specificity = TN/N 
% Precision = TP/(TP+FP)

%% maint
%Data 
load dico_recon3D.mat
load intest_skin.mat
colnames = INTESTINE.Properties.VariableNames; %sample names

% I will run it only for the cell lines which have CRISPR data to avoid
% long running times.
load sample_info.mat
load binaryDepScores_broad_database.mat

rownames = binaryDepScores_broad_database.Gene;
crc = sampleinfo((find(contains(sampleinfo.lineage,'colorectal'))),:);

crc_binaryDepScores = [];
binaryDepScores = [];
for x = 2:size(binaryDepScores_broad_database.Properties.VariableNames,2)
    if ismember((erase(binaryDepScores_broad_database.Properties.VariableNames(x),'ACH')),convertStringsToChars(erase(crc.DepMap_ID,'ACH-')))
        loc = find(ismember(convertStringsToChars(erase(crc.DepMap_ID,'ACH-')),(erase(binaryDepScores_broad_database.Properties.VariableNames(x),'ACH'))));
        crc_binaryDepScores = [crc_binaryDepScores,crc.CCLE_Name(loc)];
        binaryDepScores = [binaryDepScores,binaryDepScores_broad_database(:,x)];
    end
end
binaryDepScores.Properties.VariableNames = crc_binaryDepScores;
binaryDepScores.Properties.RowNames = binaryDepScores_broad_database.Gene;


CRISPR = binaryDepScores.Properties.VariableNames;

idx_CRISPR = [];
for i=1:size(CRISPR,2)
    idx = find(ismember(colnames,CRISPR(i)))
    idx_CRISPR = [idx_CRISPR,idx];
end

CRC = [INTESTINE(:,[1:2]),INTESTINE(:,idx_CRISPR)];
fpkmIntestine = table2array(CRC(:,3:end)); %get the numbers
%% 
% Colnames

colnames = (CRC.Properties.VariableNames(3:end))'; %sample names

%In that case no different conditions
conditions = {[1:length(colnames)]};
condition_names = colnames;

%% 
% rownames

rownames = CRC.Name;

%remove NaN entries
rownames(isnan(sum(fpkmIntestine,2)),:) = [];
fpkmIntestine(isnan(sum(fpkmIntestine,2)),:) = [];

%remove 0 entries
rownames(sum(fpkmIntestine,2) == 0,:) = [];
fpkmIntestine(sum(fpkmIntestine,2) == 0,:) = [];
rownames = regexprep(rownames,'\.\d*',''); %remove transcripts

% get symbols for ENSG
[Lia,Locb] = ismember(rownames,recon3D_genes.ENSG);
rownames_SYMBOL = cell(size(rownames));
rownames_SYMBOL(Lia) = recon3D_genes.SYMBOL(Locb(find(Locb)));

%convert empty Symbols to char
for i=1:numel(rownames_SYMBOL)
    if ~ischar(rownames_SYMBOL{i})
        rownames_SYMBOL{i} = char(rownames_SYMBOL{i});
    end
end

% get ENTREZ for ENSG
[Lia,Locb] = ismember(rownames,recon3D_genes.ENSG);
rownames_ENTREZ = cell(size(rownames));
rownames_ENTREZ(Lia) = recon3D_genes.ENTREZ(Locb(find(Locb)));

%convert empty ENTREZ to char
for i=1:numel(rownames_ENTREZ)
    if ~ischar(rownames_ENTREZ{i})
        rownames_ENTREZ{i} = char(rownames_ENTREZ{i});
    end
end

[rownames, rownames_SYMBOL,rownames_ENTREZ]
clear ans i Lia Locb

load KO_generic_biomass_CCLE_medium_maint_input.mat 
essential_genes = zeros(size(grRatio_biomass_single,1),size(grRatio_biomass_single,2));
for i=1:size(grRatio_biomass_single,2)
threshold = 0.5;
B = grRatio_biomass_single(:,i) <= threshold;
essential_genes(B,i)=1;
end

EG_final=cell.empty(0, numel(colnames));
for i=1:size(essential_genes,2)
    x =(genelist(find(essential_genes(:,i)==1)));
    [~, ia,ib] = intersect(x,recon3D_genes.ENTREZ); % find the intersection and keep indices
    EG_final(i)= {recon3D_genes.SYMBOL(ib)};
end

confusion_table_maint = True_False_positives(EG_final,colnames,'Recon3D');
confusion_table_maint = sortrows(confusion_table_maint,'Cell_Line','ascend');

Sensitivity_Maintenance_input = [];
Specificity_Maintenance_input = [];
Precision_Maintenance_input = [];
for i = 1:size(confusion_table_maint,1)
    sensitivity = cell2mat(confusion_table_maint.Known_EG(i))/cell2mat(confusion_table_maint.Metabolic_cancer_genes(i));
    specificity = cell2mat(confusion_table_maint.Known_non_EG(i))/cell2mat(confusion_table_maint.non_EG_in_model(i));
    precision = cell2mat(confusion_table_maint.Known_EG(i))/cell2mat(confusion_table_maint.Predicted_EG(i)); %FP + TP is the same as the predicted genes
    Sensitivity_Maintenance_input = [Sensitivity_Maintenance_input,sensitivity];
    Specificity_Maintenance_input = [Specificity_Maintenance_input,specificity];
    Precision_Maintenance_input = [Precision_Maintenance_input,precision];
end

%% noTrTr
load KO_generic_biomass_CCLE_medium_noTrTr_input.mat 
essential_genes = zeros(size(grRatio_biomass_single,1),size(grRatio_biomass_single,2));
for i=1:size(grRatio_biomass_single,2)
threshold = 0.5;
B = grRatio_biomass_single(:,i) <= threshold;
essential_genes(B,i)=1;
end

EG_final=cell.empty(0, numel(colnames));
for i=1:size(essential_genes,2)
    x =(genelist(find(essential_genes(:,i)==1)));
    [~, ia,ib] = intersect(x,recon3D_genes.ENTREZ); % find the intersection and keep indices
    EG_final(i)= {recon3D_genes.SYMBOL(ib)};
end

confusion_table_noTrTr = True_False_positives(EG_final,colnames,'Recon3D');
confusion_table_noTrTr = sortrows(confusion_table_noTrTr,'Cell_Line','ascend');

Sensitivity_noTrTr_input = [];
Specificity_noTrTr_input = [];
Precision_noTrTr_input = [];
for i = 1:size(confusion_table_noTrTr,1)
    sensitivity = cell2mat(confusion_table_noTrTr.Known_EG(i))/cell2mat(confusion_table_noTrTr.Metabolic_cancer_genes(i));
    specificity = cell2mat(confusion_table_noTrTr.Known_non_EG(i))/cell2mat(confusion_table_noTrTr.non_EG_in_model(i));
    precision = cell2mat(confusion_table_noTrTr.Known_EG(i))/cell2mat(confusion_table_noTrTr.Predicted_EG(i)); %FP + TP is the same as the predicted genes
    Sensitivity_noTrTr_input = [Sensitivity_noTrTr_input,sensitivity];
    Specificity_noTrTr_input = [Specificity_noTrTr_input,specificity];
    Precision_noTrTr_input = [Precision_noTrTr_input,precision];
end

%% generic R3
load KO_GenericR3_biomass_CCLE.mat
essential_genes = zeros(size(grRatio_biomass_single,1),size(grRatio_biomass_single,2));
for i=1:size(grRatio_biomass_single,2)
threshold = 0.5;
B = grRatio_biomass_single(:,i) <= threshold;
essential_genes(B,i)=1;
end

EG_final=cell.empty(0, numel(colnames));
for i=1:size(essential_genes,2)
    x =(genelist(find(essential_genes(:,i)==1)));
    [~, ia,ib] = intersect(x,recon3D_genes.ENTREZ); % find the intersection and keep indices
    EG_final(i)= {recon3D_genes.SYMBOL(ib)};
end
confusion_table_genR3 = True_False_positives(EG_final,colnames,'Recon3D');
confusion_table_genR3 = sortrows(confusion_table_genR3,'Cell_Line','ascend');

Sensitivity_Generic_R3 = [];
Specificity_Generic_R3 = [];
Precision_Generic_R3 = [];
for i = 1:size(confusion_table_genR3,1)
    sensitivity = cell2mat(confusion_table_genR3.Known_EG(i))/cell2mat(confusion_table_genR3.Metabolic_cancer_genes(i));
    specificity = cell2mat(confusion_table_genR3.Known_non_EG(i))/cell2mat(confusion_table_genR3.non_EG_in_model(i));
    precision = cell2mat(confusion_table_genR3.Known_EG(i))/cell2mat(confusion_table_genR3.Predicted_EG(i)); %FP + TP is the same as the predicted genes
    Sensitivity_Generic_R3 = [Sensitivity_Generic_R3,sensitivity];
    Specificity_Generic_R3 = [Specificity_Generic_R3,specificity];
    Precision_Generic_R3 = [Precision_Generic_R3,precision];
end

%% Renal R3
load KO_generic_biomass_CCLE_medium_Renal_input.mat
essential_genes = zeros(size(grRatio_biomass_single,1),size(grRatio_biomass_single,2));
for i=1:size(grRatio_biomass_single,2)
threshold = 0.5;
B = grRatio_biomass_single(:,i) <= threshold;
essential_genes(B,i)=1;
end

EG_final=cell.empty(0, numel(colnames));
for i=1:size(essential_genes,2)
    x =(genelist(find(essential_genes(:,i)==1)));
    [~, ia,ib] = intersect(x,recon3D_genes.ENTREZ); % find the intersection and keep indices
    EG_final(i)= {recon3D_genes.SYMBOL(ib)};
end
confusion_table_Renal = True_False_positives(EG_final,colnames,'Recon3D');
confusion_table_Renal = sortrows(confusion_table_Renal,'Cell_Line','ascend');

Sensitivity_Renal = [];
Specificity_Renal = [];
Precision_Renal = [];
for i = 1:size(confusion_table_Renal,1)
    sensitivity = cell2mat(confusion_table_Renal.Known_EG(i))/cell2mat(confusion_table_Renal.Metabolic_cancer_genes(i));
    specificity = cell2mat(confusion_table_Renal.Known_non_EG(i))/cell2mat(confusion_table_Renal.non_EG_in_model(i));
    precision = cell2mat(confusion_table_Renal.Known_EG(i))/cell2mat(confusion_table_Renal.Predicted_EG(i)); %FP + TP is the same as the predicted genes
    Sensitivity_Renal = [Sensitivity_Renal,sensitivity];
    Specificity_Renal = [Specificity_Renal,specificity];
    Precision_Renal = [Precision_Renal,precision];
end
%% iHsa
%Data 
load dico_iHsa.mat
load intest_skin.mat
colnames = INTESTINE.Properties.VariableNames; %sample names

% I will run it only for the cell lines which have CRISPR data to avoid
% long running times.
load sample_info.mat
load binaryDepScores_broad_database.mat

rownames = binaryDepScores_broad_database.Gene;
crc = sampleinfo((find(contains(sampleinfo.lineage,'colorectal'))),:);

crc_binaryDepScores = [];
binaryDepScores = [];
for x = 2:size(binaryDepScores_broad_database.Properties.VariableNames,2)
    if ismember((erase(binaryDepScores_broad_database.Properties.VariableNames(x),'ACH')),convertStringsToChars(erase(crc.DepMap_ID,'ACH-')))
        loc = find(ismember(convertStringsToChars(erase(crc.DepMap_ID,'ACH-')),(erase(binaryDepScores_broad_database.Properties.VariableNames(x),'ACH'))));
        crc_binaryDepScores = [crc_binaryDepScores,crc.CCLE_Name(loc)];
        binaryDepScores = [binaryDepScores,binaryDepScores_broad_database(:,x)];
    end
end
binaryDepScores.Properties.VariableNames = crc_binaryDepScores;
binaryDepScores.Properties.RowNames = binaryDepScores_broad_database.Gene;


CRISPR = binaryDepScores.Properties.VariableNames;

idx_CRISPR = [];
for i=1:size(CRISPR,2)
    idx = find(ismember(colnames,CRISPR(i)))
    idx_CRISPR = [idx_CRISPR,idx];
end

CRC = [INTESTINE(:,[1:2]),INTESTINE(:,idx_CRISPR)];
fpkmIntestine = table2array(CRC(:,3:end)); %get the numbers
%% 
% Colnames

colnames = (CRC.Properties.VariableNames(3:end))'; %sample names

%In that case no different conditions
conditions = {[1:length(colnames)]};
condition_names = colnames;

%% 
% rownames

rownames = CRC.Name;

%remove NaN entries
rownames(isnan(sum(fpkmIntestine,2)),:) = [];
fpkmIntestine(isnan(sum(fpkmIntestine,2)),:) = [];

%remove 0 entries
rownames(sum(fpkmIntestine,2) == 0,:) = [];
fpkmIntestine(sum(fpkmIntestine,2) == 0,:) = [];
rownames = regexprep(rownames,'\.\d*',''); %remove transcripts

% get symbols for ENSG
[Lia,Locb] = ismember(rownames,iHsa_genes.ENSG);
rownames_SYMBOL = cell(size(rownames));
rownames_SYMBOL(Lia) = iHsa_genes.SYMBOL(Locb(find(Locb)));

%convert empty Symbols to char
for i=1:numel(rownames_SYMBOL)
    if ~ischar(rownames_SYMBOL{i})
        rownames_SYMBOL{i} = char(rownames_SYMBOL{i});
    end
end

% get ENTREZ for ENSG
[Lia,Locb] = ismember(rownames,iHsa_genes.ENSG);
rownames_ENTREZ = cell(size(rownames));
rownames_ENTREZ(Lia) = iHsa_genes.ENTREZ(Locb(find(Locb)));

%convert empty ENTREZ to char
for i=1:numel(rownames_ENTREZ)
    if ~ischar(rownames_ENTREZ{i})
        rownames_ENTREZ{i} = char(rownames_ENTREZ{i});
    end
end

[rownames, rownames_SYMBOL,rownames_ENTREZ]
clear ans i Lia Locb

load KO_generic_biomass_CCLE_medium_iHsa_input.mat
essential_genes = zeros(size(grRatio_biomass_single,1),size(grRatio_biomass_single,2));
for i=1:size(grRatio_biomass_single,2)
threshold = 0.5;
B = grRatio_biomass_single(:,i) <= threshold;
essential_genes(B,i)=1;
end

EG_final=cell.empty(0, numel(colnames));
for i=1:size(essential_genes,2)
    x =(genelist(find(essential_genes(:,i)==1)));
    [~, ia,ib] = intersect(x,iHsa_genes.ENTREZ); % find the intersection and keep indices
    EG_final(i)= {iHsa_genes.SYMBOL(ib)};
end

confusion_table_iHsa = True_False_positives(EG_final,colnames,'iHsa');
confusion_table_iHsa = sortrows(confusion_table_iHsa,'Cell_Line','ascend');

Sensitivity_iHsa_input = [];
Specificity_iHsa_input = [];
Precision_iHsa_input = [];
for i = 1:size(confusion_table_iHsa,1)
    sensitivity = cell2mat(confusion_table_iHsa.Known_EG(i))/cell2mat(confusion_table_iHsa.Metabolic_cancer_genes(i));
    specificity = cell2mat(confusion_table_iHsa.Known_non_EG(i))/cell2mat(confusion_table_iHsa.non_EG_in_model(i));
    precision = cell2mat(confusion_table_iHsa.Known_EG(i))/cell2mat(confusion_table_iHsa.Predicted_EG(i));
    Sensitivity_iHsa_input = [Sensitivity_iHsa_input,sensitivity];
    Specificity_iHsa_input = [Specificity_iHsa_input,specificity];
    Precision_iHsa_input = [Precision_iHsa_input,precision];
end


%% Human1
%Data 
load dico_Human1.mat
load intest_skin.mat
colnames = INTESTINE.Properties.VariableNames; %sample names

% I will run it only for the cell lines which have CRISPR data to avoid
% long running times.
load sample_info.mat
load binaryDepScores_broad_database.mat

rownames = binaryDepScores_broad_database.Gene;
crc = sampleinfo((find(contains(sampleinfo.lineage,'colorectal'))),:);

crc_binaryDepScores = [];
binaryDepScores = [];
for x = 2:size(binaryDepScores_broad_database.Properties.VariableNames,2)
    if ismember((erase(binaryDepScores_broad_database.Properties.VariableNames(x),'ACH')),convertStringsToChars(erase(crc.DepMap_ID,'ACH-')))
        loc = find(ismember(convertStringsToChars(erase(crc.DepMap_ID,'ACH-')),(erase(binaryDepScores_broad_database.Properties.VariableNames(x),'ACH'))));
        crc_binaryDepScores = [crc_binaryDepScores,crc.CCLE_Name(loc)];
        binaryDepScores = [binaryDepScores,binaryDepScores_broad_database(:,x)];
    end
end
binaryDepScores.Properties.VariableNames = crc_binaryDepScores;
binaryDepScores.Properties.RowNames = binaryDepScores_broad_database.Gene;


CRISPR = binaryDepScores.Properties.VariableNames;

idx_CRISPR = [];
for i=1:size(CRISPR,2)
    idx = find(ismember(colnames,CRISPR(i)))
    idx_CRISPR = [idx_CRISPR,idx];
end

CRC = [INTESTINE(:,[1:2]),INTESTINE(:,idx_CRISPR)];
fpkmIntestine = table2array(CRC(:,3:end)); %get the numbers
%% 
% Colnames

colnames = (CRC.Properties.VariableNames(3:end))'; %sample names

%In that case no different conditions
conditions = {[1:length(colnames)]};
condition_names = colnames;

%% 
% rownames

rownames = CRC.Name;

%remove NaN entries
rownames(isnan(sum(fpkmIntestine,2)),:) = [];
fpkmIntestine(isnan(sum(fpkmIntestine,2)),:) = [];

%remove 0 entries
rownames(sum(fpkmIntestine,2) == 0,:) = [];
fpkmIntestine(sum(fpkmIntestine,2) == 0,:) = [];
rownames = regexprep(rownames,'\.\d*',''); %remove transcripts

% get symbols for ENSG
[Lia,Locb] = ismember(rownames,human1_genes.ENSG);
rownames_SYMBOL = cell(size(rownames));
rownames_SYMBOL(Lia) = human1_genes.SYMBOL(Locb(find(Locb)));

%convert empty Symbols to char
for i=1:numel(rownames_SYMBOL)
    if ~ischar(rownames_SYMBOL{i})
        rownames_SYMBOL{i} = char(rownames_SYMBOL{i});
    end
end

% get ENTREZ for ENSG
[Lia,Locb] = ismember(rownames,human1_genes.ENSG);
rownames_ENTREZ = cell(size(rownames));
rownames_ENTREZ(Lia) = human1_genes.ENTREZ(Locb(find(Locb)));

%convert empty ENTREZ to char
for i=1:numel(rownames_ENTREZ)
    if ~ischar(rownames_ENTREZ{i})
        rownames_ENTREZ{i} = char(rownames_ENTREZ{i});
    end
end

[rownames, rownames_SYMBOL,rownames_ENTREZ]
clear ans i Lia Locb

load KO_generic_biomass_human1_input_medium.mat
essential_genes = zeros(size(grRatio_biomass_single,1),size(grRatio_biomass_single,2));
for i=1:size(grRatio_biomass_single,2)
threshold = 0.5;
B = grRatio_biomass_single(:,i) <= threshold;
essential_genes(B,i)=1;
end

EG_final=cell.empty(0, numel(colnames));
for i=1:size(essential_genes,2)
    x =(genelist(find(essential_genes(:,i)==1)));
    [~, ia,ib] = intersect(x,human1_genes.ENSG); % find the intersection and keep indices
    EG_final(i)= {human1_genes.SYMBOL(ib)};
end

confusion_table_Human1 = True_False_positives(EG_final,colnames,'Human1');
confusion_table_Human1 = sortrows(confusion_table_Human1,'Cell_Line','ascend');

Sensitivity_Human1_input = [];
Specificity_Human1_input = [];
Precision_Human1_input = [];
for i = 1:size(confusion_table_Human1,1)
    sensitivity = cell2mat(confusion_table_Human1.Known_EG(i))/cell2mat(confusion_table_Human1.Metabolic_cancer_genes(i));
    specificity = cell2mat(confusion_table_Human1.Known_non_EG(i))/cell2mat(confusion_table_Human1.non_EG_in_model(i));
    precision = cell2mat(confusion_table_Human1.Known_EG(i))/cell2mat(confusion_table_Human1.Predicted_EG(i));
    Sensitivity_Human1_input = [Sensitivity_Human1_input,sensitivity];
    Specificity_Human1_input = [Specificity_Human1_input,specificity];
    Precision_Human1_input = [Precision_Human1_input,precision];
end


%% HMR
%Data 
load dico_hmr.mat
load intest_skin.mat
colnames = INTESTINE.Properties.VariableNames; %sample names

% I will run it only for the cell lines which have CRISPR data to avoid
% long running times.
load sample_info.mat
load binaryDepScores_broad_database.mat

rownames = binaryDepScores_broad_database.Gene;
crc = sampleinfo((find(contains(sampleinfo.lineage,'colorectal'))),:);

crc_binaryDepScores = [];
binaryDepScores = [];
for x = 2:size(binaryDepScores_broad_database.Properties.VariableNames,2)
    if ismember((erase(binaryDepScores_broad_database.Properties.VariableNames(x),'ACH')),convertStringsToChars(erase(crc.DepMap_ID,'ACH-')))
        loc = find(ismember(convertStringsToChars(erase(crc.DepMap_ID,'ACH-')),(erase(binaryDepScores_broad_database.Properties.VariableNames(x),'ACH'))));
        crc_binaryDepScores = [crc_binaryDepScores,crc.CCLE_Name(loc)];
        binaryDepScores = [binaryDepScores,binaryDepScores_broad_database(:,x)];
    end
end
binaryDepScores.Properties.VariableNames = crc_binaryDepScores;
binaryDepScores.Properties.RowNames = binaryDepScores_broad_database.Gene;


CRISPR = binaryDepScores.Properties.VariableNames;

idx_CRISPR = [];
for i=1:size(CRISPR,2)
    idx = find(ismember(colnames,CRISPR(i)))
    idx_CRISPR = [idx_CRISPR,idx];
end

CRC = [INTESTINE(:,[1:2]),INTESTINE(:,idx_CRISPR)];
fpkmIntestine = table2array(CRC(:,3:end)); %get the numbers
%% 
% Colnames

colnames = (CRC.Properties.VariableNames(3:end))'; %sample names

%In that case no different conditions
conditions = {[1:length(colnames)]};
condition_names = colnames;

%% 
% rownames

rownames = CRC.Name;

%remove NaN entries
rownames(isnan(sum(fpkmIntestine,2)),:) = [];
fpkmIntestine(isnan(sum(fpkmIntestine,2)),:) = [];

%remove 0 entries
rownames(sum(fpkmIntestine,2) == 0,:) = [];
fpkmIntestine(sum(fpkmIntestine,2) == 0,:) = [];
rownames = regexprep(rownames,'\.\d*',''); %remove transcripts

% get symbols for ENSG
[Lia,Locb] = ismember(rownames,hmr_genes.ENSG);
rownames_SYMBOL = cell(size(rownames));
rownames_SYMBOL(Lia) = hmr_genes.SYMBOL(Locb(find(Locb)));

%convert empty Symbols to char
for i=1:numel(rownames_SYMBOL)
    if ~ischar(rownames_SYMBOL{i})
        rownames_SYMBOL{i} = char(rownames_SYMBOL{i});
    end
end

% get ENTREZ for ENSG
[Lia,Locb] = ismember(rownames,hmr_genes.ENSG);
rownames_ENTREZ = cell(size(rownames));
rownames_ENTREZ(Lia) = hmr_genes.ENTREZ(Locb(find(Locb)));

%convert empty ENTREZ to char
for i=1:numel(rownames_ENTREZ)
    if ~ischar(rownames_ENTREZ{i})
        rownames_ENTREZ{i} = char(rownames_ENTREZ{i});
    end
end

[rownames, rownames_SYMBOL,rownames_ENTREZ]
clear ans i Lia Locb

load KO_generic_biomass_HMR_input_medium.mat
essential_genes = zeros(size(grRatio_biomass_single,1),size(grRatio_biomass_single,2));
for i=1:size(grRatio_biomass_single,2)
threshold = 0.5;
B = grRatio_biomass_single(:,i) <= threshold;
essential_genes(B,i)=1;
end

EG_final=cell.empty(0, numel(colnames));
for i=1:size(essential_genes,2)
    x =(genelist(find(essential_genes(:,i)==1)));
    [~, ia,ib] = intersect(x,hmr_genes.ENSG); % find the intersection and keep indices
    EG_final(i)= {hmr_genes.SYMBOL(ib)};
end

confusion_table_HMR = True_False_positives(EG_final,colnames,'HMR');
confusion_table_HMR = sortrows(confusion_table_HMR,'Cell_Line','ascend');

Sensitivity_HMR_input = [];
Specificity_HMR_input = [];
Precision_HMR_input = [];
for i = 1:size(confusion_table_HMR,1)
    sensitivity = cell2mat(confusion_table_HMR.Known_EG(i))/cell2mat(confusion_table_HMR.Metabolic_cancer_genes(i));
    specificity = cell2mat(confusion_table_HMR.Known_non_EG(i))/cell2mat(confusion_table_HMR.non_EG_in_model(i));
    precision = cell2mat(confusion_table_HMR.Known_EG(i))/cell2mat(confusion_table_HMR.Predicted_EG(i));
    Sensitivity_HMR_input = [Sensitivity_HMR_input,sensitivity];
    Specificity_HMR_input = [Specificity_HMR_input,specificity];
    Precision_HMR_input = [Precision_HMR_input,precision];
end


%% Generic R2
%Data 
load dico_recon.mat
load intest_skin.mat
colnames = INTESTINE.Properties.VariableNames; %sample names

% I will run it only for the cell lines which have CRISPR data to avoid
% long running times.
load sample_info.mat
load binaryDepScores_broad_database.mat

rownames = binaryDepScores_broad_database.Gene;
crc = sampleinfo((find(contains(sampleinfo.lineage,'colorectal'))),:);

crc_binaryDepScores = [];
binaryDepScores = [];
for x = 2:size(binaryDepScores_broad_database.Properties.VariableNames,2)
    if ismember((erase(binaryDepScores_broad_database.Properties.VariableNames(x),'ACH')),convertStringsToChars(erase(crc.DepMap_ID,'ACH-')))
        loc = find(ismember(convertStringsToChars(erase(crc.DepMap_ID,'ACH-')),(erase(binaryDepScores_broad_database.Properties.VariableNames(x),'ACH'))));
        crc_binaryDepScores = [crc_binaryDepScores,crc.CCLE_Name(loc)];
        binaryDepScores = [binaryDepScores,binaryDepScores_broad_database(:,x)];
    end
end
binaryDepScores.Properties.VariableNames = crc_binaryDepScores;
binaryDepScores.Properties.RowNames = binaryDepScores_broad_database.Gene;


CRISPR = binaryDepScores.Properties.VariableNames;

idx_CRISPR = [];
for i=1:size(CRISPR,2)
    idx = find(ismember(colnames,CRISPR(i)))
    idx_CRISPR = [idx_CRISPR,idx];
end

CRC = [INTESTINE(:,[1:2]),INTESTINE(:,idx_CRISPR)];
fpkmIntestine = table2array(CRC(:,3:end)); %get the numbers
%% 
% Colnames

colnames = (CRC.Properties.VariableNames(3:end))'; %sample names

%In that case no different conditions
conditions = {[1:length(colnames)]};
condition_names = colnames;

%% 
% rownames

rownames = CRC.Name;

%remove NaN entries
rownames(isnan(sum(fpkmIntestine,2)),:) = [];
fpkmIntestine(isnan(sum(fpkmIntestine,2)),:) = [];

%remove 0 entries
rownames(sum(fpkmIntestine,2) == 0,:) = [];
fpkmIntestine(sum(fpkmIntestine,2) == 0,:) = [];
rownames = regexprep(rownames,'\.\d*',''); %remove transcripts

% get symbols for ENSG
[Lia,Locb] = ismember(rownames,dico_RECON.ENSG);
rownames_SYMBOL = cell(size(rownames));
rownames_SYMBOL(Lia) = dico_RECON.SYMBOL(Locb(find(Locb)));

%convert empty Symbols to char
for i=1:numel(rownames_SYMBOL)
    if ~ischar(rownames_SYMBOL{i})
        rownames_SYMBOL{i} = char(rownames_SYMBOL{i});
    end
end

% get ENTREZ for ENSG
[Lia,Locb] = ismember(rownames,dico_RECON.ENSG);
rownames_ENTREZ = cell(size(rownames));
rownames_ENTREZ(Lia) = dico_RECON.ENTREZ(Locb(find(Locb)));

%convert empty ENTREZ to char
for i=1:numel(rownames_ENTREZ)
    if ~ischar(rownames_ENTREZ{i})
        rownames_ENTREZ{i} = char(rownames_ENTREZ{i});
    end
end

[rownames, rownames_SYMBOL,rownames_ENTREZ]
clear ans i Lia Locb

load KO_generic_biomass_CCLE_medium.mat
essential_genes = zeros(size(grRatio_biomass_single,1),size(grRatio_biomass_single,2));
for i=1:size(grRatio_biomass_single,2)
threshold = 0.5;
B = grRatio_biomass_single(:,i) <= threshold;
essential_genes(B,i)=1;
end

EG_final=cell.empty(0, numel(colnames));
for i=1:size(essential_genes,2)
    x =(genelist(find(essential_genes(:,i)==1)));
    [~, ia,ib] = intersect(x,dico_RECON.ENTREZ); % find the intersection and keep indices
    EG_final(i)= {dico_RECON.SYMBOL(ib)};
end

confusion_table_GenR2 = True_False_positives(EG_final,colnames,'Recon');
confusion_table_GenR2 = sortrows(confusion_table_GenR2,'Cell_Line','ascend');

Sensitivity_Generic_R = [];
Specificity_Generic_R = [];
Precision_Generic_R = [];
for i = 1:size(confusion_table_GenR2,1)
    sensitivity = cell2mat(confusion_table_GenR2.Known_EG(i))/cell2mat(confusion_table_GenR2.Metabolic_cancer_genes(i));
    specificity = cell2mat(confusion_table_GenR2.Known_non_EG(i))/cell2mat(confusion_table_GenR2.non_EG_in_model(i));
    precision = cell2mat(confusion_table_GenR2.Known_EG(i))/cell2mat(confusion_table_GenR2.Predicted_EG(i));
    Sensitivity_Generic_R = [Sensitivity_Generic_R,sensitivity];
    Specificity_Generic_R = [Specificity_Generic_R,specificity];
    Precision_Generic_R = [Precision_Generic_R,precision];

end

%% Plot the results
sensitivity_per_cellLine = [];
for i = 1:size(Sensitivity_Generic_R,2)
    generic = Sensitivity_Generic_R(i);
    hmri = Sensitivity_HMR_input(i);
    humani = Sensitivity_Human1_input(i);
    ihsai = Sensitivity_iHsa_input(i);
    mainti = Sensitivity_Maintenance_input(i);
    %renali = Sensitivity_Renal(i);
    noTri = Sensitivity_noTrTr_input(i);
    genericR3 = Sensitivity_Generic_R3(i);
    sensitivity_per_cellLine(i,:) = [generic,hmri,humani,ihsai,mainti,noTri,genericR3];
end

original_order = [{'R_generic_R2'},{'HMR_gen_HMR'},{'Human1_gen_H1'},{'iHsa_gen_iHsa'},{'R3_main_R3'},{'R3_noTrTr_R3'},{'R3_Renal_R3'},{'R_generic_R3'}];

original_order = [{'R_generic_R2'},{'HMR_gen_HMR'},{'Human1_gen_H1'},{'iHsa_gen_iHsa'},{'R3_main_R3'},{'R3_noTrTr_R3'},{'R_generic_R3'}];
sensitivity_per_cellLine(isnan(sensitivity_per_cellLine))=0;

test = array2table(zeros(size(sensitivity_per_cellLine,2),size(sensitivity_per_cellLine,1)));
for i = 1:size(sensitivity_per_cellLine,1)
    test(:,i) = cell2table([num2cell(sensitivity_per_cellLine(i,:))']);
end

imagesc(table2array(test)')
set(gca,'Ytick',[1:18],'YTickLabel',(erase(confusion_table_GenR2.Cell_Line,'_LARGE_INTESTINE')))
set(gca,'Xtick',[1:size(original_order,2)],'XTickLabel',original_order)
colorbar
colormap(altcolor)

Specificity_per_cellLine = [];
for i = 1:size(Specificity_Generic_R,2)
    generic = Specificity_Generic_R(i);
    hmri = Specificity_HMR_input(i);
    humani = Specificity_Human1_input(i);
    ihsai = Specificity_iHsa_input(i);
    mainti = Specificity_Maintenance_input(i);
    noTri = Specificity_noTrTr_input(i);
    renali = Specificity_Renal(i);
    genericR3 = Specificity_Generic_R3(i);
    Specificity_per_cellLine(i,:) = [generic,hmri,humani,ihsai,mainti,noTri,renali,genericR3];
end

original_order = [{'R_generic_R2'},{'HMR_gen_HMR'},{'Human1_gen_H1'},{'iHsa_gen_iHsa'},{'R3_main_R3'},{'R3_noTrTr_R3'},{'R3_Renal_R3'},{'R_generic_R3'}];
Specificity_per_cellLine(isnan(Specificity_per_cellLine))=0.95; %Remove NaN to avoid future problems

test = array2table(zeros(size(Specificity_per_cellLine,2),size(Specificity_per_cellLine,1)));
for i = 1:size(Specificity_per_cellLine,1)
    test(:,i) = cell2table([num2cell(Specificity_per_cellLine(i,:))']);
end
figure
imagesc(table2array(test)')
set(gca,'Ytick',[1:18],'YTickLabel',(erase(confusion_table_GenericB_R.Cell_Line,'_LARGE_INTESTINE')))
set(gca,'Xtick',[1:size(original_order,2)],'XTickLabel',original_order)
colorbar
colormap(altcolor)

precision_per_cellLine = [];
for i = 1:size(Precision_Generic_R,2)
    generic = Precision_Generic_R(i);
    hmri = Precision_HMR_input(i);
    humani = Precision_Human1_input(i);
    ihsai = Precision_iHsa_input(i);
    mainti = Precision_Maintenance_input(i);
    noTri = Precision_noTrTr_input(i);
    renali = Precision_Renal(i);
    genericR3 = Precision_Generic_R3(i);
    precision_per_cellLine(i,:) = [generic,hmri,humani,ihsai,mainti,noTri,renali,genericR3];
end

original_order = [{'R_generic_R2'},{'HMR_gen_HMR'},{'Human1_gen_H1'},{'iHsa_gen_iHsa'},{'R3_main_R3'},{'R3_noTrTr_R3'},{'R3_Renal_R3'},{'R_generic_R3'}];
precision_per_cellLine(isnan(precision_per_cellLine))=0;

test = array2table(zeros(size(precision_per_cellLine,2),size(precision_per_cellLine,1)));
for i = 1:size(precision_per_cellLine,1)
    test(:,i) = cell2table([num2cell(precision_per_cellLine(i,:))']);
end

imagesc(table2array(test)')
set(gca,'Ytick',[1:18],'YTickLabel',(erase(confusion_table_GenericB_R.Cell_Line,'_LARGE_INTESTINE')))
set(gca,'Xtick',[1:size(original_order,2)],'XTickLabel',original_order)
colorbar
colormap(altcolor)