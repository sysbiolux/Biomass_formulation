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

load KO_iHsa_biomass_CCLE_medium.mat 
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
    idx = find(ismember(colnames,CRISPR(i)));
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

[rownames, rownames_SYMBOL,rownames_ENTREZ];
clear ans i Lia Locb

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

confusion_table_iHsa_Input = True_False_positives(EG_final,colnames,'Recon');
confusion_table_iHsa_Input = sortrows(confusion_table_iHsa_Input,'Cell_Line','ascend');

Sensitivity_iHsa_input = [];
Specificity_iHsa_input = [];
Precision_iHsa_input = [];
for i = 1:size(confusion_table_iHsa_Input,1)
    sensitivity = cell2mat(confusion_table_iHsa_Input.Known_EG(i))/cell2mat(confusion_table_iHsa_Input.Metabolic_cancer_genes(i));
    specificity = cell2mat(confusion_table_iHsa_Input.Known_non_EG(i))/cell2mat(confusion_table_iHsa_Input.non_EG_in_model(i));
    precision = cell2mat(confusion_table_iHsa_Input.Known_EG(i))/cell2mat(confusion_table_iHsa_Input.Predicted_EG(i));
    Sensitivity_iHsa_input = [Sensitivity_iHsa_input,sensitivity];
    Specificity_iHsa_input = [Specificity_iHsa_input,specificity];
    Precision_iHsa_input = [Precision_iHsa_input,precision];
end

%% Sensitivity and specificity analysis Human1
load KO_human1_biomass_CCLE_medium.mat
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

confusion_table_Human1_Input = True_False_positives(EG_final,colnames,'Recon');
confusion_table_Human1_Input = sortrows(confusion_table_Human1_Input,'Cell_Line','ascend');

Sensitivity_Human1_input = [];
Specificity_Human1_input = [];
Precision_Human1_input = [];
for i = 1:size(confusion_table_Human1_Input,1)
    sensitivity = cell2mat(confusion_table_Human1_Input.Known_EG(i))/cell2mat(confusion_table_Human1_Input.Metabolic_cancer_genes(i));
    specificity = cell2mat(confusion_table_Human1_Input.Known_non_EG(i))/cell2mat(confusion_table_Human1_Input.non_EG_in_model(i));
    precision = cell2mat(confusion_table_Human1_Input.Known_EG(i))/cell2mat(confusion_table_Human1_Input.Predicted_EG(i));
    Sensitivity_Human1_input = [Sensitivity_Human1_input,sensitivity];
    Specificity_Human1_input = [Specificity_Human1_input,specificity];
    Precision_Human1_input = [Precision_Human1_input,precision];
end

%% Sensitivity and specificity analysis Renal
load KO_Renal_biomass_CCLE_medium.mat
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

confusion_table_Renal_Input = True_False_positives(EG_final,colnames,'Recon');
confusion_table_Renal_Input = sortrows(confusion_table_Renal_Input,'Cell_Line','ascend');

Sensitivity_Renal_input = [];
Specificity_Renal_input = [];
Precision_Renal_input = [];
for i = 1:size(confusion_table_Renal_Input,1)
    sensitivity = cell2mat(confusion_table_Renal_Input.Known_EG(i))/cell2mat(confusion_table_Renal_Input.Metabolic_cancer_genes(i));
    specificity = cell2mat(confusion_table_Renal_Input.Known_non_EG(i))/cell2mat(confusion_table_Renal_Input.non_EG_in_model(i));
    precision = cell2mat(confusion_table_Renal_Input.Known_EG(i))/cell2mat(confusion_table_Renal_Input.Predicted_EG(i));
    Sensitivity_Renal_input = [Sensitivity_Renal_input,sensitivity];
    Specificity_Renal_input = [Specificity_Renal_input,specificity];
    Precision_Renal_input = [Precision_Renal_input, precision];
end

%% Sensitivity and specificity analysis maintenance
load KO_maint_biomass_CCLE_medium.mat
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

confusion_table_Maintenance_Input = True_False_positives(EG_final,colnames,'Recon');
confusion_table_Maintenance_Input = sortrows(confusion_table_Maintenance_Input,'Cell_Line','ascend');

Sensitivity_Maintenance_input = [];
Specificity_Maintenance_input = [];
Precision_Maintenance_input = [];
for i = 1:size(confusion_table_Maintenance_Input,1)
    sensitivity = cell2mat(confusion_table_Maintenance_Input.Known_EG(i))/cell2mat(confusion_table_Maintenance_Input.Metabolic_cancer_genes(i));
    specificity = cell2mat(confusion_table_Maintenance_Input.Known_non_EG(i))/cell2mat(confusion_table_Maintenance_Input.non_EG_in_model(i));
    precision = cell2mat(confusion_table_Maintenance_Input.Known_EG(i))/cell2mat(confusion_table_Maintenance_Input.Predicted_EG(i));
    Sensitivity_Maintenance_input = [Sensitivity_Maintenance_input,sensitivity];
    Specificity_Maintenance_input = [Specificity_Maintenance_input,specificity];
    Precision_Maintenance_input = [Precision_Maintenance_input,precision];
end

%% Sensitivity and specificity analysis noTrTr
load KO_noTrTr_biomass_CCLE_medium.mat
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

confusion_table_noTrTr_Input = True_False_positives(EG_final,colnames,'Recon');
confusion_table_noTrTr_Input = sortrows(confusion_table_noTrTr_Input,'Cell_Line','ascend');

Sensitivity_noTrTr_input = [];
Specificity_noTrTr_input = [];
Precision_noTrTr_input = [];
for i = 1:size(confusion_table_noTrTr_Input,1)
    sensitivity = cell2mat(confusion_table_noTrTr_Input.Known_EG(i))/cell2mat(confusion_table_noTrTr_Input.Metabolic_cancer_genes(i));
    specificity = cell2mat(confusion_table_noTrTr_Input.Known_non_EG(i))/cell2mat(confusion_table_noTrTr_Input.non_EG_in_model(i));
    precision = cell2mat(confusion_table_noTrTr_Input.Known_EG(i))/cell2mat(confusion_table_noTrTr_Input.Predicted_EG(i));
    Sensitivity_noTrTr_input = [Sensitivity_noTrTr_input,sensitivity];
    Specificity_noTrTr_input = [Specificity_noTrTr_input,specificity];
    Precision_noTrTr_input = [Precision_noTrTr_input,precision];
end

%% Sensitivity and specificity analysis Generic
load KO_HMR_biomass_CCLE_medium.mat
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

confusion_table_HMR_Input = True_False_positives(EG_final,colnames,'Recon');
confusion_table_HMR_Input = sortrows(confusion_table_HMR_Input,'Cell_Line','ascend');

Sensitivity_HMR_input = [];
Specificity_HMR_input = [];
Precision_HMR_input = [];
for i = 1:size(confusion_table_HMR_Input,1)
    sensitivity = cell2mat(confusion_table_HMR_Input.Known_EG(i))/cell2mat(confusion_table_HMR_Input.Metabolic_cancer_genes(i));
    specificity = cell2mat(confusion_table_HMR_Input.Known_non_EG(i))/cell2mat(confusion_table_HMR_Input.non_EG_in_model(i));
    precision = cell2mat(confusion_table_HMR_Input.Known_EG(i))/cell2mat(confusion_table_HMR_Input.Predicted_EG(i));
    Sensitivity_HMR_input = [Sensitivity_HMR_input,sensitivity];
    Specificity_HMR_input = [Specificity_HMR_input,specificity];
    Precision_HMR_input = [Precision_HMR_input,precision];
end

%% Sensitivity and specificity analysis Generic
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

confusion_table_GenericB_R = True_False_positives(EG_final,colnames,'Recon');
confusion_table_GenericB_R = sortrows(confusion_table_GenericB_R,'Cell_Line','ascend');

Sensitivity_Generic_R = [];
Specificity_Generic_R = [];
Precision_Generic_R = [];
for i = 1:size(confusion_table_GenericB_R,1)
    sensitivity = cell2mat(confusion_table_GenericB_R.Known_EG(i))/cell2mat(confusion_table_GenericB_R.Metabolic_cancer_genes(i));
    specificity = cell2mat(confusion_table_GenericB_R.Known_non_EG(i))/cell2mat(confusion_table_GenericB_R.non_EG_in_model(i));
    precision = cell2mat(confusion_table_GenericB_R.Known_EG(i))/cell2mat(confusion_table_GenericB_R.Predicted_EG(i));
    Sensitivity_Generic_R = [Sensitivity_Generic_R,sensitivity];
    Specificity_Generic_R = [Specificity_Generic_R,specificity];
    Precision_Generic_R = [Precision_Generic_R,precision];
end

sensitivity_per_cellLine = [];
for i = 1:size(Sensitivity_Generic_R,2)
    generic = Sensitivity_Generic_R(i);
    hmri = Sensitivity_HMR_input(i);
    humani = Sensitivity_Human1_input(i);
    ihsai = Sensitivity_iHsa_input(i);
    mainti = Sensitivity_Maintenance_input(i);
    noTri = Sensitivity_noTrTr_input(i);
    renali = Sensitivity_Renal_input(i);
    sensitivity_per_cellLine(i,:) = [generic,hmri,humani,ihsai,mainti,noTri,renali];
end

original_order = [{'R_generic_R2'},{'HMR_gen_R2'},{'Human1_gen_R2'},{'iHsa_gen_R2'},{'R3_main_R2'},{'R3_noTrTr_R2'},{'HMR_renal_R2'}];
sensitivity_per_cellLine(isnan(sensitivity_per_cellLine))=0;

test = array2table(zeros(size(sensitivity_per_cellLine,2),size(sensitivity_per_cellLine,1)));
for i = 1:size(sensitivity_per_cellLine,1)
    test(:,i) = cell2table([num2cell(sensitivity_per_cellLine(i,:))']);
end

imagesc(table2array(test)')
set(gca,'Ytick',[1:18],'YTickLabel',(erase(confusion_table_GenericB_R.Cell_Line,'_LARGE_INTESTINE')))
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
    renali = Specificity_Renal_input(i);
    Specificity_per_cellLine(i,:) = [generic,hmri,humani,ihsai,mainti,noTri,renali];
end

original_order = [{'R_generic_R2'},{'HMR_gen_R2'},{'Human1_gen_R2'},{'iHsa_gen_R2'},{'R3_main_R2'},{'R3_noTrTr_R2'},{'HMR_renal_R2'}];
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

Precision_per_cellLine = [];
for i = 1:size(Precision_Generic_R,2)
    generic = Precision_Generic_R(i);
    hmri = Precision_HMR_input(i);
    humani = Precision_Human1_input(i);
    ihsai = Precision_iHsa_input(i);
    mainti = Precision_Maintenance_input(i);
    noTri = Precision_noTrTr_input(i);
    renali = Precision_Renal_input(i);
    Precision_per_cellLine(i,:) = [generic,hmri,humani,ihsai,mainti,noTri,renali];
end

original_order = [{'R_generic_R2'},{'HMR_gen_R2'},{'Human1_gen_R2'},{'iHsa_gen_R2'},{'R3_main_R2'},{'R3_noTrTr_R2'},{'HMR_renal_R2'}];
Precision_per_cellLine(isnan(Precision_per_cellLine))=0.95; %Remove NaN to avoid future problems

test = array2table(zeros(size(Precision_per_cellLine,2),size(Precision_per_cellLine,1)));
for i = 1:size(Precision_per_cellLine,1)
    test(:,i) = cell2table([num2cell(Precision_per_cellLine(i,:))']);
end
figure
imagesc(table2array(test)')
set(gca,'Ytick',[1:18],'YTickLabel',(erase(confusion_table_GenericB_R.Cell_Line,'_LARGE_INTESTINE')))
set(gca,'Xtick',[1:size(original_order,2)],'XTickLabel',original_order)
colorbar
colormap(altcolor)