altcolor= [255 255 255;255 204 204; 255 153 153; 255 102 102; 255 51 51;...
    255 0 0; 204 0 0; 152 0 0; 102 0 0;  51 0 0; 0 0 0]/255; %some red color
set(0,'defaultTextInterpreter','none')
set(groot, 'DefaultAxesTickLabelInterpreter', 'none')
set(groot, 'defaultLegendInterpreter','none');

load Recon204.mat
Cmodel = CmodelR204; clear CmodelR204
%% Data 
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


%% iHsa
load KO_iHsa_biomass_CCLE_medium.mat 
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
%% Drugs targeting specific cell lines
load dico_drugs_final_recon_filtered.mat

grRatio = grRatio_biomass_single(sum(grRatio_biomass_single,2)<=(numel(colnames)-0.001),:);
gene_names = genelist(sum(grRatio_biomass_single,2)<=(numel(colnames)-0.001),:);

idx = []
for i = 1:size(grRatio,1)
    if any(grRatio(i,:)<= 0.5)
        idx = [idx,i];
    end
end
grRatio_keep = grRatio(idx,:);
EG_names = gene_names(idx);

%get pathways for each gene, gene -> reactions(s) -> pathway(s)
for i=1:numel(Cmodel.genes)
    pathways_per_gene{i,1} = strjoin(unique(Cmodel.subSystems(find(Cmodel.rxnGeneMat(:,i)))),'|');
end
[gene_info, ia, ~]  = unique(regexprep(Cmodel.genes,'\..*',''));
pathways_per_gene   = pathways_per_gene(ia);
[~, Lia, Locb]      = intersect(gene_info, dico_RECON.ENTREZ);
gene_info(Lia,2)    = dico_RECON.SYMBOL(Locb);
gene_info(:,3)      = pathways_per_gene;

clear i pathways_per_gene ia Lia Locb

gene_drug_interactions = dico_final_recon_filtered(ismember(dico_final_recon_filtered.ENTREZ,gene_info(:,1)),:);
gene_drug_interactions(cellfun(@isempty, gene_drug_interactions.DBID),:) = [];
gene_drug_interactions = unique(gene_drug_interactions);
gene_drug_interactions = gene_drug_interactions(~cellfun(@isempty,(strfind(gene_drug_interactions.Action, 'inhibitor'))),:);

drugs = unique(gene_drug_interactions.DrugName);

drugsXsample = zeros(numel(drugs),size(grRatio_keep,2));

for i=1:numel(drugs)
    idx = find(ismember(gene_drug_interactions.DrugName,drugs(i)));
    idx2 = find(ismember(EG_names(:,1), gene_drug_interactions.ENTREZ(idx)));
    
    drugsXsample(i,:) = sum(grRatio_keep(idx2,:),1)/numel(idx2);
end

%remove NaN entries
non_effective = ~isnan(sum(drugsXsample,2));
drugsXsample(isnan(sum(drugsXsample,2)),:) = [];
drugs = drugs(non_effective);
table = [drugs;drugsXsample];

drugs_iHsa = drugs;
drugsXsample_iHsa = drugsXsample;

%% Human1
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
%% Drugs targeting specific cell lines
load dico_drugs_final_recon_filtered.mat

grRatio = grRatio_biomass_single(sum(grRatio_biomass_single,2)<=(numel(colnames)-0.001),:);
gene_names = genelist(sum(grRatio_biomass_single,2)<=(numel(colnames)-0.001),:);

idx = []
for i = 1:size(grRatio,1)
    if any(grRatio(i,:)<= 0.5)
        idx = [idx,i];
    end
end
grRatio_keep = grRatio(idx,:);
EG_names = gene_names(idx);

%get pathways for each gene, gene -> reactions(s) -> pathway(s)
for i=1:numel(Cmodel.genes)
    pathways_per_gene{i,1} = strjoin(unique(Cmodel.subSystems(find(Cmodel.rxnGeneMat(:,i)))),'|');
end
[gene_info, ia, ~]  = unique(regexprep(Cmodel.genes,'\..*',''));
pathways_per_gene   = pathways_per_gene(ia);
[~, Lia, Locb]      = intersect(gene_info, dico_RECON.ENTREZ);
gene_info(Lia,2)    = dico_RECON.SYMBOL(Locb);
gene_info(:,3)      = pathways_per_gene;

clear i pathways_per_gene ia Lia Locb

gene_drug_interactions = dico_final_recon_filtered(ismember(dico_final_recon_filtered.ENTREZ,gene_info(:,1)),:);
gene_drug_interactions(cellfun(@isempty, gene_drug_interactions.DBID),:) = [];
gene_drug_interactions = unique(gene_drug_interactions);
gene_drug_interactions = gene_drug_interactions(~cellfun(@isempty,(strfind(gene_drug_interactions.Action, 'inhibitor'))),:);

drugs = unique(gene_drug_interactions.DrugName);

drugsXsample = zeros(numel(drugs),size(grRatio_keep,2));

for i=1:numel(drugs)
    idx = find(ismember(gene_drug_interactions.DrugName,drugs(i)));
    idx2 = find(ismember(EG_names(:,1), gene_drug_interactions.ENTREZ(idx)));
    
    drugsXsample(i,:) = sum(grRatio_keep(idx2,:),1)/numel(idx2);
end

%remove NaN entries
non_effective = ~isnan(sum(drugsXsample,2));
drugsXsample(isnan(sum(drugsXsample,2)),:) = [];
drugs = drugs(non_effective);
table = [drugs;drugsXsample];

drugs_Human1 = drugs;
drugsXsample_Human1 = drugsXsample;

%% Generic
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
%% Drugs targeting specific cell lines
load dico_drugs_final_recon_filtered.mat

grRatio = grRatio_biomass_single(sum(grRatio_biomass_single,2)<=(numel(colnames)-0.001),:);
gene_names = genelist(sum(grRatio_biomass_single,2)<=(numel(colnames)-0.001),:);

idx = []
for i = 1:size(grRatio,1)
    if any(grRatio(i,:)<= 0.5)
        idx = [idx,i];
    end
end
grRatio_keep = grRatio(idx,:);
EG_names = gene_names(idx);

%get pathways for each gene, gene -> reactions(s) -> pathway(s)
for i=1:numel(Cmodel.genes)
    pathways_per_gene{i,1} = strjoin(unique(Cmodel.subSystems(find(Cmodel.rxnGeneMat(:,i)))),'|');
end
[gene_info, ia, ~]  = unique(regexprep(Cmodel.genes,'\..*',''));
pathways_per_gene   = pathways_per_gene(ia);
[~, Lia, Locb]      = intersect(gene_info, dico_RECON.ENTREZ);
gene_info(Lia,2)    = dico_RECON.SYMBOL(Locb);
gene_info(:,3)      = pathways_per_gene;

clear i pathways_per_gene ia Lia Locb

gene_drug_interactions = dico_final_recon_filtered(ismember(dico_final_recon_filtered.ENTREZ,gene_info(:,1)),:);
gene_drug_interactions(cellfun(@isempty, gene_drug_interactions.DBID),:) = [];
gene_drug_interactions = unique(gene_drug_interactions);
gene_drug_interactions = gene_drug_interactions(~cellfun(@isempty,(strfind(gene_drug_interactions.Action, 'inhibitor'))),:);

drugs = unique(gene_drug_interactions.DrugName);

drugsXsample = zeros(numel(drugs),size(grRatio_keep,2));

for i=1:numel(drugs)
    idx = find(ismember(gene_drug_interactions.DrugName,drugs(i)));
    idx2 = find(ismember(EG_names(:,1), gene_drug_interactions.ENTREZ(idx)));
    
    drugsXsample(i,:) = sum(grRatio_keep(idx2,:),1)/numel(idx2);
end

%remove NaN entries
non_effective = ~isnan(sum(drugsXsample,2));
drugsXsample(isnan(sum(drugsXsample,2)),:) = [];
drugs = drugs(non_effective);
table = [drugs;drugsXsample];

drugs_generic = drugs;
drugsXsample_generic = drugsXsample;

%% Maint
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
%% Drugs targeting specific cell lines
load dico_drugs_final_recon_filtered.mat

grRatio = grRatio_biomass_single(sum(grRatio_biomass_single,2)<=(numel(colnames)-0.001),:);
gene_names = genelist(sum(grRatio_biomass_single,2)<=(numel(colnames)-0.001),:);

idx = []
for i = 1:size(grRatio,1)
    if any(grRatio(i,:)<= 0.5)
        idx = [idx,i];
    end
end
grRatio_keep = grRatio(idx,:);
EG_names = gene_names(idx);

%get pathways for each gene, gene -> reactions(s) -> pathway(s)
for i=1:numel(Cmodel.genes)
    pathways_per_gene{i,1} = strjoin(unique(Cmodel.subSystems(find(Cmodel.rxnGeneMat(:,i)))),'|');
end
[gene_info, ia, ~]  = unique(regexprep(Cmodel.genes,'\..*',''));
pathways_per_gene   = pathways_per_gene(ia);
[~, Lia, Locb]      = intersect(gene_info, dico_RECON.ENTREZ);
gene_info(Lia,2)    = dico_RECON.SYMBOL(Locb);
gene_info(:,3)      = pathways_per_gene;

clear i pathways_per_gene ia Lia Locb

gene_drug_interactions = dico_final_recon_filtered(ismember(dico_final_recon_filtered.ENTREZ,gene_info(:,1)),:);
gene_drug_interactions(cellfun(@isempty, gene_drug_interactions.DBID),:) = [];
gene_drug_interactions = unique(gene_drug_interactions);
gene_drug_interactions = gene_drug_interactions(~cellfun(@isempty,(strfind(gene_drug_interactions.Action, 'inhibitor'))),:);

drugs = unique(gene_drug_interactions.DrugName);

drugsXsample = zeros(numel(drugs),size(grRatio_keep,2));

for i=1:numel(drugs)
    idx = find(ismember(gene_drug_interactions.DrugName,drugs(i)));
    idx2 = find(ismember(EG_names(:,1), gene_drug_interactions.ENTREZ(idx)));
    
    drugsXsample(i,:) = sum(grRatio_keep(idx2,:),1)/numel(idx2);
end

%remove NaN entries
non_effective = ~isnan(sum(drugsXsample,2));
drugsXsample(isnan(sum(drugsXsample,2)),:) = [];
drugs = drugs(non_effective);
table = [drugs;drugsXsample];

drugs_maint = drugs;
drugsXsample_maint = drugsXsample;

%% noTrTr
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
%% Drugs targeting specific cell lines
load dico_drugs_final_recon_filtered.mat

grRatio = grRatio_biomass_single(sum(grRatio_biomass_single,2)<=(numel(colnames)-0.001),:);
gene_names = genelist(sum(grRatio_biomass_single,2)<=(numel(colnames)-0.001),:);

idx = []
for i = 1:size(grRatio,1)
    if any(grRatio(i,:)<= 0.5)
        idx = [idx,i];
    end
end
grRatio_keep = grRatio(idx,:);
EG_names = gene_names(idx);

%get pathways for each gene, gene -> reactions(s) -> pathway(s)
for i=1:numel(Cmodel.genes)
    pathways_per_gene{i,1} = strjoin(unique(Cmodel.subSystems(find(Cmodel.rxnGeneMat(:,i)))),'|');
end
[gene_info, ia, ~]  = unique(regexprep(Cmodel.genes,'\..*',''));
pathways_per_gene   = pathways_per_gene(ia);
[~, Lia, Locb]      = intersect(gene_info, dico_RECON.ENTREZ);
gene_info(Lia,2)    = dico_RECON.SYMBOL(Locb);
gene_info(:,3)      = pathways_per_gene;

clear i pathways_per_gene ia Lia Locb

gene_drug_interactions = dico_final_recon_filtered(ismember(dico_final_recon_filtered.ENTREZ,gene_info(:,1)),:);
gene_drug_interactions(cellfun(@isempty, gene_drug_interactions.DBID),:) = [];
gene_drug_interactions = unique(gene_drug_interactions);
gene_drug_interactions = gene_drug_interactions(~cellfun(@isempty,(strfind(gene_drug_interactions.Action, 'inhibitor'))),:);

drugs = unique(gene_drug_interactions.DrugName);

drugsXsample = zeros(numel(drugs),size(grRatio_keep,2));

for i=1:numel(drugs)
    idx = find(ismember(gene_drug_interactions.DrugName,drugs(i)));
    idx2 = find(ismember(EG_names(:,1), gene_drug_interactions.ENTREZ(idx)));
    
    drugsXsample(i,:) = sum(grRatio_keep(idx2,:),1)/numel(idx2);
end

%remove NaN entries
non_effective = ~isnan(sum(drugsXsample,2));
drugsXsample(isnan(sum(drugsXsample,2)),:) = [];
drugs = drugs(non_effective);
table = [drugs;drugsXsample];

drugs_noTrTr = drugs;
drugsXsample_noTrTr = drugsXsample;

%% HMR
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
%% Drugs targeting specific cell lines
load dico_drugs_final_recon_filtered.mat

grRatio = grRatio_biomass_single(sum(grRatio_biomass_single,2)<=(numel(colnames)-0.001),:);
gene_names = genelist(sum(grRatio_biomass_single,2)<=(numel(colnames)-0.001),:);

idx = []
for i = 1:size(grRatio,1)
    if any(grRatio(i,:)<= 0.5)
        idx = [idx,i];
    end
end
grRatio_keep = grRatio(idx,:);
EG_names = gene_names(idx);

%get pathways for each gene, gene -> reactions(s) -> pathway(s)
for i=1:numel(Cmodel.genes)
    pathways_per_gene{i,1} = strjoin(unique(Cmodel.subSystems(find(Cmodel.rxnGeneMat(:,i)))),'|');
end
[gene_info, ia, ~]  = unique(regexprep(Cmodel.genes,'\..*',''));
pathways_per_gene   = pathways_per_gene(ia);
[~, Lia, Locb]      = intersect(gene_info, dico_RECON.ENTREZ);
gene_info(Lia,2)    = dico_RECON.SYMBOL(Locb);
gene_info(:,3)      = pathways_per_gene;

clear i pathways_per_gene ia Lia Locb

gene_drug_interactions = dico_final_recon_filtered(ismember(dico_final_recon_filtered.ENTREZ,gene_info(:,1)),:);
gene_drug_interactions(cellfun(@isempty, gene_drug_interactions.DBID),:) = [];
gene_drug_interactions = unique(gene_drug_interactions);
gene_drug_interactions = gene_drug_interactions(~cellfun(@isempty,(strfind(gene_drug_interactions.Action, 'inhibitor'))),:);

drugs = unique(gene_drug_interactions.DrugName);

drugsXsample = zeros(numel(drugs),size(grRatio_keep,2));

for i=1:numel(drugs)
    idx = find(ismember(gene_drug_interactions.DrugName,drugs(i)));
    idx2 = find(ismember(EG_names(:,1), gene_drug_interactions.ENTREZ(idx)));
    
    drugsXsample(i,:) = sum(grRatio_keep(idx2,:),1)/numel(idx2);
end

%remove NaN entries
non_effective = ~isnan(sum(drugsXsample,2));
drugsXsample(isnan(sum(drugsXsample,2)),:) = [];
drugs = drugs(non_effective);
table = [drugs;drugsXsample];

drugs_HMR = drugs;
drugsXsample_HMR = drugsXsample;

%% Renal
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
%% Drugs targeting specific cell lines
load dico_drugs_final_recon_filtered.mat

grRatio = grRatio_biomass_single(sum(grRatio_biomass_single,2)<=(numel(colnames)-0.001),:);
gene_names = genelist(sum(grRatio_biomass_single,2)<=(numel(colnames)-0.001),:);

idx = []
for i = 1:size(grRatio,1)
    if any(grRatio(i,:)<= 0.5)
        idx = [idx,i];
    end
end
grRatio_keep = grRatio(idx,:);
EG_names = gene_names(idx);

%get pathways for each gene, gene -> reactions(s) -> pathway(s)
for i=1:numel(Cmodel.genes)
    pathways_per_gene{i,1} = strjoin(unique(Cmodel.subSystems(find(Cmodel.rxnGeneMat(:,i)))),'|');
end
[gene_info, ia, ~]  = unique(regexprep(Cmodel.genes,'\..*',''));
pathways_per_gene   = pathways_per_gene(ia);
[~, Lia, Locb]      = intersect(gene_info, dico_RECON.ENTREZ);
gene_info(Lia,2)    = dico_RECON.SYMBOL(Locb);
gene_info(:,3)      = pathways_per_gene;

clear i pathways_per_gene ia Lia Locb

gene_drug_interactions = dico_final_recon_filtered(ismember(dico_final_recon_filtered.ENTREZ,gene_info(:,1)),:);
gene_drug_interactions(cellfun(@isempty, gene_drug_interactions.DBID),:) = [];
gene_drug_interactions = unique(gene_drug_interactions);
gene_drug_interactions = gene_drug_interactions(~cellfun(@isempty,(strfind(gene_drug_interactions.Action, 'inhibitor'))),:);

drugs = unique(gene_drug_interactions.DrugName);

drugsXsample = zeros(numel(drugs),size(grRatio_keep,2));

for i=1:numel(drugs)
    idx = find(ismember(gene_drug_interactions.DrugName,drugs(i)));
    idx2 = find(ismember(EG_names(:,1), gene_drug_interactions.ENTREZ(idx)));
    
    drugsXsample(i,:) = sum(grRatio_keep(idx2,:),1)/numel(idx2);
end

%remove NaN entries
non_effective = ~isnan(sum(drugsXsample,2));
drugsXsample(isnan(sum(drugsXsample,2)),:) = [];
drugs = drugs(non_effective);
table = [drugs;drugsXsample];

drugs_renal = drugs;
drugsXsample_renal = drugsXsample;

%% Plot the previous results
drugs = [{drugs_generic},{drugs_HMR},{drugs_Human1},{drugs_iHsa},{drugs_maint},{drugs_noTrTr},{drugs_renal}];

overlap=[];
for r = 1:numel(drugs)
    for c = 1:numel(drugs)
        overlap(r,c) = sum(ismember(drugs{r},drugs{c}));
    end
end
names = {'R_generic_R2','HMR_gen_R2','Human1_gen_R2','iHsa_gen_R2','R3_main_R2','R3_noTrTr_R2','HMR_renal_R2'};
heatmap(overlap,...
    'Colormap',altcolor,...
    'XDisplayLabels',names,...
    "YDisplayLabels",names)
% Determine how many predicted drugs are cancer drugs
drugs = [{drugs_generic},{drugs_HMR},{drugs_Human1},{drugs_iHsa},{drugs_maint},{drugs_noTrTr},{drugs_renal}];
cancer = [];
for i = 1:size(drugs,2)
    enrichment = DrugEnrichments(drugs{i});
    cancer(i) = cell2mat(enrichment.DrugList_cancer(find(ismember(enrichment.Database_website,'SEERRx_all'))));
end

%% Part b of the figure. Using home models 
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

load dico_drugs_final_Recon3D.mat
load Consistent_model_Recon3D.mat
Cmodel = consistent_model;

%Since we had issues with cell line 15 and we set all the column to zero, we will remove this column to get the essential genes. 
grRatio_biomass_single(:,15) = [];
grRatio = grRatio_biomass_single(sum(grRatio_biomass_single,2)<=(size(grRatio_biomass_single,2)-0.001),:); 
gene_names = genelist(sum(grRatio_biomass_single,2)<=(size(grRatio_biomass_single,2)-0.001),:);

idx = []
for i = 1:size(grRatio,1)
    if any(grRatio(i,:)<= 0.5)
        idx = [idx,i];
    end
end
grRatio_keep = grRatio(idx,:);
EG_names = gene_names(idx);


[gene_info, ia, ~]  = unique(regexprep(Cmodel.genes,'\..*',''));
[~, Lia, Locb]      = intersect(gene_info, recon3D_genes.ENTREZ);
gene_info(Lia,2)    = recon3D_genes.SYMBOL(Locb);

clear i pathways_per_gene ia Lia Locb

gene_drug_interactions = dico_Recon3D_final(ismember(dico_Recon3D_final.ENTREZ,gene_info(:,1)),:);
gene_drug_interactions(cellfun(@isempty, gene_drug_interactions.DBID),:) = [];
gene_drug_interactions = unique(gene_drug_interactions);
gene_drug_interactions = gene_drug_interactions(~cellfun(@isempty,(strfind(gene_drug_interactions.Action, 'inhibitor'))),:);

drugs = unique(gene_drug_interactions.DrugName);

drugsXsample = zeros(numel(drugs),size(grRatio_keep,2));

for i=1:numel(drugs)
    idx = find(ismember(gene_drug_interactions.DrugName,drugs(i)));
    idx2 = find(ismember(EG_names(:,1), gene_drug_interactions.ENTREZ(idx)));
    
    drugsXsample(i,:) = sum(grRatio_keep(idx2,:),1)/numel(idx2);
end

%remove NaN entries
non_effective = ~isnan(sum(drugsXsample,2));
drugsXsample(isnan(sum(drugsXsample,2)),:) = [];
drugs = drugs(non_effective);
table = [drugs;drugsXsample];

drugs_maint_input = drugs;
drugsXsample_maint_input = drugsXsample;

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

%Since we had issues with cell line 15 and we set all the column to zero, we will remove this column to get the essential genes. 
grRatio_biomass_single(:,15) = [];
grRatio = grRatio_biomass_single(sum(grRatio_biomass_single,2)<=(size(grRatio_biomass_single,2)-0.001),:);
gene_names = genelist(sum(grRatio_biomass_single,2)<=(size(grRatio_biomass_single,2)-0.001),:);

idx = []
for i = 1:size(grRatio,1)
    if any(grRatio(i,:)<= 0.5)
        idx = [idx,i];
    end
end
grRatio_keep = grRatio(idx,:);
EG_names = gene_names(idx);


[gene_info, ia, ~]  = unique(regexprep(Cmodel.genes,'\..*',''));
[~, Lia, Locb]      = intersect(gene_info, recon3D_genes.ENTREZ);
gene_info(Lia,2)    = recon3D_genes.SYMBOL(Locb);

clear i pathways_per_gene ia Lia Locb

gene_drug_interactions = dico_Recon3D_final(ismember(dico_Recon3D_final.ENTREZ,gene_info(:,1)),:);
gene_drug_interactions(cellfun(@isempty, gene_drug_interactions.DBID),:) = [];
gene_drug_interactions = unique(gene_drug_interactions);
gene_drug_interactions = gene_drug_interactions(~cellfun(@isempty,(strfind(gene_drug_interactions.Action, 'inhibitor'))),:);

drugs = unique(gene_drug_interactions.DrugName);

drugsXsample = zeros(numel(drugs),size(grRatio_keep,2));

for i=1:numel(drugs)
    idx = find(ismember(gene_drug_interactions.DrugName,drugs(i)));
    idx2 = find(ismember(EG_names(:,1), gene_drug_interactions.ENTREZ(idx)));
    
    drugsXsample(i,:) = sum(grRatio_keep(idx2,:),1)/numel(idx2);
end

%remove NaN entries
non_effective = ~isnan(sum(drugsXsample,2));
drugsXsample(isnan(sum(drugsXsample,2)),:) = [];
drugs = drugs(non_effective);
table = [drugs;drugsXsample];

drugs_noTrTr_input = drugs;
drugsXsample_noTrTr_input = drugsXsample;

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

grRatio_biomass_single(:,15) = [];
grRatio = grRatio_biomass_single(sum(grRatio_biomass_single,2)<=(size(grRatio_biomass_single,2)-0.001),:);
gene_names = genelist(sum(grRatio_biomass_single,2)<=(size(grRatio_biomass_single,2)-0.001),:);

idx = []
for i = 1:size(grRatio,1)
    if any(grRatio(i,:)<= 0.5)
        idx = [idx,i];
    end
end
grRatio_keep = grRatio(idx,:);
EG_names = gene_names(idx);


[gene_info, ia, ~]  = unique(regexprep(Cmodel.genes,'\..*',''));
[~, Lia, Locb]      = intersect(gene_info, recon3D_genes.ENTREZ);
gene_info(Lia,2)    = recon3D_genes.SYMBOL(Locb);

clear i pathways_per_gene ia Lia Locb

gene_drug_interactions = dico_Recon3D_final(ismember(dico_Recon3D_final.ENTREZ,gene_info(:,1)),:);
gene_drug_interactions(cellfun(@isempty, gene_drug_interactions.DBID),:) = [];
gene_drug_interactions = unique(gene_drug_interactions);
gene_drug_interactions = gene_drug_interactions(~cellfun(@isempty,(strfind(gene_drug_interactions.Action, 'inhibitor'))),:);

drugs = unique(gene_drug_interactions.DrugName);

drugsXsample = zeros(numel(drugs),size(grRatio_keep,2));

for i=1:numel(drugs)
    idx = find(ismember(gene_drug_interactions.DrugName,drugs(i)));
    idx2 = find(ismember(EG_names(:,1), gene_drug_interactions.ENTREZ(idx)));
    
    drugsXsample(i,:) = sum(grRatio_keep(idx2,:),1)/numel(idx2);
end

%remove NaN entries
non_effective = ~isnan(sum(drugsXsample,2));
drugsXsample(isnan(sum(drugsXsample,2)),:) = [];
drugs = drugs(non_effective);
table = [drugs;drugsXsample];

drugs_GenericR3_input = drugs;
drugsXsample_GenericR3_input = drugsXsample;
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

grRatio_biomass_single(:,15) = [];
grRatio = grRatio_biomass_single(sum(grRatio_biomass_single,2)<=(size(grRatio_biomass_single,2)-0.001),:);
gene_names = genelist(sum(grRatio_biomass_single,2)<=(size(grRatio_biomass_single,2)-0.001),:);

idx = []
for i = 1:size(grRatio,1)
    if any(grRatio(i,:)<= 0.5)
        idx = [idx,i];
    end
end
grRatio_keep = grRatio(idx,:);
EG_names = gene_names(idx);


[gene_info, ia, ~]  = unique(regexprep(Cmodel.genes,'\..*',''));
[~, Lia, Locb]      = intersect(gene_info, recon3D_genes.ENTREZ);
gene_info(Lia,2)    = recon3D_genes.SYMBOL(Locb);

clear i pathways_per_gene ia Lia Locb

gene_drug_interactions = dico_Recon3D_final(ismember(dico_Recon3D_final.ENTREZ,gene_info(:,1)),:);
gene_drug_interactions(cellfun(@isempty, gene_drug_interactions.DBID),:) = [];
gene_drug_interactions = unique(gene_drug_interactions);
gene_drug_interactions = gene_drug_interactions(~cellfun(@isempty,(strfind(gene_drug_interactions.Action, 'inhibitor'))),:);

drugs = unique(gene_drug_interactions.DrugName);

drugsXsample = zeros(numel(drugs),size(grRatio_keep,2));

for i=1:numel(drugs)
    idx = find(ismember(gene_drug_interactions.DrugName,drugs(i)));
    idx2 = find(ismember(EG_names(:,1), gene_drug_interactions.ENTREZ(idx)));
    
    drugsXsample(i,:) = sum(grRatio_keep(idx2,:),1)/numel(idx2);
end

%remove NaN entries
non_effective = ~isnan(sum(drugsXsample,2));
drugsXsample(isnan(sum(drugsXsample,2)),:) = [];
drugs = drugs(non_effective);
table = [drugs;drugsXsample];

drugs_Renal_input = drugs;
drugsXsample_Renal_input = drugsXsample;
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

load dico_drugs_final_iHsa.mat
load iHsa_consistent.mat
Cmodel = consistent_model;

grRatio = grRatio_biomass_single(sum(grRatio_biomass_single,2)<=(numel(colnames)-0.001),:);
gene_names = genelist(sum(grRatio_biomass_single,2)<=(numel(colnames)-0.001),:);

idx = []
for i = 1:size(grRatio,1)
    if any(grRatio(i,:)<= 0.5)
        idx = [idx,i];
    end
end
grRatio_keep = grRatio(idx,:);
EG_names = gene_names(idx);


[gene_info, ia, ~]  = unique(regexprep(Cmodel.genes,'\..*',''));
[~, Lia, Locb]      = intersect(gene_info,iHsa_genes.ENTREZ);
gene_info(Lia,2)    = iHsa_genes.SYMBOL(Locb);

clear i pathways_per_gene ia Lia Locb

gene_drug_interactions = dico_iHsa_final(ismember(dico_iHsa_final.ENTREZ,gene_info(:,1)),:);
gene_drug_interactions(cellfun(@isempty, gene_drug_interactions.DBID),:) = [];
gene_drug_interactions = unique(gene_drug_interactions);
gene_drug_interactions = gene_drug_interactions(~cellfun(@isempty,(strfind(gene_drug_interactions.Action, 'inhibitor'))),:);

drugs = unique(gene_drug_interactions.DrugName);

drugsXsample = zeros(numel(drugs),size(grRatio_keep,2));

for i=1:numel(drugs)
    idx = find(ismember(gene_drug_interactions.DrugName,drugs(i)));
    idx2 = find(ismember(EG_names(:,1), gene_drug_interactions.ENTREZ(idx)));
    
    drugsXsample(i,:) = sum(grRatio_keep(idx2,:),1)/numel(idx2);
end

%remove NaN entries
non_effective = ~isnan(sum(drugsXsample,2));
drugsXsample(isnan(sum(drugsXsample,2)),:) = [];
drugs = drugs(non_effective);
table = [drugs;drugsXsample];

drugs_iHsa_input = drugs;
drugsXsample_iHsa_input = drugsXsample;

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

load dico_drugs_final_human1.mat
load Consistent_model_Human1.mat
Cmodel = consistent_model;

grRatio = grRatio_biomass_single(sum(grRatio_biomass_single,2)<=(numel(colnames)-0.001),:);
gene_names = genelist(sum(grRatio_biomass_single,2)<=(numel(colnames)-0.001),:);

idx = []
for i = 1:size(grRatio,1)
    if any(grRatio(i,:)<= 0.5)
        idx = [idx,i];
    end
end
grRatio_keep = grRatio(idx,:);
EG_names = gene_names(idx);


[gene_info, ia, ~]  = unique(regexprep(Cmodel.genes,'\..*',''));
[~, Lia, Locb]      = intersect(gene_info,human1_genes.ENSG);
gene_info(Lia,2)    = human1_genes.SYMBOL(Locb);

clear i pathways_per_gene ia Lia Locb

gene_drug_interactions = dico_Human1_final(ismember(dico_Human1_final.ENSG,gene_info(:,1)),:);
gene_drug_interactions(cellfun(@isempty, gene_drug_interactions.DBID),:) = [];
gene_drug_interactions = unique(gene_drug_interactions);
gene_drug_interactions = gene_drug_interactions(~cellfun(@isempty,(strfind(gene_drug_interactions.Action, 'inhibitor'))),:);

drugs = unique(gene_drug_interactions.DrugName);

drugsXsample = zeros(numel(drugs),size(grRatio_keep,2));

for i=1:numel(drugs)
    idx = find(ismember(gene_drug_interactions.DrugName,drugs(i)));
    idx2 = find(ismember(EG_names(:,1), gene_drug_interactions.ENSG(idx)));
    
    drugsXsample(i,:) = sum(grRatio_keep(idx2,:),1)/numel(idx2);
end

%remove NaN entries
non_effective = ~isnan(sum(drugsXsample,2));
drugsXsample(isnan(sum(drugsXsample,2)),:) = [];
drugs = drugs(non_effective);
table = [drugs;drugsXsample];

drugs_Human1_input = drugs;
drugsXsample_Human1_input = drugsXsample;


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
    [~, ia,ib] = intersect(x,hmr_genes.ENTREZ); % find the intersection and keep indices
    EG_final(i)= {hmr_genes.SYMBOL(ib)};
end

load dico_drugs_final_HMR.mat
load Consistent_model_HMR20.mat
Cmodel = consistent_model;

grRatio = grRatio_biomass_single(sum(grRatio_biomass_single,2)<=(numel(colnames)-0.001),:);
gene_names = genelist(sum(grRatio_biomass_single,2)<=(numel(colnames)-0.001),:);

idx = []
for i = 1:size(grRatio,1)
    if any(grRatio(i,:)<= 0.5)
        idx = [idx,i];
    end
end
grRatio_keep = grRatio(idx,:);
EG_names = gene_names(idx);


[gene_info, ia, ~]  = unique(regexprep(Cmodel.genes,'\..*',''));
[~, Lia, Locb]      = intersect(gene_info,hmr_genes.ENSG);
gene_info(Lia,2)    = hmr_genes.SYMBOL(Locb);

clear i pathways_per_gene ia Lia Locb

gene_drug_interactions = dico_HMR_final(ismember(dico_HMR_final.ENSG,gene_info(:,1)),:);
gene_drug_interactions(cellfun(@isempty, gene_drug_interactions.DBID),:) = [];
gene_drug_interactions = unique(gene_drug_interactions);
gene_drug_interactions = gene_drug_interactions(~cellfun(@isempty,(strfind(gene_drug_interactions.Action, 'inhibitor'))),:);

drugs = unique(gene_drug_interactions.DrugName);

drugsXsample = zeros(numel(drugs),size(grRatio_keep,2));

for i=1:numel(drugs)
    idx = find(ismember(gene_drug_interactions.DrugName,drugs(i)));
    idx2 = find(ismember(EG_names(:,1), gene_drug_interactions.ENSG(idx)));
    
    drugsXsample(i,:) = sum(grRatio_keep(idx2,:),1)/numel(idx2);
end

%remove NaN entries
non_effective = ~isnan(sum(drugsXsample,2));
drugsXsample(isnan(sum(drugsXsample,2)),:) = [];
drugs = drugs(non_effective);
table = [drugs;drugsXsample];

drugs_HMR_input = drugs;
drugsXsample_HMR_input = drugsXsample;

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

load dico_drugs_final_recon_filtered.mat
load Recon204.mat
Cmodel = CmodelR204;

grRatio = grRatio_biomass_single(sum(grRatio_biomass_single,2)<=(numel(colnames)-0.001),:);
gene_names = genelist(sum(grRatio_biomass_single,2)<=(numel(colnames)-0.001),:);

idx = []
for i = 1:size(grRatio,1)
    if any(grRatio(i,:)<= 0.5)
        idx = [idx,i];
    end
end
grRatio_keep = grRatio(idx,:);
EG_names = gene_names(idx);


[gene_info, ia, ~]  = unique(regexprep(Cmodel.genes,'\..*',''));
[~, Lia, Locb]      = intersect(gene_info,dico_RECON.ENTREZ);
gene_info(Lia,2)    = dico_RECON.SYMBOL(Locb);

clear i pathways_per_gene ia Lia Locb

gene_drug_interactions = dico_final_recon_filtered(ismember(dico_final_recon_filtered.ENTREZ,gene_info(:,1)),:);
gene_drug_interactions(cellfun(@isempty, gene_drug_interactions.DBID),:) = [];
gene_drug_interactions = unique(gene_drug_interactions);
gene_drug_interactions = gene_drug_interactions(~cellfun(@isempty,(strfind(gene_drug_interactions.Action, 'inhibitor'))),:);

drugs = unique(gene_drug_interactions.DrugName);

drugsXsample = zeros(numel(drugs),size(grRatio_keep,2));

for i=1:numel(drugs)
    idx = find(ismember(gene_drug_interactions.DrugName,drugs(i)));
    idx2 = find(ismember(EG_names(:,1), gene_drug_interactions.ENTREZ(idx)));
    
    drugsXsample(i,:) = sum(grRatio_keep(idx2,:),1)/numel(idx2);
end

%remove NaN entries
non_effective = ~isnan(sum(drugsXsample,2));
drugsXsample(isnan(sum(drugsXsample,2)),:) = [];
drugs = drugs(non_effective);
table = [drugs;drugsXsample];

drugs_GenericR2_input = drugs;
drugsXsample_GenericR2_input = drugsXsample;

%% Plot the previous results
drugs = [{drugs_GenericR2_input},{drugs_HMR_input},{drugs_Human1_input},{drugs_iHsa_input},{drugs_maint_input},{drugs_noTrTr_input},{drugs_Renal_input},{drugs_GenericR3_input}];

overlap=[];
for r = 1:numel(drugs)
    for c = 1:numel(drugs)
        overlap(r,c) = sum(ismember(drugs{r},drugs{c}));
    end
end
names = {'R_generic_R2','HMR_gen_HMR','Human1_gen_H1','iHsa_gen_iHsa','R3_main_R3','R3_noTrTr_R3','R3_Renal_R3','R_generic_R3'};
heatmap(overlap,...
    'Colormap',altcolor,...
    'XDisplayLabels',names,...
    "YDisplayLabels",names)
% Determine how many predicted drugs are cancer drugs
drugs = [{drugs_GenericR2_input},{drugs_HMR_input},{drugs_Human1_input},{drugs_iHsa_input},{drugs_maint_input},{drugs_noTrTr_input},{drugs_Renal_input},{drugs_GenericR3_input}];
cancer = [];
for i = 1:size(drugs,2)
    enrichment = DrugEnrichments(drugs{i});
    cancer(i) = cell2mat(enrichment.DrugList_cancer(find(ismember(enrichment.Database_website,'SEERRx_all'))));
end
