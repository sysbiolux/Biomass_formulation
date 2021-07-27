%% Initialise

cd('/mnt/irisgpfs/users/mamoscardo/') %working folder
addpath(genpath(pwd)) %add all files and folders in the working folder to matlab
addpath(genpath('/mnt/irisgpfs/users/mamoscardo/cobratoolbox/')) %path to rFASTCORMICS
addpath(genpath('/mnt/irisgpfs/apps/resif/data/production/v1.2-20191021/default/software/math/CPLEX/12.9-foss-2019a/cplex/matlab/x86-64_linux/')) %add cplex path

altcolor= [255 255 255;255 204 204; 255 153 153; 255 102 102; 255 51 51;...
    255 0 0; 204 0 0; 152 0 0; 102 0 0;  51 0 0; 0 0 0]/255; %some red color
set(0,'defaultTextInterpreter','none')
set(groot, 'DefaultAxesTickLabelInterpreter', 'none')
set(groot, 'defaultLegendInterpreter','none');

changeCobraSolver('ibm_cplex')
load dico_Human1.mat
%% Import and clean Data

%data = readtable('CCLE_RNAseq_genes_rpkm_20180929.txt');

% colnames = data.Properties.VariableNames
% INTESTINE_idx = ~cellfun(@isempty,regexpi(colnames, 'Intest'))
% SKIN_idx =  ~cellfun(@isempty,regexpi(colnames, 'skin'))
% 
% colnames(INTESTINE_idx) %check colnames
% 
% INTESTINE = [data(:,[1:2]),data(:,INTESTINE_idx)];
% SKIN = [data(:,[1:2]),data(:,SKIN_idx)];
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

discretized = discretize_FPKM(fpkmIntestine, colnames);

%% Create models
%% Initialise
load Consistent_model_Human1.mat
consistent_model.grRules = strrep(consistent_model.grRules,'or','|');
consistent_model.grRules = strrep(consistent_model.grRules,'and','&');
consistent_model.rules = consistent_model.grRules;
% Obtain rules compatible with rFASTCORMICS
for i = 1:size(consistent_model.rules,1)
    b = char(consistent_model.rules(i));
    C = strtrim(strsplit(b,{' '}));
    for a = 1:size(C,2)
        if (contains(C(a),'ENSG'))
            [z,matches] = strsplit(char(C(a)),{'(',')'});
            idx = find(ismember(consistent_model.genes,z(~cellfun(@isempty,strsplit(char(C(a)),{'(',')'})))));
            if contains(matches,'(')
                C{a} = strcat(char(matches),'x(',num2str(idx),')'); 
            elseif contains(matches,')')
                C{a} = strcat('x(',num2str(idx),')',char(matches)); 
            elseif isempty(matches)
                C{a} = strcat('x(',num2str(idx),')'); 
            end
        else 
            C{a} = char(C{a});
        end
    end
    consistent_model.rules(i) = {strjoin(C)};
end

Cmodel = consistent_model;

load dico_medium_culture.mat
A_final = zeros(numel(Cmodel.rxns),numel(colnames))
col = cell2table(colnames);
medium_CellLines = cellstr(medium_CellLines{:, :});

for a = 1:size(Cmodel.subSystems,1)
     Cmodel.subSystems{a}=char(Cmodel.subSystems{a}); %I did this becaus otherwise it was not a cell of character arrays
end

for i=1:numel(col)   % need to define medium for the cells used here
    used_medium=medium_CellLines(i,2);
    if ~cellfun(@isempty, strfind(used_medium,'RPMI'))
         load ('RPMI_human1_CompF.mat')
    elseif ~cellfun(@isempty, strfind(used_medium,'DMEM'))
         load ('DMEM_Human1_CompF.mat')
    elseif ~cellfun(@isempty, strfind(used_medium,'EMEM'))
         load ('EMEM_Human1_CompF.mat')
    elseif ~cellfun(@isempty, strfind(used_medium,'F12K'))
         load ('F12_Human1_CompF.mat')
    elseif ~cellfun(@isempty, strfind(used_medium,'McCoy'))
         load ('McCoys_Human1_CompF.mat')
    elseif ~cellfun(@isempty, strfind(used_medium,'Leibovitz'))
         load ('L15_Human1_CompF.mat')
     elseif ~cellfun(@isempty, strfind(used_medium,'MEM'))
         load ('MEM_Human1_CompF.mat')
    end
    

epsilon = 1e-4;
consensus_proportion = 0.9; %gene has to be expressed in 90% of the cases in order to be included. Only relevant if you want to create one model from different samples
already_mapped_tag = 0;

unpenalizedSystems = {'Transport, endoplasmic reticular';
    'Transport, extracellular';
    'Transport, golgi apparatus';
    'Transport, mitochondrial';
    'Transport, peroxisomal';
    'Transport, lysosomal';
    'Transport, nuclear'};
unpenalized = Cmodel.rxns(ismember(Cmodel.subSystems,unpenalizedSystems));
not_medium_constrained = 'HMR_10041';

optional_settings.unpenalized = unpenalized;
optional_settings.func = {'biomass_human'};

optional_settings.not_medium_constrained = not_medium_constrained;
optional_settings.medium = medium;



%% Single models
[~, A] = fastcormics_RNAseq(Cmodel, discretized(:,i), rownames, human1_genes, ...
            already_mapped_tag, consensus_proportion, epsilon, optional_settings);
A_final(A,i)= 1;
    
delete *.log
poolobj = gcp('nocreate');
delete(poolobj);

models_keep_single = A_final;

save('models_keep_single_CCLE_Human1_input_biomass_medium','models_keep_single')
end


%load models_keep_single_CCLE_Human1_input_biomass_medium
optional_settings.func = {'biomass_human'};
model = Cmodel;
%load KO_generic_biomass_human1_input_medium
for i=1:size(models_keep_single,2)
        ind = find(ismember(model.rxns,optional_settings.func{1}));
        model_out = removeRxns(model,model.rxns(setdiff(1:numel(model.rxns),find(models_keep_single(:,i))))); % create model based on active reactions
        model_out = changeObjective(model_out,model.rxns(ind)); % set objective function
        %model_out = removeUnusedGenes(model_out);

        [grRatio, grRateKO, grRateWT, hasEffect, delRxns, fluxSolution, genelist] = singleGeneDeletionTamara(model_out,'FBA',[],0,1);
        
        grRatio_biomass_single(:,i)    = grRatio;
        grRateKO_biomass_single(:,i)   = grRateKO;
        grRateWT_biomass_single(:,i)   = grRateWT;
        hasEffect_biomass_single(:,i)  = hasEffect;
        delRxns_biomass_single(:,i)    = delRxns;
        fluxSolution_biomass_single{i} = fluxSolution;
        
        save('KO_generic_biomass_human1_input_medium','grRatio_biomass_single','grRateKO_biomass_single','grRateWT_biomass_single','genelist')

        end
    clear i ind grRatio grRateKO grRateWT hasEffect delRxns fluxSolution
