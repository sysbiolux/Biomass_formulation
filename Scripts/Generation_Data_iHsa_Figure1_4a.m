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
load dico_recon.mat

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

 %% Create models Recon2 with Recon3 biomass maintenance reaction
%% Introduce the biomass reaction
%% %% Introduce iHsa biomass reaction
load Recon204.mat
Cmodel = CmodelR204;
% We remove the original biomass reaction in Recon2 to avoid overlapping.
reactionR2 = find(~cellfun(@isempty,regexpi(Cmodel.rxns, 'biomass')));
Cmodel_no_Biomass = removeRxns(Cmodel,Cmodel.rxns(reactionR2));

model = addReaction(Cmodel_no_Biomass,'clpn_consume','reactionFormula','clpn_hs[c] -> ');
%model = addReaction(model,'biomass_reaction','reactionFormula','33 dna[c] + 50 rna[c] + 10 m90191[c] + 120 m90192[c] + m90193[c] + 87 m90194[c] + 42 m90195[c] + 160 m90199[c] -> m99999[c]');
model = addReaction(model,'biomass_reaction','reactionFormula','33 dna[c] + 50 rna[c] + 10 m90191[c] + 120 m90192[c] + m90193[c] + 87 m90194[c] + 42 m90195[c] + 160 m90199[c] -> m99999[c]');

model = addReaction(model,'biomass_out','reactionFormula','m99999[c] -> ');

%DNA
model = addReaction(model,'dna','reactionFormula','0.3 datp[c] + 0.2 dctp[c] + 0.2 dgtp[c] + 0.3 dttp[c] -> dna[c] + ppi[c]');

%RNA 
model = addReaction(model,'rna','reactionFormula','0.18 atp[c] + 0.3 ctp[c] + 0.34 gtp[c] + 0.18 utp[c] -> ppi[c] + rna[c]');

%Aminoacids
model = addReaction(model,'aminoacids','reactionFormula','0.05 ile_L[c] + 0.16 leu_L[c] + 0.14 lys_L[c] + 0.03 met_L[c] + 0.05 phe_L[c] + 0.3 thr_L[c] + 0.01 trp_L[c] + 0.08 tyr_L[c] + 0.18 val_L[c] -> m90191[c]');

model = addReaction(model,'non_aminoacids','reactionFormula','0.15 ala_L[c] + 0.01 arg_L[c] + 0.01 asn_L[c] + 0.3 asp_L[c] + 0.01 cys_L[c] + 0.08 glu_L[c] + 0.25 gln_L[c] + 0.11 gly[c] + 0.03 his_L[c] + 0.01 pro_L[c] + 0.04 ser_L[c] -> m90192[c]');

%Bile acid
model = addReaction(model,'bile_acid','reactionFormula','0.0007 m00745[c] + 0.0004 cholate[c] + 0.32 dgchol[c] + 0.163 gchola[c] + 0.17 m01989[c] + 0.0147 m02000[c] + 0.015 m02004[c] + 0.17 tdchola[c] + 0.067 tchola[c] + 0.0002 m90173[c] + 0.004 m02966[c] + 0.008 m02965[c] + 0.067 m02964[c] -> m90193[c]');

model = addReaction(model,'m00745','reactionFormula','cholate[c] + h2o[c] + nadp[c] -> m00745[c] + h[c] + nadph[c] + o2[c]');
model = addReaction(model,'m01989','reactionFormula','m00745[c] + gly[c] + h[c] -> m01989[c] + h2o[c]');
model = addReaction(model,'m02000','reactionFormula','atp[c] + gly[c] + lithocholate[c] -> amp[c] + m02000[c] + h[c] + ppi[c]');
model = addReaction(model,'lithocholate','reactionFormula','3B_hydroxy_5_cholestenal[c] + h2o[c] + nadp[c] -> 2 h[c] + lithocholate[c] + nadph[c]');
model = addReaction(model,'3B_hydroxy_5_cholestenal','reactionFormula','xol27oh[m] + nadp[c] -> 3B_hydroxy_5_cholestenal[c] + h[c] + nadph[c]');
model = addReaction(model,'m02004','reactionFormula','gly[c] + ursodeoxycholylCoA[c] -> coa[c] + m02004[c]');
model = addReaction(model,'ursodeoxycholylCoA','reactionFormula','atp[c] + coa[c] + ursodeoxycholate[c] -> amp[c] + ppi[c] + ursodeoxycholylCoA[c]');
model = addReaction(model,'ursodeoxycholate','reactionFormula','7ketolitocholate[c] + h[c] + nadph[c] -> nadp[c] + ursodeoxycholate[c]');
model = addReaction(model,'7ketolitocholate','reactionFormula','C02528[c] + h[c] + nadph[c] + o2[c] -> 7ketolitocholate[c] + 2 h2o[c] + nadp[c]');
model = addReaction(model,'m90173','reactionFormula','gly[c] + m90185[c] -> coa[c] + m90173[c]');
model = addReaction(model,'m90185','reactionFormula',' -> m90185[c]');
model = addReaction(model,'m90173_out','reactionFormula','m90173[c] -> ');
model = addReaction(model,'m02966','reactionFormula','taur[c] + ursodeoxycholylCoA[c] -> coa[c] + m02966[c]');
model = addReaction(model,'m02965','reactionFormula','atp[c] + lithocholate[c] + taur[c] -> amp[c] + h[c] + ppi[c] + m02965[c]');
model = addReaction(model,'m02964','reactionFormula','taur[c] + deoxycholoylCoA[c] -> coa[c] + h[c] + m02964[c]');
model = addReaction(model,'deoxycholoylCoA','reactionFormula','m00745[c] + atp[c] + coa[c] -> amp[c] + deoxycholoylCoA[c] + ppi[c]');

model = addReaction(model,'m90194','reactionFormula','udpg[c] -> udp[c] + m90194[c]');

model = addReaction(model,'m02958','reactionFormula','dag_hs[c] + acoa[c] -> coa[c] + m02958[c]');
model = addReaction(model,'sphings','reactionFormula','crm_hs[c] + coa[c] + h[c] -> sphings[c] + acoa[c]');
model = addReaction(model,'m90195','reactionFormula','0.003 xolest2_hs[c] + 0.374 pchol_hs[c] + 0.25 pe_hs[c] + 0.106 pail_hs[c] + 0.082 ps_hs[c] + 0.081 sphmyln_hs[c] + 0.104 m02958[c] -> m90195[c]');

%Metabolite synthesis
%model = addReaction(model,'met_synthesis','reactionFormula','0.001 bhb[c] + 0.001 3mob[c] + 0.001 3pg[c] + 0.002 4abut[c] + 0.004 adp[c] + 0.004 akg[c] + 0.001 amp[c] + 0.038 ascb_L[c] + 0.088 ala_B[c] + 0.001 betaD_glucose6p[c] + 0.001 glyb[c] + 0.002 camp[c] + 0.001 carn[c] + 0.001 cdp[c] + 0.002 chol[c] + 0.001 cis_aconiate[c] + 0.015 cit[c] + 0.005 citr_L[c] + 0.001 cmp[c] + 0.001 coa[c] + 0.004 creat[c] + 0.001 crtn[c] + 0.002 D_gluconic_acid[c] + 0.001 dhap[c] + 0.001 e4p[c] + 0.001 fol[c] + 0.001 fdp[c] + 0.001 f6p[c] + 0.003 fum[c] + 0.001 gdp[c] + 0.001 g1p[c] + 0.001 gmp[c] + 0.49 gthrd[c] + 0.034 gthox[c] + 0.004 hcys_L[c] + 0.001 icit[c] + 0.002 cyst_L[c] + 0.047 lac_L[c] + 0.01 mal_L[c] + 0.022 nad[c] + 0.002 nadp[c] + 0.008 orn[c] + 0.001 pep[c] + 0.002 prpp[c] + 0.001 ptrc[c] + 0.001 pydx5p[c] + 0.001 r5p[c] + 0.001 ru5p_D[c] + 0.002 ahcys[c] + 0.003 amet[c] + 0.001 s7p[c] + 0.002 glyc3p[c] + 0.001 spmd[c] + 0.001 sprm[c] + 0.002 succ[c] + 0.16 taur[c] + 0.001 thf[c] + 0.001 4hpro_LT[c] + 0.003 q10h2[m] + 0.005 q10[m] + 0.001 udp[c] + 0.001 ump[c] -> m90199[c]');
model = addReaction(model,'met_synthesis','reactionFormula','0.001 bhb[c] + 0.001 3mob[c] + 0.001 3pg[c] + 0.002 4abut[c] + 0.004 adp[c] + 0.004 akg[c] + 0.001 amp[c] + 0.038 ascb_L[c] + 0.088 ala_B[c] + 0.001 betaD_glucose6p[c] + 0.001 glyb[c] + 0.002 camp[c] + 0.001 carn[c] + 0.001 cdp[c] + 0.002 chol[c] + 0.001 HC00342[c] + 0.015 cit[c] + 0.005 citr_L[c] + 0.001 cmp[c] + 0.001 coa[c] + 0.004 creat[c] + 0.001 crtn[c] + 0.002 glcn[c] + 0.001 dhap[c] + 0.001 e4p[c] + 0.001 fol[c] + 0.001 fdp[c] + 0.001 f6p[c] + 0.003 fum[c] + 0.001 gdp[c] + 0.001 g1p[c] + 0.001 gmp[c] + 0.49 gthrd[c] + 0.034 gthox[c] + 0.004 hcys_L[c] + 0.001 icit[c] + 0.002 cyst_L[c] + 0.047 lac_L[c] + 0.01 mal_L[c] + 0.022 nad[c] + 0.002 nadp[c] + 0.008 orn[c] + 0.001 pep[c] + 0.002 prpp[c] + 0.001 ptrc[c] + 0.001 pydx5p[c] + 0.001 r5p[c] + 0.001 ru5p_D[c] + 0.002 ahcys[c] + 0.003 amet[c] + 0.001 s7p[c] + 0.002 glyc3p[c] + 0.001 spmd[c] + 0.001 sprm[c] + 0.002 succ[c] + 0.16 taur[c] + 0.001 thf[c] + 0.001 4hpro_LT[c] + 0.003 q10h2[m] + 0.005 q10[m] + 0.001 udp[c] + 0.001 ump[c] -> m90199[c]'); %falta cis aconitate (m01580c)
model = addReaction(model,'betaD_glucose6p','reactionFormula','atp[c] + betaD_glucose6[c] -> betaD_glucose6p[c] + adp[c] + h[c]'); 
model = addReaction(model,'betaD_glucose6','reactionFormula','glc_D[c] -> betaD_glucose6[c]');
model = addReaction(model,'creat','reactionFormula','pcreat[c] -> crtn[c] + ppi[c]');
model = addReaction(model,'glcn','reactionFormula','6pgc[c] + adp[c] + h[c] -> glcn[c] + atp[c]');
model = addReaction(model,'nadp','reactionFormula',' -> nadp[c]');

%Extra reactions needed

% To set the reversible reactions
model.rev=zeros(numel(model.lb),1);
model.rev(model.lb<0)=1;
model.lb=model.lb*1000;
model.ub=model.ub*1000;

r=find(model.ub==0 & model.lb<0);
model.S(:,r)=-model.S(:,r);
tmp=model.lb(r);
model.lb(r)=-model.ub(r);
model.ub(r)=-tmp;
A = fastcc_4_rfastcormics(model,1e-4,0);

Cmodel = model;
discretized = discretize_FPKM(fpkmIntestine, colnames);
load dico_medium_culture.mat
A_final = zeros(numel(Cmodel.rxns),numel(colnames));
col = cell2table(colnames);
medium_CellLines = cellstr(medium_CellLines{:, :});

for a = 1:size(Cmodel.subSystems,1)
     Cmodel.subSystems{a}=char(Cmodel.subSystems{a}); %I did this becaus otherwise it was not a cell of character arrays
end

for i=1:numel(col)   % need to define medium for the cells used here
    used_medium=medium_CellLines(i,2);
    if ~cellfun(@isempty, strfind(used_medium,'RPMI'))
         load ('RPMI_charstrippedFBS_DHT_MTA_content.mat')
    elseif ~cellfun(@isempty, strfind(used_medium,'DMEM'))
         load ('DMEM_medium_iHsa.mat')
    elseif ~cellfun(@isempty, strfind(used_medium,'EMEM'))
         load ('EMEM_medium.mat')
    elseif ~cellfun(@isempty, strfind(used_medium,'F12K'))
         load ('F12K_medium.mat')
    elseif ~cellfun(@isempty, strfind(used_medium,'McCoy'))
         load ('McCoys_Medium_iHsa.mat')
    elseif ~cellfun(@isempty, strfind(used_medium,'Leibovitz'))
         load ('L15_medium_iHsa.mat')
    elseif ~cellfun(@isempty, strfind(used_medium,'MEM'))
         load ('MEM_medium.mat')
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
not_medium_constrained = 'EX_tag_hs(e)';

optional_settings.unpenalized = unpenalized;
optional_settings.func = {'biomass_reaction'};

optional_settings.not_medium_constrained = not_medium_constrained;
optional_settings.medium = medium;



% Single models
[~, A] = fastcormics_RNAseq(Cmodel, discretized(:,i), rownames, dico_RECON, ...
            already_mapped_tag, consensus_proportion, epsilon, optional_settings);
A_final(A,i)= 1;
    
delete *.log
poolobj = gcp('nocreate');
delete(poolobj);
models_keep_single = A_final;

save('models_keep_single_CCLE_iHsa_biomass_medium','models_keep_single')

end

model = Cmodel;
load models_keep_single_CCLE_iHsa_biomass_medium.mat
optional_settings.func = {'biomass_reaction'};

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
        end
    clear i ind grRatio grRateKO grRateWT hasEffect delRxns fluxSolution
save('KO_iHsa_biomass_CCLE_medium','grRatio_biomass_single','grRateKO_biomass_single','grRateWT_biomass_single','genelist')
