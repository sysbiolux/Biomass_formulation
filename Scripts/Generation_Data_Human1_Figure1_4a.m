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

%% Biomass integration
load Recon204.mat
model = CmodelR204;
% We remove the original biomass reaction in Recon2 to avoid overlapping.
reactionR2 = find(~cellfun(@isempty,regexpi(model.rxns, 'biomass')));
model_no_Biomass = removeRxns(model,model.rxns(reactionR2));

%The origina human1 biomass reaction is (obtained from Supplementary Data 1): 
model = addReaction(model_no_Biomass,'biomass_human1','reactionFormula','0.0012 cofactor_pool_biomass[c] + 0.1124 rna[c] + 0.0267 dna[c] + 0.2212 lipid_pool_biomass[c] + 0.4062 glycogen[c] + 0.4835 metabolite_pool_biomass[c] + 45 atp[c] + 45 h2o[c] -> biomass[c] + 45 pi[c] + 45 h[c]');

%% DNA
model = addReaction(model,'dna','reactionFormula','0.3 datp[c] + 0.2 dctp[c] + 0.2 dgtp[c] + 0.3 dttp[c] -> dna[c] + ppi[c]');

%% RNA
model = addReaction(model,'rna','reactionFormula','0.18 atp[c] + 0.3 ctp[c] + 0.34 gtp[c] + 0.18 utp[c] -> ppi[c] + rna[c]');

%% Lipid pool
model = addReaction(model,'lipid_pool','reactionFormula','0.1155 chsterol[c] + 0.5029 pchol_hs[c] + 0.1905 pe_hs[c] + 0.0692 pail_hs[c] + 0.019 ps_hs[c] + 0.0096 pglyc_hs[c] + 0.0205 clpn_hs[c] + 0.0613 sphmyln_hs[c] + 0.0115 xolest2_hs[c] -> lipid_pool_biomass[c]');							
							
%% glycogen
model = addReaction(model,'glycogen','reactionFormula','udpg[c] -> h[c] + udp[c] + glycogen[c]');

%% Metabolite pool biomass
model = addReaction(model,'metabolite_pool','reactionFormula','0.3694 glu_L[c] + 0.0998 gln_L[c] + 0.0865 asp_L[c] + 0.052 uacgam[c] + 0.0404 ala_L[c] + 0.0388 thr_L[c] + 0.034 pyr[c] + 0.0282 ser_L[c] + 0.0271 atp[c] + 0.0234 coa[c] + 0.0215 gly[c] + 0.0179 gthrd[c] + 0.0102 phpyr[c] + 0.0102 utp[c] + 0.0102 ile_L[c] + 0.0102 leu_L[c] + 0.0094 dhap[c] + 0.0088 udpg[c] + 0.0088 fdp[c] + 0.0088 val_L[c] + 0.0081 mal_L[c] + 0.0071 pro_L[c] + 0.0054 tyr_L[c] + 0.0052 ctp[c] + 0.0049 phe_L[c] + 0.0046 akg[c] + 0.0043 dhor_S[c] + 0.0039 gtp[c] + 0.0039 g6p[c] + 0.0037 met_L[c] + 0.0034 cit[c] + 0.0033 adp[c] + 0.0029 lys_L[c] + 0.0029 nad[c] + 0.0028 fum[c] + 0.0025 pser_L[c] + 0.0024 his_L[c] + 0.0022 3pg[c] + 0.0015 arg_L[c] + 0.0014 23dpg[c] + 0.0012 asn_L[c] + 0.001 trp_L[c] + 0.0008 g3p[c] + 0.0008 udp[c] + 0.0006 udpglcur[c] + 0.0006 f6p[c] + 0.0005 cys_L[c] + 0.0004 nadh[c] + 0.0004 nadph[c] + 0.0002 amp[c] + 0.0002 dcmp[c] + 0.0002 icit[c] + 0.0002 gdp[c] + 0.0002 accoa[c] + 0.0002 nadp[c] + 0.0002 r5p[c] + 0.0001 gmp[c] + 0.0001 gthox[c] + 0.0001 damp[c] -> metabolite_pool_biomass[c]');							
model = addReaction(model,'phpyr','reactionFormula','akg[c] + phe_L[c] -> glu_L[c] + phpyr[c]');
% integrating the metabolite_pool reaction, nadph was giving inconsistency.
% Hence I decided to add more nadph to the system, since could be that the
% addition of extra reaction lead to not cofactors enough.
model = addReaction(model,'nadph','reactionFormula',' -> nadph[c]');
%% Cofactor pool
model = addReaction(model,'cofactor_pool','reactionFormula','0.0345 btn[c] + 0.0345 adocbl[c] + 0.0345 pheme[m] + 0.0345 fadh2[c] + 0.0345 crn[c] + 0.0345 ribflv[c] + 0.0345 thbpt[c] + 0.0345 thf[c] + 0.0345 q10h2[m] + 0.0345 retinol_cis_11[c] + 0.0172 hydroxy_all_trans_retinoate_18[c] + 0.0172 hydroxy_all_trans_retinoate_18[r] + 0.0172 hydroxyvitamin_A1_4[c] + 0.0172 hydroxyvitamin_A1_4[r] + 0.0172 OH_13_cis_retinal_4[c] + 0.0172 OH_13_cis_retinal_4[r] + 0.0172 OH_retinal_4[c] + 0.0172 OH_retinal_4[r] + 0.0345 oretn[c] + 0.0172 oxo_9_cis_retinal_4[c] + 0.0172 oxo_9_cis_retinal_4[r] + 0.0172 oxo_9_cis_retinoyl_beta_glucuronide_4[c] + 0.0172 oxo_9_cis_retinoyl_beta_glucuronide_4[r] + 0.0172 oxo_all_trans_retinoate_4[c] + 0.0172 oxo_all_trans_retinoate_4[r] + 0.0345 oxoretinol_4[c] + 0.0345 epoxy_13_cis_retinoate_5_8[c] + 0.0345 retinol_9_cis[c] + 0.0172 CE2206[c] + 0.0172 CE2206[m] + 0.0172 CE2204[c] + 0.0172 CE2204[m] + 0.0345 CE1925[c] + 0.0345 CE5854[c] + 0.0345 thmpp[c] + 0.0345 pydx5p[c] -> cofactor_pool_biomass[c]'); %CE2207[c] (m00611c) and [r],  N-retinylidene-N-retinylethanolamine (m02624c), lipoic acid (m02394c) were not included in the biomass reaction							

%The addition of adocbl led to the model inconsistency so I decided to
%add adocbl to the model. 
model = addReaction(model,'adocbl_in','reactionFormula',' -> adocbl[c]');

% Similarly, the addition of thbpt was leading to model inconsistency
model = addReaction(model,'dhbpt_in',' -> dhbpt[c]');

% Similarly, the addition of CE1925 was leading to model inconsistency
model = addReaction(model,'CE1925_in',' -> CE1925[c]');

% Similarly, the addition of CE5854 was leading to model inconsistency
model = addReaction(model,'CE5854_in',' -> CE5854[c]');


model = addReaction(model, '18hydroxy_all_trans_retinoate_c', 'reactionFormula', 'h[c] + nadph[c] + o2[c] + retn[c] -> hydroxy_all_trans_retinoate_18[c] + h2o[c] + nadp[c]');
model = addReaction(model, '18hydroxy_all_trans_retinoate_r', 'reactionFormula', 'h[r] + nadph[r] + o2[r] + retn[r] -> hydroxy_all_trans_retinoate_18[r] + h2o[r] + nadp[r]');
model = addReaction(model, '4_hydroxyvitamin_A1_c', 'reactionFormula', 'h[c] + nadph[c] + o2[c] + retinol[c] -> hydroxyvitamin_A1_4[c] + h2o[c] + nadp[c]');
model = addReaction(model, '4_hydroxyvitamin_A1_r', 'reactionFormula', 'h[r] + nadph[r] + o2[r] + retinol[c] -> hydroxyvitamin_A1_4[r] + h2o[r] + nadp[r]'); %There is no retinol[r], we can either use retinol[c] or add a transporter retinol[c] -> retinol[r]
model = addReaction(model, '4_OH_cis_retinal_c', 'reactionFormula', 'retinal_cis_13[c] + h[c] + nadph[c] + o2[c] -> OH_13_cis_retinal_4[c] + h2o[c] + nadp[c]');
model = addReaction(model, '4_OH_cis_retinal_r', 'reactionFormula', 'retinal_cis_13[r] + h[r] + nadph[r] + o2[r] -> OH_13_cis_retinal_4[r] + h2o[r] + nadp[r]');
model = addReaction(model,'4_OH_retinal_c', 'reactionFormula', 'retinal[c] + h[c] + nadph[c] + o2[c] -> OH_retinal_4[c] + h2o[c] + nadp[c]');
model = addReaction(model,'4_OH_retinal_r', 'reactionFormula', 'retinal[r] + h[r] + nadph[r] + o2[r] -> OH_retinal_4[r] + h2o[r] + nadp[r]');
model = addReaction(model,'OH_9_cis_retinal_4_c', 'reactionFormula', 'retinal_cis_9[c] + h[c] + nadph[c] + o2[c] -> OH_9_cis_retinal_4[c] + h2o[c] + nadp[c]');
model = addReaction(model,'oxo_9_cis_retinal_4_c','reactionFormula', 'OH_9_cis_retinal_4[c] + h[c] + nadph[c] + o2[c] -> oxo_9_cis_retinal_4[c] + 2 h2o[c] + nadp[c]');
model = addReaction(model,'OH_9_cis_retinal_4_r', 'reactionFormula', 'retinal_cis_9[r] + h[r] + nadph[r] + o2[r] -> OH_9_cis_retinal_4[r] + h2o[r] + nadp[r]');
model = addReaction(model,'oxo_9_cis_retinal_4_r','reactionFormula', 'OH_9_cis_retinal_4[r] + h[r] + nadph[r] + o2[r] -> oxo_9_cis_retinal_4[r] + 2 h2o[r] + nadp[r]');
model = addReaction(model, 'oxo_9_cis_retinoate_4_c', 'reactionFormula','CE1617[r] + nadph[c] + o2[c] -> oxo_9_cis_retinoate_4[c] + h[c] + h2o[c] + nadp[c]'); %CE1617 is only present in the [r] form 
model = addReaction(model,'oxo_9_cis_retinoyl_beta_glucuronide_c','reactionFormula', 'oxo_9_cis_retinoate_4[c] + udpglcur[c] -> oxo_9_cis_retinoyl_beta_glucuronide_4[c] + udp[c]');
model = addReaction(model, 'oxo_9_cis_retinoate_4_r', 'reactionFormula','CE1617[r] + nadph[r] + o2[r] -> oxo_9_cis_retinoate_4[r] + h[r] + h2o[r] + nadp[r]'); 
model = addReaction(model,'oxo_9_cis_retinoyl_beta_glucuronide_r','reactionFormula', 'oxo_9_cis_retinoate_4[r] + udpglcur[r] -> oxo_9_cis_retinoyl_beta_glucuronide_4[r] + udp[r]'); 
model = addReaction(model,'oxo_all_trans_retinoate_4_c', 'reactionFormula', 'hretn[c] + h[c] + nadph[c] + o2[c] -> oxo_all_trans_retinoate_4[c] + 2 h2o[c] + nadp[c]');
model = addReaction(model,'oxo_all_trans_retinoate_4_r', 'reactionFormula', 'hretn[c] + h[r] + nadph[r] + o2[r] -> oxo_all_trans_retinoate_4[r] + 2 h2o[r] + nadp[r]');
model = addReaction(model,'oxoretinol_4_c','reactionFormula', 'nadph[c] + o2[c] + retinol[c] -> oxoretinol_4[c] + h[c] + h2o[c] + nadp[c]');
model = addReaction(model, 'epoxy_13_cis_retinoate_5_8_c', '13_cis_retn[c] + h[c] + nadph[c] + o2[c] -> epoxy_13_cis_retinoate_5_8[c] + h2o[c] + nadp[c]');
model = addReaction(model,'calcitetrol_c','reactionFormula', '2425dhvitd3[c] + h[c] + nadh[c] + o2[c] -> h2o[c] + nad[c] + calcitetrol[c]');
model = addReaction(model,'CE2206_c','reactionFormula','calcitetrol[c] + h[c] + nadph[c] + o2[c] -> CE2206[c] + 2 h2o[c] + nadp[c]');
model = addReaction(model,'calcitetrol_m','reactionFormula','2425dhvitd3[m] + h[m] + nadh[m] + o2[m] -> h2o[m] + nad[m] + calcitetrol[m]');
model = addReaction(model,'CE2206_m','reactionFormula','calcitetrol[m] + h[m] + nadph[m] + o2[m] -> CE2206[m] + 2 h2o[m] + nadp[m]');						
model = addReaction(model,'CE2201_c','reactionFormula','25hvitd3[c] + h[c] + nadph[c] + o2[c] -> CE2201[c] + h2o[c] + nadp[c]');
model = addReaction(model,'CE2202_c','reactionFormula','CE2201[c] + h[c] + nadh[c] + o2[c] -> CE2202[c] + h2o[c] + nad[c]');
model = addReaction(model,'CE2203_c','reactionFormula','CE2202[c] + o2[c] -> CE2203[c] + h2o2[c]');
model = addReaction(model,'CE2204_c','reactionFormula','CE2203[c] + o2[c] -> CE2204[c] + h2o2[c]');							
model = addReaction(model,'CE2201_m','reactionFormula','25hvitd3[m] + h[m] + nadph[m] + o2[m] -> CE2201[m] + h2o[m] + nadp[m]');
model = addReaction(model,'CE2202_m','reactionFormula','CE2201[m] + h[m] + nadh[m] + o2[m] -> CE2202[m] + h2o[m] + nad[m]');
model = addReaction(model,'CE2203_m','reactionFormula','CE2202[m] + o2[m] -> CE2203[m] + h2o2[m]');
model = addReaction(model,'CE2204_m','reactionFormula','CE2203[m] + o2[m] -> CE2204[m] + h2o2[m]');					

model = addReaction(model,'biomass out','reactionFormula','biomass[c] -> ');

A = fastcc_4_rfastcormics(model, 1e-4,0);

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

discretized = discretize_FPKM(fpkmIntestine, colnames);
load dico_medium_culture.mat
A_final = zeros(numel(model.rxns),numel(colnames));
col = cell2table(colnames);
medium_CellLines = cellstr(medium_CellLines{:, :});

for a = 1:size(model.subSystems,1)
     model.subSystems{a}=char(model.subSystems{a}); %I did this becaus otherwise it was not a cell of character arrays
end

for i=1:numel(col)   % need to define medium for the cells used here
    used_medium=medium_CellLines(i,2);
    if ~cellfun(@isempty, strfind(used_medium,'RPMI'))
         load ('RPMI_medium_88_essentials_human1.mat')
    elseif ~cellfun(@isempty, strfind(used_medium,'DMEM'))
         load ('DMEM_medium_Human1.mat')
    elseif ~cellfun(@isempty, strfind(used_medium,'EMEM'))
         load ('EMEM_medium_Human1.mat')
    elseif ~cellfun(@isempty, strfind(used_medium,'F12K'))
         load ('F12k_medium_Human1.mat')
    elseif ~cellfun(@isempty, strfind(used_medium,'McCoy'))
         load ('McCoys_medium_Human1.mat')
    elseif ~cellfun(@isempty, strfind(used_medium,'Leibovitz'))
         load ('L15_medium_Human1.mat')
    elseif ~cellfun(@isempty, strfind(used_medium,'MEM'))
         load ('MEM_medium_Human1.mat')
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
unpenalized = model.rxns(ismember(model.subSystems,unpenalizedSystems));
not_medium_constrained = 'EX_tag_hs(e)';

optional_settings.unpenalized = unpenalized;
optional_settings.func = {'biomass_human1'};

optional_settings.not_medium_constrained = not_medium_constrained;
optional_settings.medium = medium;



%% Single models
[~, A] = fastcormics_RNAseq(model, discretized(:,i), rownames, dico_RECON, ...
            already_mapped_tag, consensus_proportion, epsilon, optional_settings);
A_final(A,i)= 1;
    
delete *.log
poolobj = gcp('nocreate');
delete(poolobj);
models_keep_single = A_final;

save('models_keep_single_CCLE_human1_biomass_medium','models_keep_single')

end

load models_keep_single_CCLE_human1_biomass_medium
optional_settings.func = {'biomass_human1'};

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
save('KO_human1_biomass_CCLE_medium','grRatio_biomass_single','grRateKO_biomass_single','grRateWT_biomass_single','genelist')

