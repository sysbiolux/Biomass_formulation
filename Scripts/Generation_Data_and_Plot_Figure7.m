%% HeLa model
%% Initialise

cd('C:\Users\maria.moscardo\Desktop\Internship') %working folder
addpath(genpath(pwd)) %add all files and folders in the working folder to matlab
addpath(genpath('C:\Users\maria.moscardo\Desktop\rFASTCORMICS')) %path to rFASTCORMICS

altcolor= [255 255 255;255 204 204; 255 153 153; 255 102 102; 255 51 51;...
    255 0 0; 204 0 0; 152 0 0; 102 0 0;  51 0 0; 0 0 0]/255; %some red color
set(0,'defaultTextInterpreter','none')
set(groot, 'DefaultAxesTickLabelInterpreter', 'none')
set(groot, 'defaultLegendInterpreter','none');

load('C:\Users\maria.moscardo\Desktop\Internship\dico_recon.mat','dico_RECON')
changeCobraSolver('ibm_cplex')
%% Biomass integration
load Recon204.mat
Cmodel = CmodelR204;
% Remove biomass reaction
reactionR2 = find(~cellfun(@isempty,regexpi(Cmodel.rxns, 'biomass')));
Cmodel_no_Biomass = removeRxns(Cmodel,Cmodel.rxns(reactionR2));

model = addReaction(Cmodel_no_Biomass,'biomass_reaction','reactionFormula','0.000275265 ksi_pre29[g] + 0.00180761 5dhf[c] + 0.000362808 ksi_deg16[l] + 0.000536138 ksi_deg26[l] + 0.001822393 pgp_hs[c] + 0.009570256 man[r] + 0.015376488 tsul[m] + 0.005540551 grdp[c] + 0.002515013 dpcoa[m] + 0.00102564 m2gacpail_hs[r] + 0.00765966 glu5p[m] + 0.000237341 ksi_deg2[l] + 0.000973747 l2fn2m2masn[g] + 0.00078694 g2m8masn[r] + 0.003829277 25aics[c] + 0.012403509 cbp[c] + 0.000235716 ksi_pre36[g] + 0.005368125 cmp[n] + 0.000420047 ksi_deg20[l] + 0.008697272 mepi[c] + 0.001822393 orot5p[c] + 0.00032681 ksi_deg12[l] + 0.006312948 6pgc[c] + 0.006061608 retinal_cis_13[c] + 0.004177796 44mzym[r] + 0.001117845 m6masnA[g] + 0.015518514 csn[c] + 0.000622104 ksi_deg29[l] + 0.010974487 dhor_S[c] + 0.001984475 coa[c] + 0.002287307 nad[c] + 0.002283832 nadh[c] + 0.001880916 accoa[m] + 0.009399267 crn[c] + 0.004389131 amp[c] + 0.821254091 gly[c] + 0.223700398 gln_L[c] + 0.024831694 hco3[c] + 0.490869599 ser_L[c] + 0.093009907 cys_L[c] + 0.729477177 ala_L[c] + 0.154539641 phe_L[c] + 0.015772543 so4[l] + 0.511875474 leu_L[c] + 0.340702092 pro_L[c] + 0.111233516 tyr_L[c] + 0.004946375 gthrd[c] + 0.125934209 met_L[c] + 0.212415656 asn_L[c] + 0.044544047 h2o2[x] + 0.433249561 val_L[c] + 0.291429234 ile_L[c] + 0.032615329 trp_L[c] + 0.243935128 arg_L[c] + 0.3327661 thr_L[c] + 0.100558832 his_L[c] + 0.399707699 glu_L[c] + 0.402488643 lys_L[c] + 0.00864883 citr_L[c] + 0.011377666 orn[c] + 0.006849449 acgam[l] + 0.083995614 nh4[c] + 0.010515878 akg[c] + 0.017404708 pyr[c] + 0.042736906 cl[c] + 0.009512811 datp[n] + 0.009795616 dttp[n] + 0.00762316 dctp[n] + 0.006693507 dgtp[n] + 0.100674101 utp[c] + 0.065907989 ctp[c] + 0.067239471 gtp[c] + 0.321320883 asp_L[c] + 0.057460661 pail_hs[c] + 0.028715639 pe_hs[c] + 0.006473115 pglyc_hs[c] + 0.042673104 ps_hs[c] + 0.016848685 sphmyln_hs[c] + 45.08253367 atp[c] + 38.62691683 h2o[c] -> 45 adp[c] + 45 pi[c] + 45 h[c] + 0.349980323 ppi[c]');
model = addReaction(model, 'clpn_hs_out','clpn_hs[c] -> ');

model.rev=zeros(numel(model.lb),1);
model.rev(model.lb<0 & model.ub>0)=1;

model.lb=model.lb*1000;
model.ub=model.ub*1000;
r=find(model.ub==0 & model.lb<0);
model.S(:,r)=-model.S(:,r);
tmp=model.lb(r);
model.lb(r)=-model.ub(r);
model.ub(r)=-tmp;

%% Import and clean Data
load HeLa_fpkm_values.mat
colnames = {'HeLa'}; %sample names

fpkm = table2array(HeLafpkmvalues(:,2)); %get the numbers
%% 
% Colnames

colnames = (colnames)'; %sample names

%In that case no different conditions
conditions = {[1:length(colnames)]};
condition_names = colnames;

%% 
% rownames

rownames = HeLafpkmvalues.Gene;

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

%% Model reconstruction
% Set bounds to avoid having a inconsistent core model
load RPMI_medium_77_Essentials_BOFdat.mat
load dico_recon.mat

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

for a = 1:size(model.subSystems,1)
     model.subSystems{a}=char(model.subSystems{a}); %I did this becaus otherwise it was not a cell of character arrays
end
unpenalized = model.rxns(ismember(model.subSystems,unpenalizedSystems));
not_medium_constrained = 'EX_tag_hs(e)';

optional_settings.unpenalized = unpenalized;
optional_settings.func = {'biomass_reaction'};

optional_settings.not_medium_constrained = not_medium_constrained;
optional_settings.medium = medium;

if ~exist('models_keep_generic')
    discretized = discretize_FPKM(fpkm, colnames);
    A_final = cell(1,numel(conditions));
    for i=1:numel(conditions)
        [~, A_keep] = fastcormics_RNAseq(model, discretized(:,conditions{i}), rownames, dico_RECON, ...
            already_mapped_tag, consensus_proportion, epsilon, optional_settings);
        A_final{i} = A_keep;
    end
    delete *.log
    poolobj = gcp('nocreate');
    delete(poolobj);
    models_keep_generic = zeros(numel(model.rxns),numel(conditions));
    for i=1:numel(conditions)
        models_keep_generic(A_final{i},i) = 1;
    end
end

save('models_keep_generic_BOFdat_HeLa_biomass','models_keep_generic')

%% Single gene deletion
load models_keep_generic_BOFdat_HeLa_biomass

if ~exist('grRatio_ATP_generic')
    changeCobraSolver('ibm_cplex')
     for i=1:size(models_keep_generic,2)
        ind = find(~cellfun(@isempty, regexp(model.rxns,optional_settings.func{1})));
        model_out = removeRxns(model,model.rxns(setdiff(1:numel(model.rxns),find(models_keep_generic(:,i))))); % create model based on active reactions
        model_out = changeObjective(model_out,model.rxns(ind)); % set objective function
        %model_out = removeUnusedGenes(model_out);
        model_out.lb = model_out.lb/1000;
        model_out.ub = model_out.ub/1000;
        
        [grRatio, grRateKO, grRateWT, hasEffect, delRxns, fluxSolution, genelist] = singleGeneDeletionTamara(model_out,'FBA',[],0,1);
        
        grRatio_biomass_generic(:,i)    = grRatio;
        grRateKO_biomass_generic(:,i)   = grRateKO;
        grRateWT_biomass_generic(:,i)   = grRateWT;
        hasEffect_biomass_generic(:,i)  = hasEffect;
        delRxns_biomass_generic(:,i)    = delRxns;
        fluxSolution_biomass_generic{i} = fluxSolution;
        end
    clear i ind grRatio grRateKO grRateWT hasEffect delRxns fluxSolution
end

% Obtain essential genes
    essential_genes = zeros(size(grRatio_biomass_generic,1),size(grRatio_biomass_generic,2));
    for s=1:size(grRatio_biomass_generic,2)
        threshold = 0.5;
        B = grRatio_biomass_generic(:,s) <= threshold;
        essential_genes(B,s)=1;
    end

    EG_final=cell.empty(0, numel(colnames));
    for e=1:size(essential_genes,2)
        x =(genelist(find(essential_genes(:,e)==1)));
        [~, ia,ib] = intersect(x,dico_RECON.ENTREZ); % find the intersection and keep indices
        EG_final(e)= {dico_RECON.SYMBOL(ib)};
    end 
     
    % Enrichment test
     [enrichment] = GeneEnrichments_HeLa(EG_final,condition_names);


 save('Enrichment_results_HeLa_BOFdat_biomass','enrichment','EG_final','grRatio_biomass_generic')
 save('KO_BOFdat_HeLa_biomass','grRatio_biomass_generic','grRateKO_biomass_generic','grRateWT_biomass_generic')

 %% Enrichment with generic biomass
Cmodel = CmodelR204;

load DMEM_medium.mat
load dico_recon.mat

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

for a = 1:size(Cmodel.subSystems,1)
     Cmodel.subSystems{a}=char(Cmodel.subSystems{a}); %I did this becaus otherwise it was not a cell of character arrays
end
unpenalized = Cmodel.rxns(ismember(Cmodel.subSystems,unpenalizedSystems));
not_medium_constrained = 'EX_tag_hs(e)';

optional_settings.unpenalized = unpenalized;
optional_settings.func = {'biomass_reaction'};

optional_settings.not_medium_constrained = not_medium_constrained;
optional_settings.medium = medium;

if ~exist('models_keep_generic')
    discretized = discretize_FPKM(fpkm, colnames);
    A_final = cell(1,numel(conditions));
    for i=1:numel(conditions)
        [~, A_keep] = fastcormics_RNAseq(Cmodel, discretized(:,conditions{i}), rownames, dico_RECON, ...
            already_mapped_tag, consensus_proportion, epsilon, optional_settings);
        A_final{i} = A_keep;
    end
    delete *.log
    poolobj = gcp('nocreate');
    delete(poolobj);
    models_keep_generic = zeros(numel(Cmodel.rxns),numel(conditions));
    for i=1:numel(conditions)
        models_keep_generic(A_final{i},i) = 1;
    end
end

save('models_keep_generic_Generic_HeLa_biomass','models_keep_generic')

%% Single gene deletion
load models_keep_generic_Generic_HeLa_biomass

if ~exist('grRatio_ATP_generic')
    changeCobraSolver('ibm_cplex')
     for i=1:size(models_keep_generic,2)
        ind = find(~cellfun(@isempty, regexp(Cmodel.rxns,optional_settings.func{1})));
        model_out = removeRxns(Cmodel,Cmodel.rxns(setdiff(1:numel(Cmodel.rxns),find(models_keep_generic(:,i))))); % create model based on active reactions
        model_out = changeObjective(model_out,Cmodel.rxns(ind)); % set objective function
        %model_out = removeUnusedGenes(model_out);
        model_out.lb = model_out.lb/1000;
        model_out.ub = model_out.ub/1000;
        
        [grRatio, grRateKO, grRateWT, hasEffect, delRxns, fluxSolution, genelist] = singleGeneDeletionTamara(model_out,'FBA',[],0,1);
        
        grRatio_biomass_generic(:,i)    = grRatio;
        grRateKO_biomass_generic(:,i)   = grRateKO;
        grRateWT_biomass_generic(:,i)   = grRateWT;
        hasEffect_biomass_generic(:,i)  = hasEffect;
        delRxns_biomass_generic(:,i)    = delRxns;
        fluxSolution_biomass_generic{i} = fluxSolution;
        end
    clear i ind grRatio grRateKO grRateWT hasEffect delRxns fluxSolution
end

load Enrichment_results_HeLa_generic_biomass
% Obtain essential genes
    essential_genes = zeros(size(grRatio_biomass_generic,1),size(grRatio_biomass_generic,2));
    for s=1:size(grRatio_biomass_generic,2)
        threshold = 0.5;
        B = grRatio_biomass_generic(:,s) <= threshold;
        essential_genes(B,s)=1;
    end

    EG_final=cell.empty(0, numel(colnames));
    for e=1:size(essential_genes,2)
        x =(genelist(find(essential_genes(:,e)==1)));
        [~, ia,ib] = intersect(x,dico_RECON.ENTREZ); % find the intersection and keep indices
        EG_final(e)= {dico_RECON.SYMBOL(ib)};
    end 
     
    % Enrichment test
     [enrichment] = GeneEnrichments_HeLa(EG_final,condition_names);

 save('Enrichment_results_HeLa_generic_biomass','enrichment','EG_final','grRatio_biomass_generic')
 save('KO_generic_HeLa_biomass','grRatio_biomass_generic','grRateKO_biomass_generic','grRateWT_biomass_generic')

 %% TCGA dataset
 %% Biomass integration
load Recon204.mat
Cmodel = CmodelR204;
% Remove biomass reaction
reactionR2 = find(~cellfun(@isempty,regexpi(Cmodel.rxns, 'biomass')));
Cmodel_no_Biomass = removeRxns(Cmodel,Cmodel.rxns(reactionR2));

model = addReaction(Cmodel_no_Biomass,'biomass_reaction','reactionFormula','0.000275265 ksi_pre29[g] + 0.00180761 5dhf[c] + 0.000362808 ksi_deg16[l] + 0.000536138 ksi_deg26[l] + 0.001822393 pgp_hs[c] + 0.009570256 man[r] + 0.015376488 tsul[m] + 0.005540551 grdp[c] + 0.002515013 dpcoa[m] + 0.00102564 m2gacpail_hs[r] + 0.00765966 glu5p[m] + 0.000237341 ksi_deg2[l] + 0.000973747 l2fn2m2masn[g] + 0.00078694 g2m8masn[r] + 0.003829277 25aics[c] + 0.012403509 cbp[c] + 0.000235716 ksi_pre36[g] + 0.005368125 cmp[n] + 0.000420047 ksi_deg20[l] + 0.008697272 mepi[c] + 0.001822393 orot5p[c] + 0.00032681 ksi_deg12[l] + 0.006312948 6pgc[c] + 0.006061608 retinal_cis_13[c] + 0.004177796 44mzym[r] + 0.001117845 m6masnA[g] + 0.015518514 csn[c] + 0.000622104 ksi_deg29[l] + 0.010974487 dhor_S[c] + 0.001984475 coa[c] + 0.002287307 nad[c] + 0.002283832 nadh[c] + 0.001880916 accoa[m] + 0.009399267 crn[c] + 0.004389131 amp[c] + 0.821254091 gly[c] + 0.223700398 gln_L[c] + 0.024831694 hco3[c] + 0.490869599 ser_L[c] + 0.093009907 cys_L[c] + 0.729477177 ala_L[c] + 0.154539641 phe_L[c] + 0.015772543 so4[l] + 0.511875474 leu_L[c] + 0.340702092 pro_L[c] + 0.111233516 tyr_L[c] + 0.004946375 gthrd[c] + 0.125934209 met_L[c] + 0.212415656 asn_L[c] + 0.044544047 h2o2[x] + 0.433249561 val_L[c] + 0.291429234 ile_L[c] + 0.032615329 trp_L[c] + 0.243935128 arg_L[c] + 0.3327661 thr_L[c] + 0.100558832 his_L[c] + 0.399707699 glu_L[c] + 0.402488643 lys_L[c] + 0.00864883 citr_L[c] + 0.011377666 orn[c] + 0.006849449 acgam[l] + 0.083995614 nh4[c] + 0.010515878 akg[c] + 0.017404708 pyr[c] + 0.042736906 cl[c] + 0.009512811 datp[n] + 0.009795616 dttp[n] + 0.00762316 dctp[n] + 0.006693507 dgtp[n] + 0.100674101 utp[c] + 0.065907989 ctp[c] + 0.067239471 gtp[c] + 0.321320883 asp_L[c] + 0.057460661 pail_hs[c] + 0.028715639 pe_hs[c] + 0.006473115 pglyc_hs[c] + 0.042673104 ps_hs[c] + 0.016848685 sphmyln_hs[c] + 45.08253367 atp[c] + 38.62691683 h2o[c] -> 45 adp[c] + 45 pi[c] + 45 h[c] + 0.349980323 ppi[c]');
model = addReaction(model, 'clpn_hs_out','clpn_hs[c] -> ');

model.rev=zeros(numel(model.lb),1);
model.rev(model.lb<0 & model.ub>0)=1;

model.lb=model.lb*1000;
model.ub=model.ub*1000;
r=find(model.ub==0 & model.lb<0);
model.S(:,r)=-model.S(:,r);
tmp=model.lb(r);
model.lb(r)=-model.ub(r);
model.ub(r)=-tmp;

 %% Import and clean Data
load COAD_fpkm_data_for_Maria.mat
fpkm = fpkm_COAD;

%% 
% Colnames

colnames = table2array(colnames_COAD(:,1));

%In that case we have two conditions: cancer and control samples.
Cancer = [1:483];
Healthy = [484:524];

conditions = {Cancer,Healthy};
condition_names = {'Cancer','Healthy'};

%% 
% rownames

rownames = rownames;

%remove NaN entries
rownames(isnan(sum(fpkm,2))) = [];
fpkm(isnan(sum(fpkm,2)),:) = [];

%remove 0 entries
rownames(sum(fpkm,2) == 0,:) = [];
fpkm(sum(fpkm,2) == 0,:) = [];
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

% Set bounds to avoid having a inconsistent core model
load RPMI_medium_77_Essentials_BOFdat.mat
load dico_recon.mat

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

for a = 1:size(model.subSystems,1)
     model.subSystems{a}=char(model.subSystems{a}); %I did this becaus otherwise it was not a cell of character arrays
end
unpenalized = model.rxns(ismember(model.subSystems,unpenalizedSystems));
not_medium_constrained = 'EX_tag_hs(e)';

optional_settings.unpenalized = unpenalized;
optional_settings.func = {'biomass_reaction'};

optional_settings.not_medium_constrained = not_medium_constrained;
optional_settings.medium = medium;

if ~exist('models_keep_generic')
    discretized = discretize_FPKM(fpkm, colnames);
    A_final = cell(1,numel(conditions));
    for i=1:numel(conditions)
        [~, A_keep] = fastcormics_RNAseq(model, discretized(:,conditions{i}), rownames, dico_RECON, ...
            already_mapped_tag, consensus_proportion, epsilon, optional_settings);
        A_final{i} = A_keep;
    end
    delete *.log
    poolobj = gcp('nocreate');
    delete(poolobj);
    models_keep_generic = zeros(numel(model.rxns),numel(conditions));
    for i=1:numel(conditions)
        models_keep_generic(A_final{i},i) = 1;
    end
end
save('models_keep_generic_TCGA_Hela_BOFdat','models_keep_generic')
%% Single gene deletion
if ~exist('grRatio_ATP_generic')
    changeCobraSolver('ibm_cplex')
     for i=1:size(models_keep_generic,2)
        ind = find(~cellfun(@isempty, regexp(model.rxns,optional_settings.func{1})));
        model_out = removeRxns(model,model.rxns(setdiff(1:numel(model.rxns),find(models_keep_generic(:,i))))); % create model based on active reactions
        model_out = changeObjective(model_out,model.rxns(ind)); % set objective function
        %model_out = removeUnusedGenes(model_out);
        model_out.lb = model_out.lb/1000;
        model_out.ub = model_out.ub/1000;
        
        [grRatio, grRateKO, grRateWT, hasEffect, delRxns, fluxSolution, genelist] = singleGeneDeletionTamara(model_out,'FBA',[],0,1);
        
        grRatio_biomass_generic(:,i)    = grRatio;
        grRateKO_biomass_generic(:,i)   = grRateKO;
        grRateWT_biomass_generic(:,i)   = grRateWT;
        hasEffect_biomass_generic(:,i)  = hasEffect;
        delRxns_biomass_generic(:,i)    = delRxns;
        fluxSolution_biomass_generic{i} = fluxSolution;
        end
    clear i ind grRatio grRateKO grRateWT hasEffect delRxns fluxSolution
end

load Enrichment_results_BOFdat_biomass_TCGA_HeLa
% Obtain essential genes
    essential_genes = zeros(size(grRatio_biomass_generic,1),size(grRatio_biomass_generic,2));
    for s=1:size(grRatio_biomass_generic,2)
        threshold = 0.5;
        B = grRatio_biomass_generic(:,s) <= threshold;
        essential_genes(B,s)=1;
    end

    EG_final=cell.empty(0, numel(colnames));
    for e=1:size(essential_genes,2)
        x =(genelist(find(essential_genes(:,e)==1)));
        [~, ia,ib] = intersect(x,dico_RECON.ENTREZ); % find the intersection and keep indices
        EG_final(e)= {dico_RECON.SYMBOL(ib)};
    end 
     
    % Enrichment test
     [enrichment_TCGA] = GeneEnrichments_HeLa(EG_final,{'HeLa'});


 save('Enrichment_results_BOFdat_biomass_TCGA_HeLa','enrichment_TCGA','EG_final','grRatio_biomass_generic')
 save('KO_BOFdat_biomass_TCGA_HeLa','grRatio_biomass_generic','grRateKO_biomass_generic','grRateWT_biomass_generic')
 
%% Enrichment TCGA with generic biomass
Cmodel = CmodelR204;

load RPMI_charstrippedFBS_DHT_MTA_content.mat
load dico_recon.mat

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

for a = 1:size(Cmodel.subSystems,1)
     Cmodel.subSystems{a}=char(Cmodel.subSystems{a}); %I did this becaus otherwise it was not a cell of character arrays
end
unpenalized = Cmodel.rxns(ismember(Cmodel.subSystems,unpenalizedSystems));
not_medium_constrained = 'EX_tag_hs(e)';

optional_settings.unpenalized = unpenalized;
optional_settings.func = {'biomass_reaction'};

optional_settings.not_medium_constrained = not_medium_constrained;
optional_settings.medium = medium;

if ~exist('models_keep_generic')
    discretized = discretize_FPKM(fpkm, colnames);
    A_final = cell(1,numel(conditions));
    for i=1:numel(conditions)
        [~, A_keep] = fastcormics_RNAseq(Cmodel, discretized(:,conditions{i}), rownames, dico_RECON, ...
            already_mapped_tag, consensus_proportion, epsilon, optional_settings);
        A_final{i} = A_keep;
    end
    delete *.log
    poolobj = gcp('nocreate');
    delete(poolobj);
    models_keep_generic = zeros(numel(Cmodel.rxns),numel(conditions));
    for i=1:numel(conditions)
        models_keep_generic(A_final{i},i) = 1;
    end
end

save('models_keep_generic_Generic_TCGA_HeLa_biomass','models_keep_generic')

%% Single gene deletion
load models_keep_generic_Generic_TCGA_HeLa_biomass

if ~exist('grRatio_ATP_generic')
    changeCobraSolver('ibm_cplex')
     for i=1:size(models_keep_generic,2)
        ind = find(~cellfun(@isempty, regexp(Cmodel.rxns,optional_settings.func{1})));
        model_out = removeRxns(Cmodel,Cmodel.rxns(setdiff(1:numel(Cmodel.rxns),find(models_keep_generic(:,i))))); % create model based on active reactions
        model_out = changeObjective(model_out,Cmodel.rxns(ind)); % set objective function
        %model_out = removeUnusedGenes(model_out);
        model_out.lb = model_out.lb/1000;
        model_out.ub = model_out.ub/1000;
        
        [grRatio, grRateKO, grRateWT, hasEffect, delRxns, fluxSolution, genelist] = singleGeneDeletionTamara(model_out,'FBA',[],0,1);
        
        grRatio_biomass_generic(:,i)    = grRatio;
        grRateKO_biomass_generic(:,i)   = grRateKO;
        grRateWT_biomass_generic(:,i)   = grRateWT;
        hasEffect_biomass_generic(:,i)  = hasEffect;
        delRxns_biomass_generic(:,i)    = delRxns;
        fluxSolution_biomass_generic{i} = fluxSolution;
        end
    clear i ind grRatio grRateKO grRateWT hasEffect delRxns fluxSolution
end

load Enrichment_results_HeLa_generic_biomass_TCGA
% Obtain essential genes
    essential_genes = zeros(size(grRatio_biomass_generic,1),size(grRatio_biomass_generic,2));
    for s=1:size(grRatio_biomass_generic,2)
        threshold = 0.5;
        B = grRatio_biomass_generic(:,s) <= threshold;
        essential_genes(B,s)=1;
    end

    EG_final=cell.empty(0, numel(colnames));
    for e=1:size(essential_genes,2)
        x =(genelist(find(essential_genes(:,e)==1)));
        [~, ia,ib] = intersect(x,dico_RECON.ENTREZ); % find the intersection and keep indices
        EG_final(e)= {dico_RECON.SYMBOL(ib)};
    end 
     
    % Enrichment test
     [enrichment] = GeneEnrichments_HeLa(EG_final,{'HeLa'});

 save('Enrichment_results_HeLa_generic_biomass_TCGA','enrichment','EG_final','grRatio_biomass_generic')
 save('KO_generic_HeLa_biomass_TCGA','grRatio_biomass_generic','grRateKO_biomass_generic','grRateWT_biomass_generic')
 
 %% Plot results 
 load Enrichment_results_HeLa_BOFdat_biomass
Sensitivity_BOFdat_Hela = [];
Specificity_BOFdat_Hela = [];
precision_BOFdat_Hela = [];
for i = 1:size(enrichment_BOFdat,1)
    sensitivity = cell2mat(enrichment_BOFdat.Known_EG(i))/cell2mat(enrichment_BOFdat.Metabolic_cancer_genes(i));
    specificity = cell2mat(enrichment_BOFdat.Known_non_EG(i))/cell2mat(enrichment_BOFdat.non_EG_in_model(i));
    precision = cell2mat(enrichment_BOFdat.Known_EG(i))/cell2mat(enrichment_BOFdat.Predicted_EG(i));
    Sensitivity_BOFdat_Hela = [Sensitivity_BOFdat_Hela,sensitivity];
    Specificity_BOFdat_Hela = [Specificity_BOFdat_Hela,specificity];
    precision_BOFdat_Hela = [precision_BOFdat_Hela,precision];
end
load Enrichment_results_HeLa_generic_biomass
Sensitivity_Generic_Hela = [];
Specificity_Generic_Hela = [];
precision_generic_Hela = [];
for i = 1:size(enrichment_generic,1)
    sensitivity = cell2mat(enrichment_generic.Known_EG(i))/cell2mat(enrichment_generic.Metabolic_cancer_genes(i));
    specificity = cell2mat(enrichment_generic.Known_non_EG(i))/cell2mat(enrichment_generic.non_EG_in_model(i));
    precision = cell2mat(enrichment_generic.Known_EG(i))/cell2mat(enrichment_generic.Predicted_EG(i));
    Sensitivity_Generic_Hela = [Sensitivity_Generic_Hela,sensitivity];
    Specificity_Generic_Hela = [Specificity_Generic_Hela,specificity];
    precision_generic_Hela = [precision_generic_Hela,precision];
end
load Enrichment_results_HeLa_generic_biomass_TCGA.mat
Sensitivity_Generic_TCGA = [];
Specificity_Generic_TCGA = [];
precision_generic_TCGA = [];
for i = 1:size(enrichment_generic_TCGA,1)
    sensitivity = cell2mat(enrichment_generic_TCGA.Known_EG(i))/cell2mat(enrichment_generic_TCGA.Metabolic_cancer_genes(i));
    specificity = cell2mat(enrichment_generic_TCGA.Known_non_EG(i))/cell2mat(enrichment_generic_TCGA.non_EG_in_model(i));
    precision = cell2mat(enrichment_generic_TCGA.Known_EG(i))/cell2mat(enrichment_generic_TCGA.Predicted_EG(i));
    Sensitivity_Generic_TCGA = [Sensitivity_Generic_TCGA,sensitivity];
    Specificity_Generic_TCGA = [Specificity_Generic_TCGA,specificity];
    precision_generic_TCGA = [precision_generic_TCGA,precision];
end
load Enrichment_results_BOFdat_biomass_TCGA_HeLa.mat
Sensitivity_BOFdat_TCGA = [];
Specificity_BOFdat_TCGA = [];
precision_BOFdat_TCGA = [];
for i = 1:size(enrichment_BOFdat_TCGA,1)
    sensitivity = cell2mat(enrichment_BOFdat_TCGA.Known_EG(i))/cell2mat(enrichment_BOFdat_TCGA.Metabolic_cancer_genes(i));
    specificity = cell2mat(enrichment_BOFdat_TCGA.Known_non_EG(i))/cell2mat(enrichment_BOFdat_TCGA.non_EG_in_model(i));
    precision = cell2mat(enrichment_BOFdat_TCGA.Known_EG(i))/cell2mat(enrichment_BOFdat_TCGA.Predicted_EG(i));
    Sensitivity_BOFdat_TCGA = [Sensitivity_BOFdat_TCGA,sensitivity];
    Specificity_BOFdat_TCGA = [Specificity_BOFdat_TCGA,specificity];
    precision_BOFdat_TCGA = [precision_BOFdat_TCGA,precision];
end

generic_hela = Sensitivity_Generic_Hela;
bofdat_hela = Sensitivity_BOFdat_Hela;
generic_tcga_cancer = Sensitivity_Generic_TCGA(1);
bofdat_tcga_cancer = Sensitivity_BOFdat_TCGA(1);
generic_tcga_control = Sensitivity_Generic_TCGA(2);
bofdat_tcga_control = Sensitivity_BOFdat_TCGA(2);
sensitivity_per_cellLine = [generic_hela,bofdat_hela;generic_tcga_cancer,bofdat_tcga_cancer;generic_tcga_control,bofdat_tcga_control];
names_sensitivity = [{'Generic_Biomass','BOFdat_Biomass';'Generic_Biomass','BOFdat_Biomass';'Generic_Biomass','BOFdat_Biomass'}];
data = [{'HeLa','TCGA Cancer','TCGA Control'}];

imagesc((sensitivity_per_cellLine))
for a = 1:size(names_sensitivity,1)
    for b = 1:size(names_sensitivity,2)
        text(b,a,table2array(names_sensitivity(a,b)),'FontSize',10,'HorizontalAlignment','center')
    end
end
set(gca,'Ytick',[1:3],'YTickLabel',data)
set(gca,'Xtick',[])
colorbar
colormap(altcolor)


generic_hela = Specificity_Generic_Hela;
bofdat_hela = Specificity_BOFdat_Hela;
generic_tcga_cancer = Specificity_Generic_TCGA(1);
bofdat_tcga_cancer = Specificity_BOFdat_TCGA(1);
generic_tcga_control = Specificity_Generic_TCGA(2);
bofdat_tcga_control = Specificity_Generic_TCGA(2);
specificity_per_cellLine = [generic_hela,bofdat_hela;generic_tcga_cancer,bofdat_tcga_cancer;generic_tcga_control,bofdat_tcga_control];

names_specificity = [{'Generic_Biomass','BOFdat_Biomass';'Generic_Biomass','BOFdat_Biomass';'Generic_Biomass','BOFdat_Biomass'}];
data = [{'HeLa','TCGA Cancer','TCGA Control'}];

imagesc((specificity_per_cellLine))
for a = 1:size(names_specificity,1)
    for b = 1:size(names_specificity,2)
        text(b,a,table2array(names_specificity(a,b)),'FontSize',10,'HorizontalAlignment','center')
    end
end
set(gca,'Ytick',[1:3],'YTickLabel',data)
set(gca,'Xtick',[])
colorbar
colormap(altcolor)

generic_hela = precision_generic_Hela;
bofdat_hela = precision_BOFdat_Hela;
generic_tcga_cancer = precision_generic_TCGA(1);
bofdat_tcga_cancer = precision_BOFdat_TCGA(1);
generic_tcga_control = precision_generic_TCGA(2);
bofdat_tcga_control = precision_BOFdat_TCGA(2);
precision_per_cellLine = [generic_hela,bofdat_hela;generic_tcga_cancer,bofdat_tcga_cancer;generic_tcga_control,bofdat_tcga_control];
names_sensitivity = [{'Generic_R'},{'BOFdat_Biomass'}];
data = [{'HeLa','TCGA Cancer','TCGA Control'}];

test = array2table(zeros(size(precision_per_cellLine,2),size(precision_per_cellLine,1)));
for i = 1:size(precision_per_cellLine,1)
    test(:,i) = cell2table([num2cell(precision_per_cellLine(i,:))']);
end

imagesc(table2array(test)')
set(gca,'Ytick',[1:3],'YTickLabel',data)
set(gca,'Xtick',[1:size(names_sensitivity,2)],'XTickLabel',names_sensitivity)
colorbar
colormap(altcolor)