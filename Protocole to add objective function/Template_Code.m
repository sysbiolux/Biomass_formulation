%% This script is a standard template for the integration of reactions coming from other GEMs into a specific GEM. 
%It is based on the paper: Moscardó García M., Pacheco M., Sauter T.
%(2021). Integration of external biomass reactions on existing metabolic
%models. Star Protocols.
% feature astheightlimit 2000 %just to be run if line 105 leads to the
% error: “Error: The input was too complicated or too big for MATLAB to
% parse”. 
%% Before you begin
%Step 2
changeCobraSolver('ibm_cplex')

%% Step 1: Extract the reactions of interest   
%Step 1
load donor_model.mat %please change according to the model name
model = model_donor; % please change according to variable name
%Step 2
NAME_BIOMASS_RXN='MAR04413'  %please change according to reaction name
Biomass_idx = find(ismember(model.rxns,NAME_BIOMASS_RXN)) %please change according to reaction name
printRxnFormula(model,model.rxns(Biomass_idx))

%% Step 2: Obtain a Recon 2 compatible biomass reaction
% This step consists of extracting the information required in the future
% steps from Metabolic Atlas (https://metabolicatlas.org)

%% Step 3: Model consistency
%Step 5
load target_model
input_model = target_model; %please change accoding to the name of the model
NAMEREACTION='biomass_reaction'  %please change according to reaction name

idx = find(~cellfun(@isempty,regexpi(input_model.rxns, NAMEREACTION))); %idx = find(~cellfun(@isempty,regexpi(input_model.rxns, 'biomass_reaction')));
model_no_Biomass = removeRxns(input_model, input_model.rxns(idx));

%Step 6
model_with_new_biomass = addReaction (model_no_Biomass, 'biomass_human1', 'reactionFormula','0.0012 cofactor_pool[c] + 5.3375 protein_pool[c] + 0.1124 rna[c] + 0.0267 dna[c] + 0.2212 lipid_pool[c] + 0.4062 glycogen[c] + 0.4835 metabolite_pool[c] + 45 atp[c] + 45 h2o[c] -> biomass[c] + 45 pi[c] + 45 h[c]');

%Step 7
A = fastcc_4_rfastcormics(model_with_new_biomass, 1e-4,0);

%Step 8
METABOLITENAME= 'rna[c]'  %please change according to metabolite name
find(ismember(model_no_Biomass.mets,METABOLITENAME))
%Examples:
%find(ismember(model_no_Biomass.mets,'rna[c]'))
%find(ismember(model_no_Biomass.mets,'cofactor_pool[c]')) 
%find(ismember(model_no_Biomass.mets,'protein_pool[c]'))
%find(ismember(model_no_Biomass.mets,'dna[c]'))
%find(ismember(model_no_Biomass.mets,'lipid_pool[c]'))
%find(ismember(model_no_Biomass.mets,'glycogen[c]'))
%find(ismember(model_no_Biomass.mets,'atp[c]'))
%find(ismember(model_no_Biomass.mets,'h2o[c]'))

%or
Biomass_mets = {'cofactor_pool[c]','protein_pool[c]','dna[c]','rna[c]','lipid_pool[c]','glycogen[c]','atp[c]','h2o[c]'};
setdiff(Biomass_mets,model_no_Biomass.mets) 

%Step 9
%example of added reactions, others might be needed
model_with_new_biomass = addReaction(model_with_new_biomass,'RNA_Production','reactionFormula','0.18 atp[c] + 0.3 ctp[c] + 0.34 gtp[c] + 0.18 utp[c] -> ppi[c] + rna[c]')

%after the addition of all the required recations you should have a
%consistent model that we stored under 'H1Biomass_R2.mat

load('H1Biomass_R2.mat') % model with added reactions

%Please make sure that after the addition or removal of some metabolites,
%the biomass reaction is standardise. If the metabolite formula is
%available you can directly apply the formula: 1+sum(n_i*MW_j ) - sum(n_j*MW_j) 
%and then divide all the coefficients of the biomass reaction by that
%number, to readjust them. Alternatively, we refer the user to
%https://github.com/maranasgroup/BiomassMW/tree/master/MatlabCobraToolbox
%where Chan et al., 2017 provided a way to obtain the metabolite formulas
%and their corresponding MW. 

%Step 10
A = fastcc_4_rfastcormics(model_with_new_biomass, 1e-4,0);

%% Step 4: Test model performance
%Step 11
load DATAOFINTEREST.mat
%DATAOFINTEREST.mat should be provided in a table form, with column names
%containing the names of the samples, row names containing the identifiers
%of the measured genes and the FPKM values. If they are provided in a table
%form, the following lines of code could be use to generate 3 variables:
% colnames = DATAOFINTEREST.Properties.VariableNames; 
%fpkm = table2array(DATAOFINTEREST);
%rownames = DATAOFINTEREST.Name;

%Discretise the data
discretized = discretize_FPKM(fpkm, colnames);

%Step 12
load DICTIONARY.mat
%The dictionary should have at least 2 columns mapping the model gene
%identifiers and the data identifiers. 

%Step 13

NAME_OBJECTIVE_FUNCTION={'biomass_human1'} % change according to the opjective function in your model
biomass_reaction=NAME_OBJECTIVE_FUNCTION;
optional_settings.func = NAME_OBJECTIVE_FUNCTION; 
epsilon = 1e-4;
consensus_proportion = 0.9; 
already_mapped_tag = 0;
[~, A] = fastcormics_RNAseq(model_with_new_biomass, discretized, rownames, DICTIONARY,biomass_reaction, already_mapped_tag, consensus_proportion, epsilon, optional_settings); %if this line lead to the error: “Error: The input was too complicated or too big for MATLAB to parse”. Please run line 5

%Step 14
ind = find(ismember(model_with_new_biomass.rxns, optional_settings.func{1}));
model_out = removeRxns(model_with_new_biomass, model_with_new_biomass.rxns(setdiff(1:numel(model_with_new_biomass.rxns),A))); 
model_out = changeObjective(model_out, model_with_new_biomass.rxns(ind));
[grRatio, grRateKO, grRateWT, hasEffect, delRxns, fluxSolution,genelist] = singleGeneDeletion(model_out,'FBA',[],0,1);
grRatio_biomass = grRatio;
grRateKO_biomass = grRateKO;
grRateWT_biomass = grRateWT; 

%Step 15
essential_genes=zeros(size(grRatio_biomass,1),size(grRatio_biomass,2));
threshold = 0.5;

B = grRatio_biomass(:,1)<=threshold;
essential_genes(B,1) = 1;
Predicted_Essential(1).geneList=genelist(find(essential_genes));


load('EssentialGenesBinary.mat')
load('DICTIONARY_ESSENTIAL_GENES.mat')
confusion_table = computeEnrichment_tests(Predicted_Essential,colnames,data_essential,colnames_essential, rownames_essential, dico_essential,genelist);
confusion_table.Specificity
confusion_table.Precision
confusion_table.Sensitivity

%Step 16
ind = find(~cellfun(@isempty,regexp(model_with_new_biomass.rxns, optional_settings.func{1})));
model_out = removeRxns(model_with_new_biomass, model_with_new_biomass.rxns(setdiff(1:numel(model_with_new_biomass.rxns),A))); 
model_out = changeObjective(model_out, model_with_new_biomass.rxns(ind));
sol = optimizeCbModel(model_out);
optimization_results = sol.f;  

