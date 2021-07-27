%This script is equivalent to GrowthRate_CRC_CCLE. However, in this case we
%will include all the reactions present in the Zielinki et al., 2017 paper,
%to be able to constrain all of them. 
%% Initialise
cd('C:\Users\maria.moscardo\Desktop\Internship') %working folder
addpath(genpath(pwd)) %add all files and folders in the working folder to matlab
addpath(genpath('C:\Users\maria.moscardo\Desktop\rFASTCORMICS')) %path to rFASTCORMICS

altcolor= [255 255 255;255 204 204; 255 153 153; 255 102 102; 255 51 51;...
    255 0 0; 204 0 0; 152 0 0; 102 0 0;  51 0 0; 0 0 0]/255; %some red color
set(0,'defaultTextInterpreter','none')
set(groot, 'DefaultAxesTickLabelInterpreter', 'none')
set(groot, 'defaultLegendInterpreter','none');

load('dico_recon.mat','dico_RECON')

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
NCI60 = {'HCT116_LARGE_INTESTINE';'HCT15_LARGE_INTESTINE';'HT29_LARGE_INTESTINE';'SW620_LARGE_INTESTINE';'KM12_LARGE_INTESTINE'};

idx_NCI = [];
for i=1:size(NCI60)
    idx = find(ismember(colnames,NCI60(i)));
    idx_NCI = [idx_NCI,idx];
end

CRC = [INTESTINE(:,[1:2]),INTESTINE(:,idx_NCI)];
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

%% Create models Standard Recon2
%% Initialise
load('Recon204.mat','CmodelR204');
Cmodel = CmodelR204; 
Cmodel=addReaction(Cmodel, 'DM_gudac_c_', 'gudac[c] <=> '); % Add this reaction, present in the experimental data
Cmodel.subSystems(end)=Cmodel.subSystems{end};
Cmodel.rev=zeros(numel(Cmodel.rxns,1))
Cmodel.rev(Cmodel.lb<0 & Cmodel.ub>0)=1;
sens=find(Cmodel.lb<0 & Cmodel.ub==0); %Find uptake reactions
Cmodel.S(:,sens)=-Cmodel.S(:,sens);
Cmodel.ub(sens)= -Cmodel.lb(sens);
Cmodel.lb(sens)=0; 
 
load dico_recon.mat
A_final = zeros(numel(Cmodel.rxns),numel(colnames));
col = cell2table(colnames);
discretized = discretize_FPKM(fpkmIntestine, colnames);

load RPMI_charstrippedFBS_DHT_MTA_content.mat
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
optional_settings.func = {'biomass_reaction';'EX_glc(e)';'EX_lac_L(e)';'EX_gln_L(e)';'EX_glu_L(e)';'EX_asp_L(e)';'EX_asn_L(e)';'EX_pro_L(e)';'EX_arg_L(e)';'EX_ala_L(e)';'EX_ser_L(e)';'EX_gly(e)';'EX_lys_L(e)';'EX_trp_L(e)';'EX_leu_L(e)';'EX_tyr_L(e)';'EX_phe_L(e)';'EX_ile_L(e)';'EX_val__L_e';'EX_thr_L(e)';'EX_orn(e)';'EX_cit(e)';'EX_mal_L(e)';'DM_atp_c_';'EX_o2(e)'}

for i=1:size(colnames,1) 
%% Single models
[~, A] = fastcormics_RNAseq(Cmodel, discretized(:,i), rownames, dico_RECON, ...
            already_mapped_tag, consensus_proportion, epsilon, optional_settings);
 A_final(A,i)= 1;
    
delete *.log
poolobj = gcp('nocreate');
delete(poolobj);

end
models_keep_single = A_final;


% save('models_CRC_CCLE_Growth_rate_MP','models_keep_single')
%% Set the new constraints according to Zielinski et al 2017 except biomass bounds
% 
upper_bounds_data = readtable('UpperBounds_Zielinski_withoutBiomass.xlsx');
upper_bounds = table2array(upper_bounds_data(:,2:end));
reactions = table2array(upper_bounds_data(:,1));
lower_bounds_data = readtable('LowerBounds_Zielinski_withoutBiomass.xlsx');
lower_bounds = table2array(lower_bounds_data(:,2:end));
Cell_lines = erase(colnames,'_LARGE_INTESTINE');
changeCobraSolver('ibm_cplex')

%% Constrain EX reactions to the 1% of the sum of the experimental values.
optional_settings.func = {'biomass_reaction'};

production_rxns = cell(1,size(models_keep_single,2));
modified_bounds = cell(1,size(models_keep_single,2));
problem = cell(1,size(models_keep_single,2));
for i = 1:size(models_keep_single,2)
    cell_line = find(contains(upper_bounds_data.Properties.VariableNames,Cell_lines(i)));
    ind = find(~cellfun(@isempty, regexp(Cmodel.rxns,optional_settings.func{1})));
    model_out = removeRxns(Cmodel,Cmodel.rxns(setdiff(1:numel(Cmodel.rxns),find(models_keep_single(:,i))))); % create model based on active reactions
    model_out = changeObjective(model_out,Cmodel.rxns(ind));
   
    % Reset the model bounds to avoid problems identifying exchange
    % reactions
    exRxnsInd=find(sum(abs(model_out.S),1)==1); %find exchange rxns
    flip_EX_rxns = exRxnsInd(sum(model_out.S(:,exRxnsInd),1)==1); %exchange rxns written as: -> A
    model_out.S(:,flip_EX_rxns) = -model_out.S(:,flip_EX_rxns);
    tmp = model_out.lb(flip_EX_rxns);
    model_out.lb(flip_EX_rxns) = -model_out.ub(flip_EX_rxns);
    model_out.ub(flip_EX_rxns) = tmp;  
   
    for a = 1:size(reactions,1)
        if ismember(reactions(a,1),model_out.rxns)
            if (lower_bounds(a,cell_line-1)&&upper_bounds(a,cell_line-1)>0) || (upper_bounds(a,cell_line-1)>0&&lower_bounds(a,cell_line-1)==0)
                models_constraint = changeRxnBounds(model_out,reactions(a,1),0,'l');
                models_constraint = changeRxnBounds(models_constraint,reactions(a,1),upper_bounds(a,cell_line-1),'u');
                production_rxns{i} = [production_rxns{i},reactions(a,1)];
            else
                models_constraint = changeRxnBounds(model_out,reactions(a,1),0,'u');
                models_constraint = changeRxnBounds(models_constraint,reactions(a,1),lower_bounds(a,cell_line-1),'l');
            end
        
        model_out = models_constraint;
        modified_bounds{i} = [modified_bounds{i},reactions(a,1)];
        test = optimizeCbModel(models_constraint);
            if test.f==0 & test.stat~=0
                met = find(ismember(model_out.rxns,reactions(a,1)));
                metC = find(ismember(Cmodel.rxns,reactions(a,1)));
                model_out.ub(met) = Cmodel.ub(metC);
                model_out.lb(met) = Cmodel.lb(metC);
                problem{i} = [problem{i},a];
            end
        end
    end
    
    %Find other active EX reactions
    [up_exRxns_Recon, up_ex_mets_carbon_Recon,uptake_exRxnsInd_Recon]=find_uptakeCarbon_EX_Rxns(model_out,optional_settings.func);
    
    %Sum the uptake reactions
    uptake_rxns = find(contains(upper_bounds_data.Var1,setdiff(upper_bounds_data.Var1,production_rxns{i})));
    uptake_rxns(find(ismember(uptake_rxns,(find(ismember(reactions,'EX_o2(e)'))))))=[]; %Remove O2 reaction, since it is not a carbon source
    value = sum(lower_bounds(uptake_rxns,cell_line-1));
    
        % 1% of the sum will be used as lb for the sum of other uptake exchange reactions,
    % with no experimental data. To do that, a metabolite Carbon Add is
    % added to the model being the only source for the EX reactions. Then a
    % reaction producing it is added, whose lb is set to 0 and th ub is set
    % to 1% of the exp value. 
    test = addMetabolite(model_out,'CarbonAdd');
    test = addReaction(test,'EX_CarbonAdd','reactionFormula',' -> CarbonAdd');
    Ex_rxns_test = find(ismember(test.rxns,'EX_CarbonAdd'));
    test.lb(Ex_rxns_test) = 0;
    test.ub(Ex_rxns_test) = -0.01*value;
    
    reactions_ID=find(ismember(test.rxns, reactions));
    uptake_exRxnsInd_Recon=setdiff(uptake_exRxnsInd_Recon, reactions_ID);
  
    for z = 1:numel(uptake_exRxnsInd_Recon)
        if test.rev(uptake_exRxnsInd_Recon(z))==1
            %Uptake reaction
            test.rev(uptake_exRxnsInd_Recon(z))=0;
            test.lb(uptake_exRxnsInd_Recon(z)) = -1000;
            test.ub(uptake_exRxnsInd_Recon(z)) = 0;
            %Secretion reaction
            as = find(test.S(:,uptake_exRxnsInd_Recon(z))== -1);
            test.S(:,end+1)=0;
            test.S(as,end)=-1;
            test.rxns(end+1)=strcat('EX_sec_',test.mets(as));
            test.lb(end+1)=0;
            test.ub(end+1)=1000;
            test.rev(end+1)=0;
            test.c(end+1)=0;
        end    
    end
    
    CarbonAdd_met = find(ismember(test.mets,'CarbonAdd'));
    test.S(CarbonAdd_met,uptake_exRxnsInd_Recon) = 1;
    
    modified_bounds_per_model{i} = modified_bounds;
    problem_Recon2 {i} = problem;
    sol = optimizeCbModel(test);
    optimization_results_Recon2{i} = sol.f;
    
end

save('Growth_rate_Recon2_Generic_FinalScript','optimization_results_Recon2','problem')
%% Calculate growht rate based on doubling time
% growht rate = ln(2)/doubling time
annotation = readtable('Doubling_times_OConnor_1997.xlsx');

growth_rate = [];
for i = 1:size(Cell_lines)
    cell_line = find(contains(annotation.CellLine,Cell_lines(i)));
    growth = log(2)/annotation.DoublingTime(cell_line);
    growth_rate = [growth_rate,growth];
end

predicted_rates = [];
for a = 1:size(optimization_results_Recon2,2)
    predicted_rates = [predicted_rates,optimization_results_Recon2{a}];
end

figure
sz = 25;
s=scatter(categorical(colnames),growth_rate,sz,'filled')
hold on
s= scatter(categorical(colnames),predicted_rates,sz,'filled')
ylim([0.00 0.06])
legend('Measured data','FBA prediction')

%% Analyzing other active exchange reactions non-constraint
for z = 1:size(flux_ex_rxns_recon2,2)
    figure 
    x = categorical(flux_ex_rxns_recon2{1,z}(:,1));
    y = cell2mat(flux_ex_rxns_recon2{1,z}(:,2));
    bar(x,y)
end

%% Create models Recon2 with Renal biomass reaction (From Recon3D)
%% Introduce the biomass reaction
load('C:\Users\maria.moscardo\Desktop\Internship\Recon204.mat','CmodelR204')
Cmodel = CmodelR204;
% We remove the original biomass reaction in Recon2 to avoid overlapping.
reactionR2 = find(~cellfun(@isempty,regexpi(Cmodel.rxns, 'biomass')));
Cmodel_no_Biomass = removeRxns(Cmodel,Cmodel.rxns(reactionR2));

%The problem is that CE1925 and CE5854 (for the vitamin E reaction) are in loops. Additionally, they are
%only produced by the reactions which are within the loop, so we cannot 
%remove only of them. Hence we remove the reactions which are in the loop 
%and we include them through exchange reactions.
new_model = removeRxns(Cmodel_no_Biomass,'RE2404C');
new_model = removeRxns(new_model,'RE3381C');
new_model = removeRxns(new_model,'RE2405C');
new_model = removeRxns(new_model,'RE2541C');

new_model = addReaction(new_model, 'CE1925_in','reactionFormula',' -> CE1925[c]');
new_model = addReaction(new_model,'CE5854_in','reactionFormula',' -> CE5854[c]');

%% M03161
% idx = find(contains(Recon3D.mets,'M03161[c]'));
% involved_rxn = find(Recon3D.S(idx,:)~=0); %The reaction producing glycogen is Recon3D.rxns(12323), the first one in involved_rxn

Recon2_renalBiomass = addReaction(new_model,'glycogen','reactionFormula','udpg[c] -> h[c] + udp[c] + M03161[c]');

%% M01602 
% idx = find(contains(Recon3D.mets,'M01602'));
% involved_rxn = find(Recon3D.S(idx,:)~=0); %Only cofactors_vitamins reaction will be considered, since the renal_biomass reaction will be added later on

%VitaminE reaction will be included according to the HMR vitaminE reaction
%(https://www.vmh.life/#reaction/vitaminE)
Recon2_renalBiomass = addReaction(Recon2_renalBiomass,'vitaminE','reactionFormula','CE1925[c] + CE5854[c] -> M03143[c]');

%VitaminD will be represented by the metabolite vitd3[c] which is the
%result of the vitaminD3 formation reaction

%VitaminA will be recreated using the formula: 11-cis-retinol[c] + 18-hydroxy-all-trans-retinoate[c] + 18-hydroxy-all-trans-retinoate[r] + 4-hydroxyvitamin A1[c] + 4-hydroxyvitamin A1[r] + 4-OH-13-cis-retinal[c] + 4-OH-13-cis-retinal[r] + 4-OH-retinal[c] + 4-OH-retinal[r] + 4-oxo-13-cis-retinoate[r] + 4-oxo-9-cis-retinal[c]+ 4-oxo-9-cis-retinal[r] + 4-oxo-9-cis-retinoyl-beta-glucuronide[c] + 4-oxo-9-cis-retinoyl-beta-glucuronide[r] + 4-oxo-all-trans-retinoate[c] + 4-oxo-all-trans-retinoate[r] + 4-oxoretinol[c] + 5,8-epoxy-13-cis-retinoate[c] + 9-cis-retinol[c] + N-retinylidene-N-retinylethanolamine[c] ? vitamin A derivatives[c]
Recon2_renalBiomass = addReaction(Recon2_renalBiomass, '18hydroxy_all_trans_retinoate_c', 'reactionFormula', 'h[c] + nadph[c] + o2[c] + retn[c] -> hydroxy_all_trans_retinoate_18[c] + h2o[c] + nadp[c]');
Recon2_renalBiomass = addReaction(Recon2_renalBiomass, '18hydroxy_all_trans_retinoate_r', 'reactionFormula', 'h[r] + nadph[r] + o2[r] + retn[r] -> hydroxy_all_trans_retinoate_18[r] + h2o[r] + nadp[r]');
Recon2_renalBiomass = addReaction(Recon2_renalBiomass, '4_hydroxyvitamin_A1_c', 'reactionFormula', 'h[c] + nadph[c] + o2[c] + retinol[c] -> hydroxyvitamin_A1_4[c] + h2o[c] + nadp[c]');
Recon2_renalBiomass = addReaction(Recon2_renalBiomass, '4_hydroxyvitamin_A1_r', 'reactionFormula', 'h[r] + nadph[r] + o2[r] + retinol[c] -> hydroxyvitamin_A1_4[r] + h2o[r] + nadp[r]'); %There is no retinol[r], we can either use retinol[c] or add a transporter retinol[c] -> retinol[r]
Recon2_renalBiomass = addReaction(Recon2_renalBiomass, '4_OH_cis_retinal_c', 'reactionFormula', 'retinal_cis_13[c] + h[c] + nadph[c] + o2[c] -> OH_13_cis_retinal_4[c] + h2o[c] + nadp[c]');
Recon2_renalBiomass = addReaction(Recon2_renalBiomass, '4_OH_cis_retinal_r', 'reactionFormula', 'retinal_cis_13[r] + h[r] + nadph[r] + o2[r] -> OH_13_cis_retinal_4[r] + h2o[r] + nadp[r]');
Recon2_renalBiomass = addReaction(Recon2_renalBiomass,'4_OH_retinal_c', 'reactionFormula', 'retinal[c] + h[c] + nadph[c] + o2[c] -> OH_retinal_4[c] + h2o[c] + nadp[c]');
Recon2_renalBiomass = addReaction(Recon2_renalBiomass,'4_OH_retinal_r', 'reactionFormula', 'retinal[r] + h[r] + nadph[r] + o2[r] -> OH_retinal_4[r] + h2o[r] + nadp[r]');
Recon2_renalBiomass = addReaction(Recon2_renalBiomass,'OH_9_cis_retinal_4_c', 'reactionFormula', 'retinal_cis_9[c] + h[c] + nadph[c] + o2[c] -> OH_9_cis_retinal_4[c] + h2o[c] + nadp[c]');
Recon2_renalBiomass = addReaction(Recon2_renalBiomass,'oxo_9_cis_retinal_4_c','reactionFormula', 'OH_9_cis_retinal_4[c] + h[c] + nadph[c] + o2[c] -> oxo_9_cis_retinal_4[c] + 2 h2o[c] + nadp[c]');
Recon2_renalBiomass = addReaction(Recon2_renalBiomass,'OH_9_cis_retinal_4_r', 'reactionFormula', 'retinal_cis_9[r] + h[r] + nadph[r] + o2[r] -> OH_9_cis_retinal_4[r] + h2o[r] + nadp[r]');
Recon2_renalBiomass = addReaction(Recon2_renalBiomass,'oxo_9_cis_retinal_4_r','reactionFormula', 'OH_9_cis_retinal_4[r] + h[r] + nadph[r] + o2[r] -> oxo_9_cis_retinal_4[r] + 2 h2o[r] + nadp[r]');
Recon2_renalBiomass = addReaction(Recon2_renalBiomass, 'oxo_9_cis_retinoate_4_c', 'reactionFormula','CE1617[r] + nadph[c] + o2[c] -> oxo_9_cis_retinoate_4[c] + h[c] + h2o[c] + nadp[c]'); %CE1617 is only present in the [r] form 
Recon2_renalBiomass = addReaction(Recon2_renalBiomass,'oxo_9_cis_retinoyl_beta_glucuronide_c','reactionFormula', 'oxo_9_cis_retinoate_4[c] + udpglcur[c] -> oxo_9_cis_retinoyl_beta_glucuronide_4[c] + udp[c]');
Recon2_renalBiomass = addReaction(Recon2_renalBiomass, 'oxo_9_cis_retinoate_4_r', 'reactionFormula','CE1617[r] + nadph[r] + o2[r] -> oxo_9_cis_retinoate_4[r] + h[r] + h2o[r] + nadp[r]'); 
Recon2_renalBiomass = addReaction(Recon2_renalBiomass,'oxo_9_cis_retinoyl_beta_glucuronide_r','reactionFormula', 'oxo_9_cis_retinoate_4[r] + udpglcur[r] -> oxo_9_cis_retinoyl_beta_glucuronide_4[r] + udp[r]'); 
Recon2_renalBiomass = addReaction(Recon2_renalBiomass,'oxo_all_trans_retinoate_4_c', 'reactionFormula', 'hretn[c] + h[c] + nadph[c] + o2[c] -> oxo_all_trans_retinoate_4[c] + 2 h2o[c] + nadp[c]');
Recon2_renalBiomass = addReaction(Recon2_renalBiomass,'oxo_all_trans_retinoate_4_r', 'reactionFormula', 'hretn[c] + h[r] + nadph[r] + o2[r] -> oxo_all_trans_retinoate_4[r] + 2 h2o[r] + nadp[r]');
Recon2_renalBiomass = addReaction(Recon2_renalBiomass,'oxoretinol_4_c','reactionFormula', 'nadph[c] + o2[c] + retinol[c] -> oxoretinol_4[c] + h[c] + h2o[c] + nadp[c]');
Recon2_renalBiomass = addReaction(Recon2_renalBiomass, 'epoxy_13_cis_retinoate_5_8_c', '13_cis_retn[c] + h[c] + nadph[c] + o2[c] -> epoxy_13_cis_retinoate_5_8[c] + h2o[c] + nadp[c]');

Recon2_renalBiomass = addReaction(Recon2_renalBiomass, 'vitaminA_derivatives', 'reactionFormula', 'retinol_cis_11[c] + hydroxy_all_trans_retinoate_18[c] + hydroxy_all_trans_retinoate_18[r] + hydroxyvitamin_A1_4[c] + hydroxyvitamin_A1_4[r] + OH_13_cis_retinal_4[c] + OH_13_cis_retinal_4[r] + OH_retinal_4[c] + OH_retinal_4[r] + oretn[c] + oxo_9_cis_retinal_4[c] + oxo_9_cis_retinal_4[r] + oxo_9_cis_retinoyl_beta_glucuronide_4[c] + oxo_9_cis_retinoyl_beta_glucuronide_4[r] + oxo_all_trans_retinoate_4[c] + oxo_all_trans_retinoate_4[r] + oxoretinol_4[c] + epoxy_13_cis_retinoate_5_8[c] + retinol_9_cis[c] -> vitaminA[c]'); %All metabolites leading to vitamin A (oretn[c] is not in the r form, so we can either use the c form or a transporter oretn[c] -> oretn[r]

Recon2_renalBiomass = addReaction(Recon2_renalBiomass, 'Cofactors_Vitamins','reactionFormula', 'btn[c] + coa[c] + adocbl[c] + fadh2[c] + crn[c] + nadh[c] + nadph[c] + ribflv[c] + thbpt[c] + thf[c] + q10h2[m] + 0.1 M03143[c] + 0.1 vitd3[c] + 0.1 vitaminA[c] -> M01602[c]')

% Recon2_renalBiomass = addReaction(Recon2_renalBiomass,'all_out','reactionFormula','M01602[c] + M03161[c] -> ');

%% After adding the Cofactors_vitamins reaction the model is not consistent anymore, so we need to find the consistency. 
% since nadph[c] is the last i added and before it was consistent I thought
% the problem should be on nadph. 
Recon2_renalBiomass = addReaction(Recon2_renalBiomass,'nadph_in','reactionFormula',' -> nadph[c]'); %because when i added the metabolites independently nadph was giving problems 

%then i added adocbl[c] so i tried to solve it 
Recon2_renalBiomass = removeRxns(Recon2_renalBiomass,'CBLTDe');
Recon2_renalBiomass = removeRxns(Recon2_renalBiomass,'CBLtle');
Recon2_renalBiomass = addReaction(Recon2_renalBiomass,'adocbl_in','reactionFormula',' -> adocbl[c]');
%with this only EHGLATm is inconsistent, so i will try to remove the
%reaction of e4hglu[m] being in a loop
Recon2_renalBiomass = removeRxns(Recon2_renalBiomass,'r0686');

%after adding thbpt[c] to the reaction again the biomass reaction and
%EHGLATm are consistent. So we will try to fix biomass reaction by removing
%the two reactions in a loop
Recon2_renalBiomass = addReaction(Recon2_renalBiomass,'dhbpt_in',' -> dhbpt[c]'); %I got consistent the biomass reaction
%now PHCDm and EHGLATm are inconsistent.
Recon2_renalBiomass = addReaction(Recon2_renalBiomass,'1p3h5c_in','reactionFormula',' -> 1p3h5c[m]');

%% M02958
% idx = find(contains(Recon3D.mets,'M02958'));
% involved_rxn = find(Recon3D.S(idx,:)~=0);

%acoa present in the TAG reaction was not in Recon 2. Hence we include the
%acoa reaction to be able to further integrate the TAG reaction. 
Recon2_renalBiomass = addReaction(Recon2_renalBiomass, 'aCoA_pool', 'reactionFormula', 'crm_hs[c] + coa[c] + h[c] -> sphings[c] + acoa[c]');
Recon2_renalBiomass = addReaction(Recon2_renalBiomass, 'TAG-pool','reactionFormula', 'dag_hs[c] + acoa[c] -> coa[c] + M02958[c]');
% Recon2_renalBiomass = addReaction(Recon2_renalBiomass,'all_out','reactionFormula','M01602[c] + M03161[c] + M02958[c] -> ');

% RNA components coefficients were calculated considering the RNA coefficient in the
% renal biomass reaction and the ratios found on
% https://www.metabolicatlas.org/explore/gem-browser/human1/metabolite/m02847c
% Similarly, DNA components coefficients were calculated considering the RNA coefficient in the
% renal biomass reaction and the ratios found on
% https://www.metabolicatlas.org/explore/gem-browser/human1/metabolite/m01721c 

Recon2_renalBiomass = addReaction(Recon2_renalBiomass,'transporter_xolest2','reactionFormula','xolest2_hs[c] <=> xolest2_hs[r]'); % In Recon2 there is xolest2 e and c but no r. Hence we just add a transporter. 
Recon2_renalBiomass = addReaction(Recon2_renalBiomass,'biomass_reaction','reactionFormula','0.0022743 datp[n] + 0.0015162 dctp[n] + 0.0015162 dgtp[n] + 0.0022743 dttp[n] + 0.00554661 atp[c] + 0.00924435 ctp[c] + 0.01047693 gtp[c] + 0.00554661 utp[c] + 0.16786 glu_L[c] + 0.69277 asp_L[c] + 0.11871 ala_L[c] + 0.012056 asn_L[c] + 0.0012056 cys_L[c] + 0.19104 gln_L[c] + 0.13911 gly[c] + 0.037096 ser_L[c] + 0.020403 thr_L[c] + 0.0009274 arg_L[c] + 0.009274 lys_L[c] + 0.00057356 clpn_hs[c] + 0.0018548 met_L[c] + 0.0030972 pail_hs[c] + 0.017171 pchol_hs[c] + 0.020051 pe_hs[c] + 0.0082708 chsterol[c] + 0.0055644 tyr_L[c] + 0.028749 his_L[c] + 0.0037096 ile_L[c] + 0.011129 leu_L[c] + 0.0009274 trp_L[c] + 0.0037096 phe_L[c] + 0.0024143 sphmyln_hs[c] + 0.012056 val_L[c] + 0.017992 pro_L[c] + 0.005891 ps_hs[c] + 0.0031466 dag_hs[c] + 0.0038601 mag_hs[c] + 0.0012642 lpchol_hs[c] + 0.016414 pa_hs[c] + 0.0065713 xolest2_hs[r] + 0.027482 M02958[c] + 0.001 M01602[c] + 0.26551 M03161[c] -> '); % Nothing after the arrow because there is no biomass metabolite on Recon2, so just leave the flux open to the outside													

% Two reactions pop up as inconsistent: r0931 and r0932
Recon2_renalBiomass = removeRxns(Recon2_renalBiomass,'r0932');
Recon2_renalBiomass = addReaction(Recon2_renalBiomass,'arachd_l_in','reactionFormula',' -> arachd[l]');
% After solvin these two, CARPEPT1tc and EX_carn(e) pop up as inconsistent
Recon2_renalBiomass = removeRxns(Recon2_renalBiomass,'EX_carn(e)');
Recon2_renalBiomass = addReaction(Recon2_renalBiomass,'carn_e_in','reactionFormula',' -> carn[e]');


Recon2_renalBiomass=addReaction(Recon2_renalBiomass, 'DM_gudac_c_', 'gudac[c] <=> '); % Add this reaction, present in the experimental data
Recon2_renalBiomass.subSystems(end)=Recon2_renalBiomass.subSystems{end};
% To set the reversible reactions
Recon2_renalBiomass.rev=zeros(numel(Recon2_renalBiomass.rxns,1))
Recon2_renalBiomass.rev(Recon2_renalBiomass.lb<0 & Recon2_renalBiomass.ub>0)=1;
sens=find(Recon2_renalBiomass.lb<0 & Recon2_renalBiomass.ub==0); %Find uptake reactions
Recon2_renalBiomass.S(:,sens)=-Recon2_renalBiomass.S(:,sens);
Recon2_renalBiomass.ub(sens)= -Recon2_renalBiomass.lb(sens);
Recon2_renalBiomass.lb(sens)=0;

%% Model reconstruction
% load dico_recon.mat
% load dico_CRC_cell_lines.mat %The medium specification cannot be done here, because metabolites are missing in the medium to be able to have the biomass reaction active
% 
% discretized = discretize_FPKM(fpkmIntestine, colnames);
% upper_bounds_data = readtable('UpperBounds_Zielinski_withoutBiomass.xlsx');
% add_core=upper_bounds_data.Var1; %Define the reactions present in the experimental data as core reactions, to have all of them in the model 
% Cmodel = Recon2_renalBiomass;
% 
% load RPMI_medium_33_essentials.mat% need to define medium for the cells used here
% 
% epsilon = 1e-4;
% consensus_proportion = 0.9; %gene has to be expressed in 90% of the cases in order to be included. Only relevant if you want to create one model from different samples
% already_mapped_tag = 0;
% 
% unpenalizedSystems = {'Transport, endoplasmic reticular';
%     'Transport, extracellular';
%     'Transport, golgi apparatus';
%     'Transport, mitochondrial';
%     'Transport, peroxisomal';
%     'Transport, lysosomal';
%     'Transport, nuclear'};
% 
% for a = 1:size(Cmodel.subSystems,1)
%      Cmodel.subSystems{a}=char(Cmodel.subSystems{a}); %I did this because otherwise it was not a cell of character arrays
% end
% 
% unpenalized = Cmodel.rxns(ismember(Cmodel.subSystems,unpenalizedSystems));
% not_medium_constrained = 'EX_tag_hs(e)';
% 
% optional_settings.unpenalized = unpenalized;
% optional_settings.func = {'biomass_reaction'};
% 
% optional_settings.not_medium_constrained = not_medium_constrained;
% medium=union(medium, 'asp_L[e]');
% medium=union(medium, 'asn_L[e]');
% optional_settings.medium = medium;
% 
% %Single models
% A_final = cell(1,numel(colnames));
%     for i=1:numel(colnames)
%         [~, A_keep] = fastcormics_RNAseq_core(Cmodel, discretized(:,i), rownames, dico_RECON, ...
%             already_mapped_tag, consensus_proportion, epsilon, optional_settings,add_core);
%         A_final{i} = A_keep;1
%     end
%     delete *.log
%     poolobj = gcp('nocreate');
%     delete(poolobj);
%     
%     models_keep_single = zeros(numel(Cmodel.rxns),numel(colnames));
%     for i=1:numel(colnames)
%         models_keep_single(A_final{i},i) = 1;
%     end
% 
% save('models_CRC_CCLE_Growth_rate_renal_biomass_MP','models_keep_single')

%% Set the new constraints according to Zielinski et al 2017 except biomass bounds
load models_CRC_CCLE_Growth_rate_renal_biomass_MP.mat
optional_settings.func = {'biomass_reaction'};

upper_bounds_data = readtable('UpperBounds_Zielinski_withoutBiomass.xlsx');
upper_bounds = table2array(upper_bounds_data(:,2:end));
reactions = table2array(upper_bounds_data(:,1));
lower_bounds_data = readtable('LowerBounds_Zielinski_withoutBiomass.xlsx');
lower_bounds = table2array(lower_bounds_data(:,2:end));

Cell_lines = erase(colnames,'_LARGE_INTESTINE');
changeCobraSolver('ibm_cplex')
%% Constrain EX reactions to the 1% of the sum of the experimental values.
load models_CRC_CCLE_Growth_rate_renal_biomass_MP.mat
optional_settings.func = {'biomass_reaction'};

production_rxns = cell(1,size(models_keep_single,2));
modified_bounds = cell(1,size(models_keep_single,2));
problem_renalB = cell(1,size(models_keep_single,2));
for i = 1:size(models_keep_single,2)
    cell_line = find(contains(upper_bounds_data.Properties.VariableNames,Cell_lines(i)));
    ind = find(~cellfun(@isempty, regexp(Recon2_renalBiomass.rxns,optional_settings.func{1})));
    model_out = removeRxns(Recon2_renalBiomass,Recon2_renalBiomass.rxns(setdiff(1:numel(Recon2_renalBiomass.rxns),find(models_keep_single(:,i))))); % create model based on active reactions
    model_out = changeObjective(model_out,Recon2_renalBiomass.rxns(ind));
        
    % Reset the model bounds to avoid problems identifying exchange
    % reactions
    exRxnsInd=find(sum(abs(model_out.S),1)==1); %find exchange rxns
    flip_EX_rxns = exRxnsInd(sum(model_out.S(:,exRxnsInd),1)==1); %exchange rxns written as: -> A
    model_out.S(:,flip_EX_rxns) = -model_out.S(:,flip_EX_rxns);
    tmp = model_out.lb(flip_EX_rxns);
    model_out.lb(flip_EX_rxns) = -model_out.ub(flip_EX_rxns);
    model_out.ub(flip_EX_rxns) = tmp;   

    for a = 1:size(reactions,1)
        if ismember(reactions(a,1),model_out.rxns)
            if (lower_bounds(a,cell_line-1)&&upper_bounds(a,cell_line-1)>0) || (upper_bounds(a,cell_line-1)>0&&lower_bounds(a,cell_line-1)==0)
                models_constraint = changeRxnBounds(model_out,reactions(a,1),0,'l');
                models_constraint = changeRxnBounds(models_constraint,reactions(a,1),upper_bounds(a,cell_line-1),'u');
                production_rxns{i} = [production_rxns{i},reactions(a,1)];
            else
                models_constraint = changeRxnBounds(model_out,reactions(a,1),0,'u');
                models_constraint = changeRxnBounds(models_constraint,reactions(a,1),lower_bounds(a,cell_line-1),'l');
            end
        model_out = models_constraint;
        modified_bounds{i} = [modified_bounds{i},reactions(a,1)];
        test = optimizeCbModel(models_constraint);
            if test.f==0 & test.stat~=0
                met = find(ismember(model_out.rxns,reactions(a,1)));
                metC = find(ismember(Cmodel.rxns,reactions(a,1)));
                model_out.ub(met) = Cmodel.ub(metC);
                model_out.lb(met) = Cmodel.lb(metC);
                problem_renalB{i} = [problem_renalB{i},reactions(a,1)];
            end
        end
    end
    
    %Find other active EX reactions
    [up_exRxns_RenalB, up_ex_mets_carbon_RenalB,uptake_exRxnsInd_RenalB]=find_uptakeCarbon_EX_Rxns(model_out,optional_settings.func);
    
    %Sum the uptake reactions
    uptake_rxns = find(contains(upper_bounds_data.Var1,setdiff(upper_bounds_data.Var1,production_rxns{i})));
    uptake_rxns(find(ismember(uptake_rxns,(find(ismember(reactions,'EX_o2(e)'))))))=[]; %Remove O2 reaction, since it is not a carbon source
    value = sum(lower_bounds(uptake_rxns,cell_line-1));

    % 1% of the sum will be used as lb for the sum of other uptake exchange reactions,
    % with no experimental data. To do that, a metabolite Carbon Add is
    % added to the model being the only source for the EX reactions. Then a
    % reaction producing it is added, whose lb is set to 0 and th ub is set
    % to 1% of the exp value. 
    test = addMetabolite(model_out,'CarbonAdd');
    test = addReaction(test,'EX_CarbonAdd','reactionFormula',' -> CarbonAdd');
    Ex_rxns_test = find(ismember(test.rxns,'EX_CarbonAdd'));
    test.lb(Ex_rxns_test) = 0;
    test.ub(Ex_rxns_test) = -0.01*value;
    
    reactions_ID=find(ismember(test.rxns, reactions));
    uptake_exRxnsInd_RenalB=setdiff(uptake_exRxnsInd_RenalB, reactions_ID);
    
    for z = 1:numel(uptake_exRxnsInd_RenalB)
        if test.rev(uptake_exRxnsInd_RenalB(z))==1
            %Uptake reaction
            test.rev(uptake_exRxnsInd_RenalB(z))=0;
            test.lb(uptake_exRxnsInd_RenalB(z)) = -1000;
            test.ub(uptake_exRxnsInd_RenalB(z)) = 0;
            %Secretion reaction
            as = find(test.S(:,uptake_exRxnsInd_RenalB(z))== -1);
            test.S(:,end+1)=0;
            test.S(as,end)=-1;
            test.rxns(end+1)=strcat('EX_sec_',test.mets(as));
            test.lb(end+1)=0;
            test.ub(end+1)=1000;
            test.rev(end+1)=0;
            test.c(end+1)=0;
        end    
    end
    CarbonAdd_met = find(ismember(test.mets,'CarbonAdd'));
    test.S(CarbonAdd_met,uptake_exRxnsInd_RenalB) = 1;
    
    modified_bounds_per_model{i} = modified_bounds;
    sol = optimizeCbModel(test);
    optimization_results_Recon2_renalB{i} = sol.f;
    
     %Check why it is not working
    exRxnsInd=find(sum(abs(test.S),1)==1);
    ex_reactions = [];
    opt_values = [];
    for c = 1:size(exRxnsInd,2)
        if sol.x(exRxnsInd(c))~=0
            ex_reactions = [ex_reactions;test.rxns(exRxnsInd(c))];
            opt_values = [opt_values;num2cell(sol.x(exRxnsInd(c)))];
        end
    end
    table_exrxns = [ex_reactions,opt_values];
    flux_ex_rxns_renalB{i} = table_exrxns;
end
save('Growth_rate_Recon2_RenalB_FinalScript','optimization_results_Recon2_renalB','problem_renalB')

%% Calculate growht rate based on doubling time
% growht rate = ln(2)/doubling time
annotation = readtable('Doubling_times_OConnor_1997.xlsx');

growth_rate = [];
for i = 1:size(Cell_lines)
    cell_line = find(contains(annotation.CellLine,Cell_lines(i)));
    growth = log(2)/annotation.DoublingTime(cell_line);
    growth_rate = [growth_rate,growth];
end

predicted_rates = [];
for a = 1:size(optimization_results_Recon2_renalB,2)
    predicted_rates = [predicted_rates,optimization_results_Recon2_renalB{a}];
end

figure
sz = 25;
s=scatter(categorical(colnames),growth_rate,sz,'filled')
hold on
s= scatter(categorical(colnames),predicted_rates,sz,'filled')
%ylim([0.00 0.06])
legend('Measured data','FBA prediction')

%% Create models Recon2 with Human1 biomass reaction
%% Introduce the biomass reaction
load('C:\Users\maria.moscardo\Desktop\Internship\Recon204.mat','CmodelR204')
Cmodel = CmodelR204; 
% We remove the original biomass reaction in Recon2 to avoid overlapping.
reactionR2 = find(~cellfun(@isempty,regexpi(Cmodel.rxns, 'biomass')));
Cmodel_no_Biomass = removeRxns(Cmodel,Cmodel.rxns(reactionR2));

%The origina human1 biomass reaction is (obtained from Supplementary Data 1): 
model = addReaction(Cmodel_no_Biomass,'biomass_reaction','reactionFormula','0.0012 cofactor_pool_biomass[c] + 0.1124 rna[c] + 0.0267 dna[c] + 0.2212 lipid_pool_biomass[c] + 0.4062 glycogen[c] + 0.4835 metabolite_pool_biomass[c] + 45 atp[c] + 45 h2o[c] -> 45 pi[c] + 45 h[c]');

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
model = addReaction(model,'cofactor_pool','reactionFormula','0.0345 btn[c] + 0.0345 adocbl[c] + 0.0345 pheme[m] + 0.0345 fadh2[c] + 0.0345 crn[c] + 0.0345 ribflv[c] + 0.0345 thbpt[c] + 0.0345 thf[c] + 0.0345 q10h2[m] + 0.0345 retinol_cis_11[c] + 0.0172 hydroxy_all_trans_retinoate_18[c] + 0.0172 hydroxy_all_trans_retinoate_18[r] + 0.0172 hydroxyvitamin_A1_4[c] + 0.0172 hydroxyvitamin_A1_4[r] + 0.0172 OH_13_cis_retinal_4[c] + 0.0172 OH_13_cis_retinal_4[r] + 0.0172 OH_retinal_4[c] + 0.0172 OH_retinal_4[r] + 0.0345 oretn[c] + 0.0172 oxo_9_cis_retinal_4[c] + 0.0172 oxo_9_cis_retinal_4[r] + 0.0172 oxo_9_cis_retinoyl_beta_glucuronide_4[c] + 0.0172 oxo_9_cis_retinoyl_beta_glucuronide_4[r] + 0.0172 oxo_all_trans_retinoate_4[c] + 0.0172 oxo_all_trans_retinoate_4[r] + 0.0345 oxoretinol_4[c] + 0.0345 epoxy_13_cis_retinoate_5_8[c] + 0.0345 retinol_9_cis[c] + 0.0172 CE2206[c] + 0.0172 CE2206[m] + 0.0172 CE2204[c] + 0.0172 CE2204[m] + 0.0345 CE1925[c] + 0.0345 CE5854[c] + 0.0345 thmpp[c] + 0.0345 pydx5p[c] -> cofactor_pool_biomass[c]');							

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

%model = addReaction(model,'biomass out','reactionFormula','biomass[c] -> ');

% A = fastcc_4_rfastcormics(model, 1e-4,0);
model=addReaction(model, 'DM_gudac_c_', 'gudac[c] <=> '); % Add this reaction, present in the experimental data
model.subSystems(end)=model.subSystems{end};
% Set reversible reactions
model.rev=zeros(numel(model.rxns,1))
model.rev(model.lb<0 & model.ub>0)=1;
sens=find(model.lb<0 & model.ub==0); %Find uptake reactions
model.S(:,sens)=-model.S(:,sens);
model.ub(sens)= -model.lb(sens);
model.lb(sens)=0;
%% Model reconstruction
% load dico_recon.mat
% load dico_CRC_cell_lines.mat %The medium specification cannot be done here, because metabolites are missing in the medium to be able to have the biomass reaction active
% 
% discretized = discretize_FPKM(fpkmIntestine, colnames);
% upper_bounds_data = readtable('UpperBounds_Zielinski_withoutBiomass.xlsx');
% add_core=upper_bounds_data.Var1;
% 
% Cmodel = model;
% 
% load RPMI_medium_88_essentials_human1.mat % need to define medium for the cells used here
% 
% epsilon = 1e-4;
% consensus_proportion = 0.9; %gene has to be expressed in 90% of the cases in order to be included. Only relevant if you want to create one model from different samples
% already_mapped_tag = 0;
% 
% unpenalizedSystems = {'Transport, endoplasmic reticular';
%     'Transport, extracellular';
%     'Transport, golgi apparatus';
%     'Transport, mitochondrial';
%     'Transport, peroxisomal';
%     'Transport, lysosomal';
%     'Transport, nuclear'};
% 
% for a = 1:size(Cmodel.subSystems,1)
%      Cmodel.subSystems{a}=char(Cmodel.subSystems{a}); %I did this because otherwise it was not a cell of character arrays
% end
% 
% unpenalized = Cmodel.rxns(ismember(Cmodel.subSystems,unpenalizedSystems));
% not_medium_constrained = 'EX_tag_hs(e)';
% 
% optional_settings.unpenalized = unpenalized;
% optional_settings.func = {'biomass_reaction'};
% 
% optional_settings.not_medium_constrained = not_medium_constrained;
% medium=union(medium, 'asp_L[e]');
% medium=union(medium, 'asn_L[e]');
% optional_settings.medium = medium;
% 
% %Single models
% A_final = cell(1,numel(colnames));
%     for i=1:numel(colnames)
%         [~, A_keep] = fastcormics_RNAseq_core(Cmodel, discretized(:,i), rownames, dico_RECON, ...
%             already_mapped_tag, consensus_proportion, epsilon, optional_settings, add_core);
%         A_final{i} = A_keep;1
%     end
%     delete *.log
%     poolobj = gcp('nocreate');
%     delete(poolobj);
%     
%     models_keep_single = zeros(numel(Cmodel.rxns),numel(colnames));
%     for i=1:numel(colnames)
%         models_keep_single(A_final{i},i) = 1;
%     end
% 
% save('models_CRC_CCLE_Growth_rate_human1_biomass_MP','models_keep_single')
%% Set the new constraints according to Zielinski et al 2017 except biomass bounds
load models_CRC_CCLE_Growth_rate_human1_biomass_MP.mat
optional_settings.func = {'biomass_reaction'};

upper_bounds_data = readtable('UpperBounds_Zielinski_withoutBiomass.xlsx');
upper_bounds = table2array(upper_bounds_data(:,2:end));
reactions = table2array(upper_bounds_data(:,1));
lower_bounds_data = readtable('LowerBounds_Zielinski_withoutBiomass.xlsx');
lower_bounds = table2array(lower_bounds_data(:,2:end));

Cell_lines = erase(colnames,'_LARGE_INTESTINE');
changeCobraSolver('ibm_cplex')

%% Add constrain to 1% for other exchange reactions
production_rxns = cell(1,size(models_keep_single,2));
modified_bounds = cell(1,size(models_keep_single,2));
problem = cell(1,size(models_keep_single,2));

for i = 1:size(models_keep_single,2)
    cell_line = find(contains(upper_bounds_data.Properties.VariableNames,Cell_lines(i)));
    ind = find(~cellfun(@isempty, regexp(model.rxns,optional_settings.func{1})));
    model_out = removeRxns(model,model.rxns(setdiff(1:numel(model.rxns),find(models_keep_single(:,i))))); % create model based on active reactions
    model_out = changeObjective(model_out,model.rxns(ind));
    
    % Reset the model bounds to avoid problems identifying exchange
    % reactions
    exRxnsInd=find(sum(abs(model_out.S),1)==1); %find exchange rxns
    flip_EX_rxns = exRxnsInd(sum(model_out.S(:,exRxnsInd),1)==1); %exchange rxns written as: -> A
    model_out.S(:,flip_EX_rxns) = -model_out.S(:,flip_EX_rxns);
    tmp = model_out.lb(flip_EX_rxns);
    model_out.lb(flip_EX_rxns) = -model_out.ub(flip_EX_rxns);
    model_out.ub(flip_EX_rxns) = tmp;   
    
    %Constrain all the reactions.
    for a = 1:size(reactions,1)
        if ismember(reactions(a,1),model_out.rxns)
            if (lower_bounds(a,cell_line-1)&&upper_bounds(a,cell_line-1)>0) || (upper_bounds(a,cell_line-1)>0&&lower_bounds(a,cell_line-1)==0)
                models_constraint = changeRxnBounds(model_out,reactions(a,1),0,'l');
                models_constraint = changeRxnBounds(models_constraint,reactions(a,1),upper_bounds(a,cell_line-1),'u');
                production_rxns{i} = [production_rxns{i},reactions(a,1)];
            else
                models_constraint = changeRxnBounds(model_out,reactions(a,1),0,'u');
                models_constraint = changeRxnBounds(models_constraint,reactions(a,1),lower_bounds(a,cell_line-1),'l');
            end
        model_out = models_constraint;
        modified_bounds{i} = [modified_bounds{i},reactions(a,1)];
        
        %Find if any constraint is dropping the flux to 0 and reset the
        %bounds.
        test = optimizeCbModel(models_constraint);
            if test.f==0 & test.stat~=0
                met = find(ismember(model_out.rxns,reactions(a,1)));
                metC = find(ismember(Cmodel.rxns,reactions(a,1)));
                model_out.ub(met) = Cmodel.ub(metC);
                model_out.lb(met) = Cmodel.lb(metC);
                problem{i} = [problem{i},reactions(a,1)];
            end
        end
    end
    
   %Find other active EX reactions
    [up_exRxns_human1B, up_ex_mets_carbon_human1B,uptake_exRxnsInd_human1B]=find_uptakeCarbon_EX_Rxns(model_out,optional_settings.func);
    % The number of exchange reactions is higher in human 1 than in
    % the other mdoels
    %Sum the uptake reactions
    uptake_rxns = find(contains(upper_bounds_data.Var1,setdiff(upper_bounds_data.Var1,production_rxns{i})));
    uptake_rxns(find(ismember(uptake_rxns,(find(ismember(reactions,'EX_o2(e)'))))))=[]; %Remove O2 reaction, since it is not a carbon source
    value = sum(lower_bounds(uptake_rxns,cell_line-1));
    
    % 1% of the sum will be used as lb for the sum of other uptake exchange reactions,
    % with no experimental data. To do that, a metabolite Carbon Add is
    % added to the model being the only source for the EX reactions. Then a
    % reaction producing it is added, whose lb is set to 0 and th ub is set
    % to 1% of the exp value. 
    test = addMetabolite(model_out,'CarbonAdd');
    test = addReaction(test,'EX_CarbonAdd','reactionFormula',' -> CarbonAdd');
    Ex_rxns_test = find(ismember(test.rxns,'EX_CarbonAdd'));
    test.lb(Ex_rxns_test) = 0;
    test.ub(Ex_rxns_test) = -0.01*value;
    
    reactions_ID=find(ismember(test.rxns, reactions));
    uptake_exRxnsInd_human1B=setdiff(uptake_exRxnsInd_human1B, reactions_ID);
  
    for z = 1:numel(uptake_exRxnsInd_human1B)
        if test.rev(uptake_exRxnsInd_human1B(z))==1
            %Uptake reaction
            test.rev(uptake_exRxnsInd_human1B(z))=0;
            test.lb(uptake_exRxnsInd_human1B(z)) = -1000;
            test.ub(uptake_exRxnsInd_human1B(z)) = 0;
            %Secretion reaction
            as = find(test.S(:,uptake_exRxnsInd_human1B(z))== -1);
            test.S(:,end+1)=0;
            test.S(as,end)=-1;
            test.rxns(end+1)=strcat('EX_sec_',test.mets(as));
            test.lb(end+1)=0;
            test.ub(end+1)=1000;
            test.rev(end+1)=0;
            test.c(end+1)=0;
        end    
    end
    
    CarbonAdd_met = find(ismember(test.mets,'CarbonAdd'));
    test.S(CarbonAdd_met,uptake_exRxnsInd_human1B) = 1;
    
    modified_bounds_per_model{i} = modified_bounds;
    sol = optimizeCbModel(test);
    optimization_results_Recon2_human1B{i} = sol.f;
    solution_opt_human1{i} = sol;
    temp1=test.S(CarbonAdd_met,:);
    temp3_Human1{i}=[test.lb(find(temp1)), test.ub(find(temp1)), sol.x(find(temp1))];
    
     %Check why it is not working
    exRxnsInd=find(sum(abs(test.S),1)==1);
    ex_reactions = [];
    opt_values = [];
    non_0_flux = [];
    for c = 1:size(exRxnsInd,2)
        if sol.x(exRxnsInd(c))~=0
            ex_reactions = [ex_reactions;test.rxns(exRxnsInd(c))];
            opt_values = [opt_values;num2cell(sol.x(exRxnsInd(c)))];
        end
    end
    table_exrxns2 = [ex_reactions,opt_values];
    flux_ex_rxns_Human1B2{i} = table_exrxns2;
    flux_C_rxns{i} = sol.x(uptake_exRxnsInd_human1B);
end

save('Growth_rate_Recon2_Human1B_FinalScript','optimization_results_Recon2_human1B','problem')
%% Calculate growht rate based on doubling time
% growht rate = ln(2)/doubling time
annotation = readtable('Doubling_times_OConnor_1997.xlsx');

growth_rate = [];
for i = 1:size(Cell_lines)
    cell_line = find(contains(annotation.CellLine,Cell_lines(i)));
    growth = log(2)/annotation.DoublingTime(cell_line);
    growth_rate = [growth_rate,growth];
end

predicted_rates = [];
for a = 1:size(optimization_results_Recon2_human1B,2)
    predicted_rates = [predicted_rates,optimization_results_Recon2_human1B{a}];
end

figure
sz = 25;
s=scatter(categorical(colnames),growth_rate,sz,'filled')
hold on
s= scatter(categorical(colnames),predicted_rates,sz,'filled')
%ylim([0 0.9])
legend('Measured data','FBA prediction')

%% Create models Recon2 with iHsa biomass reaction
%% Introduce the biomass reaction
load('C:\Users\maria.moscardo\Desktop\Internship\Recon204.mat','CmodelR204')
Cmodel = CmodelR204; 
% We remove the original biomass reaction in Recon2 to avoid overlapping.
reactionR2 = find(~cellfun(@isempty,regexpi(Cmodel.rxns, 'biomass')));
Cmodel_no_Biomass = removeRxns(Cmodel,Cmodel.rxns(reactionR2));

% misc_biomass_metabolite_synthesis in Recon2
model = addReaction(Cmodel_no_Biomass,'misc_biomass_metabolite_synthesis','reactionFormula','0.001 bhb[c] + 0.001 3mob[c] + 0.001 3pg[c] + 0.002 4abut[c] + 0.004 adp[c] + 0.004 akg[c] + 0.001 amp[c] + 0.038 ascb_L[c] + 0.088 ala_B[c] + 0.001 betaD_glucose6p[c] + 0.001 glyb[c] + 0.002 camp[c] + 0.001 carn[c] + 0.001 cdp[c] + 0.002 chol[c] + 0.001 cis_aconiate[c] + 0.015 cit[c] + 0.005 citr_L[c] + 0.001 cmp[c] + 0.001 coa[c] + 0.004 creat[c] + 0.001 crtn[c] + 0.002 D_gluconic_acid[c] + 0.001 dhap[c] + 0.001 e4p[c] + 0.001 fol[c] + 0.001 fdp[c] + 0.001 f6p[c] + 0.003 fum[c] + 0.001 gdp[c] + 0.001 g1p[c] + 0.001 gmp[c] + 0.49 gthrd[c] + 0.034 gthox[c] + 0.004 hcys_L[c] + 0.001 icit[c] + 0.002 cyst_L[c] + 0.047 lac_L[c] + 0.01 mal_L[c] + 0.022 nad[c] + 0.002 nadp[c] + 0.008 orn[c] + 0.001 pep[c] + 0.002 prpp[c] + 0.001 ptrc[c] + 0.001 pydx5p[c] + 0.001 r5p[c] + 0.001 ru5p_D[c] + 0.002 ahcys[c] + 0.003 amet[c] + 0.001 s7p[c] + 0.002 glyc3p[c] + 0.001 spmd[c] + 0.001 sprm[c] + 0.002 succ[c] + 0.16 taur[c] + 0.001 thf[c] + 0.001 4hpro_LT[c] + 0.003 q10h2[m] + 0.005 q10[m] + 0.001 udp[c] + 0.001 ump[c] -> m90199[c]');

% include betaD_glucose_6phosphate in the system
model = addReaction(model,'betaD_g','reactionFormula','atp[c] + beta_D_glucose[c] -> adp[c] + betaD_glucose6p[c] + h[c]');
model = addReaction(model,'betaD_g6p','reactionFormula','beta_D_glucose[c] <=> glc_D[c]');

% include cis-aconiate in the system
model = addReaction(model,'cis_aconiate','reactionFormula','cis_aconiate[c] + h2o[c] <=> cit[c]');

% include creatinine
model = addReaction(model,'creatinine','reactionFormula',' -> crtn[c]');
 
% include gluconic acid in the system 
model = addReaction(model,'gluconic_acid','reactionFormula','6pgc[c] + adp[c] + h[c] -> atp[c] + D_gluconic_acid[c]');

% The model was not consistent after adding all the metabolites, so I added
% them one by one and I realized that nadp[c] was giving problems. TO solve
% the inconsistency: 
model = addReaction(model,'nadp','reactionFormula',' -> nadp[c]');

%% misc biomass metabolite demand
%model = addReaction(model,'misc biomass metabolite demand','reactionFormula','m90199[c] -> ');

%% biomass synthesis
model = addReaction(model,'biomass_reaction','reactionFormula','33 dna[c] + 50 rna[c] + 10 m90191[c] + 120 m90192[c] + m90193[c] + 87 m90194[c] + 42 m90195[c] + 160 m90199[c] -> ');

% dna is not in Recon2
model = addReaction(model,'dna','reactionFormula','0.3 datp[c] + 0.2 dctp[c] + 0.2 dgtp[c] + 0.3 dttp[c] -> dna[c] + ppi[c]');

% rna is not in Recon2
model = addReaction(model,'rna','reactionFormula','0.18 atp[c] + 0.3 ctp[c] + 0.34 gtp[c] + 0.18 utp[c] -> ppi[c] + rna[c]');
     
% m90191[c] corresponds to human free essential aminoacids
model = addReaction(model,'free_essential_aminoacids','reactionFormula','0.05 ile_L[c] + 0.16 leu_L[c] + 0.14 lys_L[c] + 0.03 met_L[c] + 0.05 phe_L[c] + 0.3 thr_L[c] + 0.01 trp_L[c] + 0.08 tyr_L[c] + 0.18 val_L[c] -> m90191[c]');
 
%m90192[c] corresponds to human free nonessential aminoacid storage
model = addReaction(model,'free_nonessential_aminoacids','reactionFormula','0.15 ala_L[c] + 0.01 arg_L[c] + 0.01 asn_L[c] + 0.3 asp_L[c] + 0.01 cys_L[c] + 0.08 glu_L[c] + 0.25 gln_L[c] + 0.11 gly[c] + 0.03 his_L[c] + 0.01 pro_L[c] + 0.04 ser_L[c] -> m90192[c]');

%m90193[c] corresponds to human bile acid synthesis
model = addReaction(model,'bile_acid_synthesis','reactionFormula','0.0007 m00745[c] + 0.0004 cholate[c] + 0.32 dgchol[c] + 0.163 gchola[c] + 0.17 m01989[c] + 0.0147 m02000[c] + 0.015 m02004[c] + 0.17 tdchola[c] + 0.067 tchola[c] + 0.067 m02964[c] + 0.008 m02965[c] + 0.004 m02966[c] + 0.0002 m90173[c] -> m90193[c]');

    %m00745 is 3alpha,12alpha-dihydroxy-5beta-cholanate
    model = addReaction(model,'m00745','reactionFormula','cholate[c] + h2o[c] + nadp[c] -> m00745[c] + h[c] + nadph[c] + o2[c]');

    %m01989 (glycodeoxycholate) is not in Recon2
    model = addReaction(model,'m01989','reactionFormula','m00745[c] + gly[c] + h[c] <=> m01989[c] + h2o[c]');

    %m02000 is not in Recon2 
    model = addReaction(model,'3B_hydroxy_5_cholestenal','reactionFormula','xol27oh[m] + nadp[c] <=> 3B_hydroxy_5_cholestenal[c] + h[c] + nadph[c]');
    model = addReaction(model,'lithocholate','reactionFormula','3B_hydroxy_5_cholestenal[c] + h2o[c] + nadp[c] -> 2 h[c] + lithocholate[c] + nadph[c]');
    model = addReaction(model,'m02000','reactionFormula','amp[c] + m02000[c] + h[c] + ppi[c] <=> atp[c] + gly[c] + lithocholate[c]');

    %m02004 is not in Recon2
    model = addReaction(model,'7ketolitocholate','reactionFormula','C02528[c] + h[c] + nadph[c] + o2[c] -> 7ketolitocholate[c] + 2 h2o[c] + nadp[c]');
    model = addReaction(model,'ursodeoxycholate','reactionFormula','7ketolitocholate[c] + h[c] + nadph[c] -> nadp[c] + ursodeoxycholate[c]');
    model = addReaction(model,'ursodeoxycholylCoA','reactionFormula','atp[c] + coa[c] + ursodeoxycholate[c] -> amp[c] + ppi[c] + ursodeoxycholylCoA[c]');
    model = addReaction(model,'m02004','reactionFormula','gly[c] + ursodeoxycholylCoA[c] -> coa[c] + m02004[c]');
    
    %m02964 is not in Recon2 
    model = addReaction(model,'m02964','reactionFormula','coa[c] + h[c] + m02964[c] <=> taur[c] + deoxycholoylCoA[c]');
    model = addReaction(model,'deoxycholylCoA','reactionFormula','m00745[c] + atp[c] + coa[c] <=> amp[c] + deoxycholoylCoA[c] + ppi[c]');
    
    %m02965 is not in Recon2
    model = addReaction(model,'m02965','reactionFormula','amp[c] + h[c] + ppi[c] + m02965[c] <=> atp[c] + lithocholate[c] + taur[c]');
    
    %m02966 is not in Recon2 
    model = addReaction(model,'m02966','reactionFormula','taur[c] + ursodeoxycholylCoA[c] -> coa[c] + m02966[c]');
    
    %m90173 is not in Recon2. It is produced by: 
    model = addReaction(model,'m90173','reactionFormula','gly[c] + m90185[c] -> coa[c] + m90173[c]');
    
    model = addReaction(model,'m90185','reactionFormula',' -> m90185[c]');
    
    model = addReaction(model,'outm90173','reactionFormula','m90173[c] -> ');
  
% m90194 is not in Recon2
model = addReaction(model,'m90194','reactionFormula','udpg[c] -> udp[c] + m90194[c]');

% m90195 is not in Recon2
model = addReaction(model,'m90195','reactionFormula','0.003 xolest2_hs[c] + 0.374 pchol_hs[c] + 0.25 pe_hs[c] + 0.106 pail_hs[c] + 0.082 ps_hs[c] + 0.081 sphmyln_hs[c] + 0.104 m02958[c] -> m90195[c]');

    % m02958 
    model = addReaction(model,'acoa','reactionFormula','crm_hs[c] + coa[c] + h[c] -> sphings[c] + acoa[c]');
    model = addReaction(model,'m02958','reactionFormula','dag_hs[c] + acoa[c] -> coa[c] + m02958[c]');
    
%% biomass demand
%model = addReaction(model,'biomass_demand','reactionFormula','m99999[c] -> ');
model = addReaction(model,'clpn_out','reactionFormula','clpn_hs[c] -> ');

model=addReaction(model, 'DM_gudac_c_', 'gudac[c] <=> '); % Add this reaction, present in the experimental data
model.subSystems(end)=model.subSystems{end};
model.rev=zeros(numel(model.rxns,1))
model.rev(model.lb<0 & model.ub>0)=1;
sens=find(model.lb<0 & model.ub==0); %Find uptake reactions
model.S(:,sens)=-model.S(:,sens);
model.ub(sens)= -model.lb(sens);
model.lb(sens)=0;

%% Model reconstruction
% load dico_recon.mat
% load dico_CRC_cell_lines.mat %The medium specification cannot be done here, because metabolites are missing in the medium to be able to have the biomass reaction active
% 
% discretized = discretize_FPKM(fpkmIntestine, colnames);
% upper_bounds_data = readtable('UpperBounds_Zielinski_withoutBiomass.xlsx');
% add_core=upper_bounds_data.Var1;
% 
% Cmodel = model;
% 
% A_final = zeros(numel(Cmodel.rxns),numel(colnames));
% col = cell2table(colnames);
% 
% for i=1:numel(col)   % need to define medium for the cells used here
%     match=ismember(dicoCRCcelllines(:,1),col(i,1));
%     used_medium=table2array(dicoCRCcelllines(match,2))
%     if ~cellfun(@isempty, strfind(used_medium,'RPMI'))
%          load ('RPMI_charstrippedFBS_DHT_MTA_content.mat')
%     elseif ~cellfun(@isempty, strfind(used_medium,'DMEM'))
%          load ('DMEM_medium.mat')
%     elseif ~cellfun(@isempty, strfind(used_medium,'EMEM'))
%          load ('EMEM_medium.mat')
%     elseif ~cellfun(@isempty, strfind(used_medium,'F12K'))
%          load ('F12K_medium.mat')
%     elseif ~cellfun(@isempty, strfind(used_medium,'McCoy'))
%          load ('McCoys_Medium.mat')
%     elseif ~cellfun(@isempty, strfind(used_medium,'Leibovitz'))
%          load ('Leibovitzs_L15_medium.mat')
%     end
%     
% 
% epsilon = 1e-4;
% consensus_proportion = 0.9; %gene has to be expressed in 90% of the cases in order to be included. Only relevant if you want to create one model from different samples
% already_mapped_tag = 0;
% 
% unpenalizedSystems = {'Transport, endoplasmic reticular';
%     'Transport, extracellular';
%     'Transport, golgi apparatus';
%     'Transport, mitochondrial';
%     'Transport, peroxisomal';
%     'Transport, lysosomal';
%     'Transport, nuclear'};
% 
% for a = 1:size(Cmodel.subSystems,1)
%      Cmodel.subSystems{a}=char(Cmodel.subSystems{a}); %I did this becaus otherwise it was not a cell of character arrays
% end
% 
% unpenalized = Cmodel.rxns(ismember(Cmodel.subSystems,unpenalizedSystems));
% not_medium_constrained = 'EX_tag_hs(e)';
% 
% optional_settings.unpenalized = unpenalized;
% optional_settings.func = {'biomass_reaction'};
% 
% optional_settings.not_medium_constrained = not_medium_constrained;
% medium=union(medium, 'asp_L[e]');
% medium=union(medium, 'asn_L[e]');
% optional_settings.medium = medium;
% 
% %% Single models
% [~, A] = fastcormics_RNAseq_core(Cmodel, discretized(:,i), rownames, dico_RECON, ...
%             already_mapped_tag, consensus_proportion, epsilon, optional_settings,add_core);
% A_final(A,i)= 1;
%     
% delete *.log
% poolobj = gcp('nocreate');
% delete(poolobj);
% 
% end
% models_keep_single = A_final;
% 
% save('models_CRC_CCLE_Growth_rate_iHsa_biomass_MP','models_keep_single')

%% Set the new constraints according to Zielinski et al 2017
load models_CRC_CCLE_Growth_rate_iHsa_biomass_MP.mat
optional_settings.func = {'biomass_reaction'};

upper_bounds_data = readtable('UpperBounds_Zielinski_withoutBiomass.xlsx');
upper_bounds = table2array(upper_bounds_data(:,2:end));
reactions = table2array(upper_bounds_data(:,1));
lower_bounds_data = readtable('LowerBounds_Zielinski_withoutBiomass.xlsx');
lower_bounds = table2array(lower_bounds_data(:,2:end));

Cell_lines = erase(colnames,'_LARGE_INTESTINE');
changeCobraSolver('ibm_cplex')

%% Add constraints to carbon exchange reactions 1%
production_rxns = cell(1,size(models_keep_single,2));
modified_bounds = cell(1,size(models_keep_single,2));
problem = cell(1,size(models_keep_single,2));

for i = 1:size(models_keep_single,2)
    cell_line = find(contains(upper_bounds_data.Properties.VariableNames,Cell_lines(i)));
    ind = find(~cellfun(@isempty,regexp(model.rxns,optional_settings.func{1})));
    model_out = removeRxns(model,model.rxns(setdiff(1:numel(model.rxns),find(models_keep_single(:,i))))); % create model based on active reactions
    model_out = changeObjective(model_out,model.rxns(ind));
    
    % Reset the model bounds to avoid problems identifying exchange
    % reactions
    exRxnsInd=find(sum(abs(model_out.S),1)==1); %find exchange rxns
    flip_EX_rxns = exRxnsInd(sum(model_out.S(:,exRxnsInd),1)==1); %exchange rxns written as: -> A
    model_out.S(:,flip_EX_rxns) = -model_out.S(:,flip_EX_rxns);
    tmp = model_out.lb(flip_EX_rxns);
    model_out.lb(flip_EX_rxns) = -model_out.ub(flip_EX_rxns);
    model_out.ub(flip_EX_rxns) = tmp;   
    
    % Constrain all the reactions
    for a = 1:size(reactions,1)
        if ismember(reactions(a,1),model_out.rxns)
            if (lower_bounds(a,cell_line-1)&&upper_bounds(a,cell_line-1)>0) || (upper_bounds(a,cell_line-1)>0&&lower_bounds(a,cell_line-1)==0)
                models_constraint = changeRxnBounds(model_out,reactions(a,1),0,'l');
                models_constraint = changeRxnBounds(models_constraint,reactions(a,1),upper_bounds(a,cell_line-1),'u');
                production_rxns{i} = [production_rxns{i},reactions(a,1)];
            else
                models_constraint = changeRxnBounds(model_out,reactions(a,1),0,'u');
                models_constraint = changeRxnBounds(models_constraint,reactions(a,1),lower_bounds(a,cell_line-1),'l');
            end
        
        model_out = models_constraint;
        modified_bounds{i} = [modified_bounds{i},reactions(a,1)];
        test = optimizeCbModel(models_constraint);
            if test.f<0.009 %Otherwise it is dropping to 0.0067
                met = find(ismember(model_out.rxns,reactions(a,1)));
                metC = find(ismember(model.rxns,reactions(a,1)));
                model_out.ub(met) = model.ub(metC);
                model_out.lb(met) = model.lb(metC);
                problem{i} = [problem{i},reactions(a,1)];
            end
        end
    end
    
    %Find other active EX reactions
    [up_exRxns_iHsaB, up_ex_mets_carbon_iHsaB,uptake_exRxnsInd_iHsaB]=find_uptakeCarbon_EX_Rxns(model_out,optional_settings.func);
    
    %Sum the uptake reactions
    uptake_rxns = find(contains(upper_bounds_data.Var1,setdiff(upper_bounds_data.Var1,production_rxns{i})));
    uptake_rxns(find(ismember(uptake_rxns,(find(ismember(reactions,'EX_o2(e)'))))))=[]; %Remove O2 reaction, since it is not a carbon source
    value = sum(lower_bounds(uptake_rxns(1:end-1),cell_line-1));
    
    % 1% of the sum will be used as lb for other uptake exchange reactions,
    % with no experimental data. 
    test = addMetabolite(model_out,'CarbonAdd');
    test = addReaction(test,'EX_CarbonAdd','reactionFormula',' -> CarbonAdd');
    Ex_rxns_test = find(ismember(test.rxns,'EX_CarbonAdd'));
    test.lb(Ex_rxns_test) = 0;
    test.ub(Ex_rxns_test) = -0.01*value;
    
    reactions_ID=find(ismember(test.rxns, reactions));
    uptake_exRxnsInd_iHsaB=setdiff(uptake_exRxnsInd_iHsaB, reactions_ID);
    
    for z = 1:numel(uptake_exRxnsInd_iHsaB)
        if test.rev(uptake_exRxnsInd_iHsaB(z))==1
            %Uptake reaction
            test.rev(uptake_exRxnsInd_iHsaB(z))=0;
            test.lb(uptake_exRxnsInd_iHsaB(z)) = -1000;
            test.ub(uptake_exRxnsInd_iHsaB(z)) = 0;
            %Secretion reaction
            as = find(test.S(:,uptake_exRxnsInd_iHsaB(z))== -1);
            test.S(:,end+1)=0;
            test.S(as,end)=-1;
            test.rxns(end+1)=strcat('EX_sec_',test.mets(as));
            test.lb(end+1)=0;
            test.ub(end+1)=1000;
            test.rev(end+1)=0;
            test.c(end+1)=0;
        end    
    end
    
    CarbonAdd_met = find(ismember(test.mets,'CarbonAdd'));
    test.S(CarbonAdd_met,uptake_exRxnsInd_iHsaB) = 1;
       
    sol = optimizeCbModel(test);
    optimization_results_iHsaB{i} = sol.f;
    EX_fluxes{i} = sol.x(uptake_exRxnsInd_iHsaB);
    
     %Check why it is not working
    exRxnsInd=find(sum(abs(model_out.S),1)==1);
    non_0_flux = [];
    for c = 1:size(exRxnsInd,2)
        if sol.x(exRxnsInd(c))~=0
            non_0_flux = [non_0_flux,model_out.rxns(exRxnsInd(c))];
        end
    end
    z = find(contains(model_out.rxns,non_0_flux));
    flux_ex_rxns_iHsaB{i} = [model_out.rxns(z),num2cell(sol.x(z))];   
end
save('Growth_rate_Recon2_iHsa_FinalScript','optimization_results_iHsaB','problem')
%% Calculate growht rate based on doubling time
% growht rate = ln(2)/doubling time
annotation = readtable('Doubling_times_OConnor_1997.xlsx');

growth_rate = [];
for i = 1:size(Cell_lines)
    cell_line = find(contains(annotation.CellLine,Cell_lines(i)));
    growth = log(2)/annotation.DoublingTime(cell_line);
    growth_rate = [growth_rate,growth];
end

predicted_rates = [];
for a = 1:size(optimization_results_iHsaB,2)
    predicted_rates = [predicted_rates,optimization_results_iHsaB{a}];
end

figure
sz = 25;
s=scatter(categorical(colnames),growth_rate,sz,'filled')
hold on
s= scatter(categorical(colnames),predicted_rates,sz,'filled')
legend('Measured data','FBA prediction')

%% Create models Recon2 with Recon3 biomass noTrTr reaction
%% Introduce the biomass reaction
load('C:\Users\maria.moscardo\Desktop\Internship\Recon204.mat','CmodelR204')
Cmodel = CmodelR204; 
% We remove the original biomass reaction in Recon2 to avoid overlapping.
reactionR2 = find(~cellfun(@isempty,regexpi(Cmodel.rxns, 'biomass')));
Cmodel_no_Biomass = removeRxns(Cmodel,Cmodel.rxns(reactionR2));

Cmodel = addReaction(Cmodel_no_Biomass,'biomass_reaction','reactionFormula','20.6508 h2o[c] + 20.7045 atp[c] + 0.023315 pail_hs[c] + 0.15446 pchol_hs[c] + 0.055374 pe_hs[c] + 0.020401 chsterol[c] + 0.002914 pglyc_hs[c] + 0.011658 clpn_hs[c] + 0.27519 g6p[c] + 0.005829 ps_hs[c] + 0.017486 sphmyln_hs[c]  -> 20.6508 h[c] + 20.6508 adp[c] + 20.6508 pi[c]');

Cmodel=addReaction(Cmodel, 'DM_gudac_c_', 'gudac[c] <=> '); % Add this reaction, present in the experimental data
Cmodel.subSystems(end)=Cmodel.subSystems{end};
Cmodel.rev=zeros(numel(Cmodel.rxns,1))
Cmodel.rev(Cmodel.lb<0 & Cmodel.ub>0)=1;
sens=find(Cmodel.lb<0 & Cmodel.ub==0); %Find uptake reactions
Cmodel.S(:,sens)=-Cmodel.S(:,sens);
Cmodel.ub(sens)= -Cmodel.lb(sens);
Cmodel.lb(sens)=0;

% load dico_recon.mat
% load dico_CRC_cell_lines.mat
% A_final = zeros(numel(Cmodel.rxns),numel(colnames))
% upper_bounds_data = readtable('UpperBounds_Zielinski_withoutBiomass.xlsx');
% add_core=upper_bounds_data.Var1;
% col = cell2table(colnames);
% discretized = discretize_FPKM(fpkmIntestine, colnames);
% 
% for i=1:numel(col)   % need to define medium for the cells used here
%     match=ismember(dicoCRCcelllines(:,1),col(i,1));
%     used_medium=table2array(dicoCRCcelllines(match,2))
%     if ~cellfun(@isempty, strfind(used_medium,'RPMI'))
%          load ('RPMI_charstrippedFBS_DHT_MTA_content.mat')
%     elseif ~cellfun(@isempty, strfind(used_medium,'DMEM'))
%          load ('DMEM_medium.mat')
%     elseif ~cellfun(@isempty, strfind(used_medium,'EMEM'))
%          load ('EMEM_medium.mat')
%     elseif ~cellfun(@isempty, strfind(used_medium,'F12K'))
%          load ('F12K_medium.mat')
%     elseif ~cellfun(@isempty, strfind(used_medium,'McCoy'))
%          load ('McCoys_Medium.mat')
%     elseif ~cellfun(@isempty, strfind(used_medium,'Leibovitz'))
%          load ('Leibovitzs_L15_medium.mat')
%     end
%     
% 
% epsilon = 1e-4;
% consensus_proportion = 0.9; %gene has to be expressed in 90% of the cases in order to be included. Only relevant if you want to create one model from different samples
% already_mapped_tag = 0;
% 
% unpenalizedSystems = {'Transport, endoplasmic reticular';
%     'Transport, extracellular';
%     'Transport, golgi apparatus';
%     'Transport, mitochondrial';
%     'Transport, peroxisomal';
%     'Transport, lysosomal';
%     'Transport, nuclear'};
% 
% for a = 1:size(Cmodel.subSystems,1)
%      Cmodel.subSystems{a}=char(Cmodel.subSystems{a}); %I did this becaus otherwise it was not a cell of character arrays
% end
% 
% unpenalized = Cmodel.rxns(ismember(Cmodel.subSystems,unpenalizedSystems));
% not_medium_constrained = 'EX_tag_hs(e)';
% 
% optional_settings.unpenalized = unpenalized;
% optional_settings.func = {'biomass_reaction'};
% 
% optional_settings.not_medium_constrained = not_medium_constrained;
% medium=union(medium, 'asp_L[e]');
% medium=union(medium, 'asn_L[e]');
% optional_settings.medium = medium;
% 
% %% Single models
% [~, A] = fastcormics_RNAseq_core(Cmodel, discretized(:,i), rownames, dico_RECON, ...
%             already_mapped_tag, consensus_proportion, epsilon, optional_settings,add_core);
% A_final(A,i)= 1;
%     
% delete *.log
% poolobj = gcp('nocreate');
% delete(poolobj);
% 
% end
% models_keep_single = A_final;
% 
% save('models_CRC_CCLE_Growth_rate_noTrTr_MP','models_keep_single')

%% Set the new constraints according to Zielinski et al 2017
load models_CRC_CCLE_Growth_rate_noTrTr_MP.mat
optional_settings.func = {'biomass_reaction'};

upper_bounds_data = readtable('UpperBounds_Zielinski_withoutBiomass.xlsx');
upper_bounds = table2array(upper_bounds_data(:,2:end));
reactions = table2array(upper_bounds_data(:,1));
lower_bounds_data = readtable('LowerBounds_Zielinski_withoutBiomass.xlsx');
lower_bounds = table2array(lower_bounds_data(:,2:end));

Cell_lines = erase(colnames,'_LARGE_INTESTINE');
changeCobraSolver('ibm_cplex')

production_rxns = cell(1,size(models_keep_single,2));
modified_bounds = cell(1,size(models_keep_single,2));
problem = cell(1,size(models_keep_single,2));

for i = 1:size(models_keep_single,2)
    cell_line = find(contains(upper_bounds_data.Properties.VariableNames,Cell_lines(i)));
    ind = find(~cellfun(@isempty, regexp(Cmodel.rxns,optional_settings.func{1})));
    model_out = removeRxns(Cmodel,Cmodel.rxns(setdiff(1:numel(Cmodel.rxns),find(models_keep_single(:,i))))); % create model based on active reactions
    model_out = changeObjective(model_out,Cmodel.rxns(ind));
    
    % Reset the model bounds to avoid problems identifying exchange
    % reactions
    exRxnsInd=find(sum(abs(model_out.S),1)==1); %find exchange rxns
    flip_EX_rxns = exRxnsInd(sum(model_out.S(:,exRxnsInd),1)==1); %exchange rxns written as: -> A
    model_out.S(:,flip_EX_rxns) = -model_out.S(:,flip_EX_rxns);
    tmp = model_out.lb(flip_EX_rxns);
    model_out.lb(flip_EX_rxns) = -model_out.ub(flip_EX_rxns);
    model_out.ub(flip_EX_rxns) = tmp;

    
    for a = 2:size(reactions,1)
        if ismember(reactions(a,1),model_out.rxns)
            if (lower_bounds(a,cell_line-1)&&upper_bounds(a,cell_line-1)>0) || (upper_bounds(a,cell_line-1)>0&&lower_bounds(a,cell_line-1)==0)
                models_constraint = changeRxnBounds(model_out,reactions(a,1),0,'l');
                models_constraint = changeRxnBounds(models_constraint,reactions(a,1),upper_bounds(a,cell_line-1),'u');
                production_rxns{i} = [production_rxns{i},reactions(a,1)];
            else
                models_constraint = changeRxnBounds(model_out,reactions(a,1),0,'u');
                models_constraint = changeRxnBounds(models_constraint,reactions(a,1),lower_bounds(a,cell_line-1),'l');
            end
        
        model_out = models_constraint;
        modified_bounds{i} = [modified_bounds{i},reactions(a,1)];
        test = optimizeCbModel(models_constraint);
            if test.f==0 & test.stat~=0
                met = find(ismember(model_out.rxns,reactions(a,1)));
                metC = find(ismember(Cmodel.rxns,reactions(a,1)));
                model_out.ub(met) = Cmodel.ub(metC);
                model_out.lb(met) = Cmodel.lb(metC);
                problem{i} = [problem{i},a];
            end
        end
    end
    
    %Find other active EX reactions
    [up_exRxns_Recon2_noTrTr, up_ex_mets_carbon_Recon2_noTrTr,uptake_exRxnsInd_Recon2_noTrTr]=find_uptakeCarbon_EX_Rxns(model_out,optional_settings.func);
    
    %Sum the uptake reactions
    uptake_rxns = find(contains(upper_bounds_data.Var1,setdiff(upper_bounds_data.Var1,production_rxns{i})));
    uptake_rxns(find(ismember(uptake_rxns,(find(ismember(reactions,'EX_o2(e)'))))))=[]; %Remove O2 reaction, since it is not a carbon source
    value = sum(lower_bounds(uptake_rxns,cell_line-1));
    
    test = addMetabolite(model_out,'CarbonAdd');
    test = addReaction(test,'EX_CarbonAdd','reactionFormula',' -> CarbonAdd');
    Ex_rxns_test = find(ismember(test.rxns,'EX_CarbonAdd'));
    test.lb(Ex_rxns_test) = 0;
    test.ub(Ex_rxns_test) = -0.01*value;
    
    reactions_ID=find(ismember(test.rxns, reactions));
    uptake_exRxnsInd_Recon2_noTrTr=setdiff(uptake_exRxnsInd_Recon2_noTrTr, reactions_ID);
  
    for z = 1:numel(uptake_exRxnsInd_Recon2_noTrTr)
        if test.rev(uptake_exRxnsInd_Recon2_noTrTr(z))==1
            %Uptake reaction
            test.rev(uptake_exRxnsInd_Recon2_noTrTr(z))=0;
            test.lb(uptake_exRxnsInd_Recon2_noTrTr(z)) = -1000;
            test.ub(uptake_exRxnsInd_Recon2_noTrTr(z)) = 0;
            %Secretion reaction
            as = find(test.S(:,uptake_exRxnsInd_Recon2_noTrTr(z))== -1);
            test.S(:,end+1)=0;
            test.S(as,end)=-1;
            test.rxns(end+1)=strcat('EX_sec_',test.mets(as));
            test.lb(end+1)=0;
            test.ub(end+1)=1000;
            test.rev(end+1)=0;
            test.c(end+1)=0;
        end    
    end
    
    CarbonAdd_met = find(ismember(test.mets,'CarbonAdd'));
    test.S(CarbonAdd_met,uptake_exRxnsInd_Recon2_noTrTr) = 1;
    
    sol = optimizeCbModel(test);
    optimization_results_Recon2_noTrTr{i} = sol.f;
    
    %Check why it is not working
    exRxnsInd=find(sum(abs(model_out.S),1)==1);
    ex_reactions = [];
    opt_values = [];
    for c = 1:size(exRxnsInd,2)
        if sol.x(exRxnsInd(c))~=0
            ex_reactions = [ex_reactions;model_out.rxns(exRxnsInd(c))];
            opt_values = [opt_values;num2cell(sol.x(exRxnsInd(c)))];
        end
    end
    table_exrxns = [ex_reactions,opt_values];
    flux_ex_rxns_notrtrB{i} = table_exrxns;
    flux_ex_rxns_notrtrB{i} = sol.x(uptake_exRxnsInd_Recon2_noTrTr);
end

save('Growth_rate_Recon2_noTrTr_FinalScript','optimization_results_Recon2_noTrTr','problem')
%% Calculate growht rate based on doubling time
% growht rate = ln(2)/doubling time
annotation = readtable('Doubling_times_OConnor_1997.xlsx');

growth_rate = [];
for i = 1:size(Cell_lines)
    cell_line = find(contains(annotation.CellLine,Cell_lines(i)));
    growth = log(2)/annotation.DoublingTime(cell_line);
    growth_rate = [growth_rate,growth];
end

predicted_rates = [];
for a = 1:size(optimization_results_Recon2_noTrTr,2)
    predicted_rates = [predicted_rates,optimization_results_Recon2_noTrTr{a}];
end

figure
sz = 25;
s=scatter(categorical(colnames),growth_rate,sz,'filled')
hold on
s= scatter(categorical(colnames),predicted_rates,sz,'filled')
legend('Measured data','FBA prediction')

%% Create models Recon2 with Recon3 biomass maintenance reaction
%% Introduce the biomass reaction
load('C:\Users\maria.moscardo\Desktop\Internship\Recon204.mat','CmodelR204')
Cmodel = CmodelR204; 
% We remove the original biomass reaction in Recon2 to avoid overlapping.
reactionR2 = find(~cellfun(@isempty,regexpi(Cmodel.rxns, 'biomass')));
Cmodel_no_Biomass = removeRxns(Cmodel,Cmodel.rxns(reactionR2));

Cmodel = addReaction(Cmodel_no_Biomass,'biomass_reaction','reactionFormula','20.6508 h2o[c] + 20.7045 atp[c] + 0.38587 glu_L[c] + 0.35261 asp_L[c] + 0.036117 gtp[c] + 0.50563 ala_L[c] + 0.27942 asn_L[c] + 0.046571 cys_L[c] + 0.326 gln_L[c] + 0.53889 gly[c] + 0.39253 ser_L[c] + 0.31269 thr_L[c] + 0.59211 lys_L[c] + 0.35926 arg_L[c] + 0.15302 met_L[c] + 0.023315 pail_hs[c] + 0.039036 ctp[c] + 0.15446 pchol_hs[c] + 0.055374 pe_hs[c] + 0.020401 chsterol[c] + 0.002914 pglyc_hs[c] + 0.011658 clpn_hs[c] + 0.053446 utp[c] + 0.27519 g6p[c] + 0.12641 his_L[c] + 0.15967 tyr_L[c] + 0.28608 ile_L[c] + 0.54554 leu_L[c] + 0.013306 trp_L[c] + 0.25947 phe_L[c] + 0.41248 pro_L[c] + 0.005829 ps_hs[c] + 0.017486 sphmyln_hs[c] + 0.35261 val_L[c]  -> 20.6508 h[c] + 20.6508 adp[c] + 20.6508 pi[c]');

Cmodel=addReaction(Cmodel, 'DM_gudac_c_', 'gudac[c] <=> '); % Add this reaction, present in the experimental data
Cmodel.subSystems(end)=Cmodel.subSystems{end};

Cmodel.rev=zeros(numel(Cmodel.rxns,1))
Cmodel.rev(Cmodel.lb<0 & Cmodel.ub>0)=1;
sens=find(Cmodel.lb<0 & Cmodel.ub==0); %Find uptake reactions
Cmodel.S(:,sens)=-Cmodel.S(:,sens);
Cmodel.ub(sens)= -Cmodel.lb(sens);
Cmodel.lb(sens)=0;
 
% load dico_recon.mat
% load dico_CRC_cell_lines.mat
% A_final = zeros(numel(Cmodel.rxns),numel(colnames));
% upper_bounds_data = readtable('UpperBounds_Zielinski_withoutBiomass.xlsx');
% add_core=upper_bounds_data.Var1;
% col = cell2table(colnames);
% discretized = discretize_FPKM(fpkmIntestine, colnames);
% 
% for i=1:numel(col)   % need to define medium for the cells used here
%     match=ismember(dicoCRCcelllines(:,1),col(i,1));
%     used_medium=table2array(dicoCRCcelllines(match,2))
%     if ~cellfun(@isempty, strfind(used_medium,'RPMI'))
%          load ('RPMI_charstrippedFBS_DHT_MTA_content.mat')
%     elseif ~cellfun(@isempty, strfind(used_medium,'DMEM'))
%          load ('DMEM_medium.mat')
%     elseif ~cellfun(@isempty, strfind(used_medium,'EMEM'))
%          load ('EMEM_medium.mat')
%     elseif ~cellfun(@isempty, strfind(used_medium,'F12K'))
%          load ('F12K_medium.mat')
%     elseif ~cellfun(@isempty, strfind(used_medium,'McCoy'))
%          load ('McCoys_Medium.mat')
%     elseif ~cellfun(@isempty, strfind(used_medium,'Leibovitz'))
%          load ('Leibovitzs_L15_medium.mat')
%     end
%     
% 
% epsilon = 1e-4;
% consensus_proportion = 0.9; %gene has to be expressed in 90% of the cases in order to be included. Only relevant if you want to create one model from different samples
% already_mapped_tag = 0;
% 
% unpenalizedSystems = {'Transport, endoplasmic reticular';
%     'Transport, extracellular';
%     'Transport, golgi apparatus';
%     'Transport, mitochondrial';
%     'Transport, peroxisomal';
%     'Transport, lysosomal';
%     'Transport, nuclear'};
% 
% for a = 1:size(Cmodel.subSystems,1)
%      Cmodel.subSystems{a}=char(Cmodel.subSystems{a}); %I did this becaus otherwise it was not a cell of character arrays
% end
% 
% unpenalized = Cmodel.rxns(ismember(Cmodel.subSystems,unpenalizedSystems));
% not_medium_constrained = 'EX_tag_hs(e)';
% 
% optional_settings.unpenalized = unpenalized;
% optional_settings.func = {'biomass_reaction'};
% 
% optional_settings.not_medium_constrained = not_medium_constrained;
% medium=union(medium, 'asp_L[e]');
% medium=union(medium, 'asn_L[e]');
% optional_settings.medium = medium;
% 
% 
% 
% %% Single models
% [~, A] = fastcormics_RNAseq_core(Cmodel, discretized(:,i), rownames, dico_RECON, ...
%             already_mapped_tag, consensus_proportion, epsilon, optional_settings,add_core);
% A_final(A,i)= 1;
%     
% delete *.log
% poolobj = gcp('nocreate');
% delete(poolobj);
% 
% end
% models_keep_single = A_final;
% 
% save('models_CRC_CCLE_Growth_rate_maintenance_MP','models_keep_single')
%% Set the new constraints according to Zielinski et al 2017
load models_CRC_CCLE_Growth_rate_maintenance_MP.mat
optional_settings.func = {'biomass_reaction'};

upper_bounds_data = readtable('UpperBounds_Zielinski_withoutBiomass.xlsx');
upper_bounds = table2array(upper_bounds_data(:,2:end));
reactions = table2array(upper_bounds_data(:,1));
lower_bounds_data = readtable('LowerBounds_Zielinski_withoutBiomass.xlsx');
lower_bounds = table2array(lower_bounds_data(:,2:end));

Cell_lines = erase(colnames,'_LARGE_INTESTINE');
changeCobraSolver('ibm_cplex')

production_rxns = cell(1,size(models_keep_single,2));
modified_bounds = cell(1,size(models_keep_single,2));
problem = cell(1,size(models_keep_single,2));

for i = 1:size(models_keep_single,2)
    cell_line = find(contains(upper_bounds_data.Properties.VariableNames,Cell_lines(i)));
    ind = find(~cellfun(@isempty, regexp(Cmodel.rxns,optional_settings.func{1})));
    model_out = removeRxns(Cmodel,Cmodel.rxns(setdiff(1:numel(Cmodel.rxns),find(models_keep_single(:,i))))); % create model based on active reactions
    model_out = changeObjective(model_out,Cmodel.rxns(ind));
    
    % Reset the model bounds to avoid problems identifying exchange
    % reactions
    exRxnsInd=find(sum(abs(model_out.S),1)==1); %find exchange rxns
    flip_EX_rxns = exRxnsInd(sum(model_out.S(:,exRxnsInd),1)==1); %exchange rxns written as: -> A
    model_out.S(:,flip_EX_rxns) = -model_out.S(:,flip_EX_rxns);
    tmp = model_out.lb(flip_EX_rxns);
    model_out.lb(flip_EX_rxns) = -model_out.ub(flip_EX_rxns);
    model_out.ub(flip_EX_rxns) = tmp; 
   
    for a = 1:size(reactions,1)
        if ismember(reactions(a,1),model_out.rxns)
            if (lower_bounds(a,cell_line-1)&&upper_bounds(a,cell_line-1)>0) || (upper_bounds(a,cell_line-1)>0&&lower_bounds(a,cell_line-1)==0)
                models_constraint = changeRxnBounds(model_out,reactions(a,1),0,'l');
                models_constraint = changeRxnBounds(models_constraint,reactions(a,1),upper_bounds(a,cell_line-1),'u');
                production_rxns{i} = [production_rxns{i},reactions(a,1)];
            else
                models_constraint = changeRxnBounds(model_out,reactions(a,1),0,'u');
                models_constraint = changeRxnBounds(models_constraint,reactions(a,1),lower_bounds(a,cell_line-1),'l');
            end
        
        model_out = models_constraint;
        modified_bounds{i} = [modified_bounds{i},reactions(a,1)];
        test = optimizeCbModel(models_constraint);
            if test.f==0 & test.stat~=0
                met = find(ismember(model_out.rxns,reactions(a,1)));
                metC = find(ismember(Cmodel.rxns,reactions(a,1)));
                model_out.ub(met) = Cmodel.ub(metC);
                model_out.lb(met) = Cmodel.lb(metC);
                problem{i} = [problem{i},a];
            end
        end
    end
    
    %Find other active EX reactions
    [up_exRxns_maintenance, up_ex_mets_carbon_maintenance,uptake_exRxnsInd_maintenance]=find_uptakeCarbon_EX_Rxns(model_out,optional_settings.func);
    
    %Sum the uptake reactions
    uptake_rxns = find(contains(upper_bounds_data.Var1,setdiff(upper_bounds_data.Var1,production_rxns{i})));
    uptake_rxns(ismember(uptake_rxns,(find(ismember(reactions,'EX_o2(e)')))))=[]; %Remove O2 reaction, since it is not a carbon source
    value = sum(lower_bounds(uptake_rxns(1:end-1),cell_line-1));
    
    
    % 1% of the sum will be used as lb for other uptake exchange reactions,
    % with no experimental data. 
    test = addMetabolite(model_out,'CarbonAdd');
    test = addReaction(test,'EX_CarbonAdd','reactionFormula',' -> CarbonAdd');
    Ex_rxns_test = find(ismember(test.rxns,'EX_CarbonAdd'));
    test.lb(Ex_rxns_test) = 0;
    test.ub(Ex_rxns_test) = -0.01*value;
    
    reactions_ID=find(ismember(test.rxns, reactions));
    uptake_exRxnsInd_maintenance=setdiff(uptake_exRxnsInd_maintenance, reactions_ID);
  
    for z = 1:numel(uptake_exRxnsInd_maintenance)
        if test.rev(uptake_exRxnsInd_maintenance(z))==1
            %Uptake reaction
            test.rev(uptake_exRxnsInd_maintenance(z))=0;
            test.lb(uptake_exRxnsInd_maintenance(z)) = -1000;
            test.ub(uptake_exRxnsInd_maintenance(z)) = 0;
            %Secretion reaction
            as = find(test.S(:,uptake_exRxnsInd_maintenance(z))== -1);
            test.S(:,end+1)=0;
            test.S(as,end)=-1;
            test.rxns(end+1)=strcat('EX_sec_',test.mets(as));
            test.lb(end+1)=0;
            test.ub(end+1)=1000;
            test.rev(end+1)=0;
            test.c(end+1)=0;
        end    
    end
    CarbonAdd_met = find(ismember(test.mets,'CarbonAdd'));
    test.S(CarbonAdd_met,uptake_exRxnsInd_maintenance) = 1;
    
    sol = optimizeCbModel(test);
    optimization_results_maintenance{i} = sol.f;
    
    
    %Check why it is not working
    exRxnsInd=find(sum(abs(test.S),1)==1);
    ex_reactions = [];
    opt_values = [];
    for c = 1:size(exRxnsInd,2)
        if sol.x(exRxnsInd(c))~=0
            ex_reactions = [ex_reactions;test.rxns(exRxnsInd(c))];
            opt_values = [opt_values;num2cell(sol.x(exRxnsInd(c)))];
        end
    end
    table_exrxns = [ex_reactions,opt_values];
    flux_ex_rxns_maintenance{i} = table_exrxns;
end

save('Growth_rate_Recon2_maintenance_FinalScript','optimization_results_maintenance')

%% Calculate growht rate based on doubling time
% growht rate = ln(2)/doubling time
annotation = readtable('Doubling_times_OConnor_1997.xlsx');

growth_rate = [];
for i = 1:size(Cell_lines)
    cell_line = find(contains(annotation.CellLine,Cell_lines(i)));
    growth = log(2)/annotation.DoublingTime(cell_line);
    growth_rate = [growth_rate,growth];
end

predicted_rates = [];
for a = 1:size(optimization_results_maintenance,2)
    predicted_rates = [predicted_rates,optimization_results_maintenance{a}];
end

figure
sz = 25;
s=scatter(categorical(colnames),growth_rate,sz,'filled')
hold on
s= scatter(categorical(colnames),predicted_rates,sz,'filled')
ylim([0.005 0.055])
legend('Measured data','FBA prediction')

%% Integrate HMR biomass reaction in Recon2
load Recon204.mat
% Cmodel = readCbModel('HMRdatabase2_00_Cobra.xml','fileType','SBML');
% %Biomass_componenets is the biomass reaction present in HMR

Cmodel = CmodelR204;
model = removeRxns(Cmodel,'biomass_reaction');
model = addReaction(model,'biomass_components','reactionFormula','ala_L[c] + arg_L[c] + asn_L[c] + asp_L[c] + chsterol[c] + xolest2_hs[c] + clpn_hs[c] + cys_L[c] + dna[c] + dna5mtc[c] + glu_L[c] + gln_L[c] + gly[c] + his_L[c] + ile_L[c] + leu_L[c] + m02392c[c] + lys_L[c] + met_L[c] + phe_L[c] + pa_hs[c] + pail_hs[c] + pro_L[c] + rna[c] + ser_L[c] + sphmyln_hs[c] + thr_L[c] + trp_L[c] + tyr_L[c] + val_L[c] + glycogen[c] -> ');
model = addReaction(model,'glycogen','reactionFormula','udpg[c] -> h[c] + udp[c] + glycogen[c]');
model = addReaction(model,'dna','reactionFormula','0.3 datp[c] + 0.2 dctp[c] + 0.2 dgtp[c] + 0.3 dttp[c] -> dna[c] + ppi[c]');
model = addReaction(model,'rna','reactionFormula','0.18 atp[c] + 0.3 ctp[c] + 0.34 gtp[c] + 0.18 utp[c] -> ppi[c] + rna[c]');
model = addReaction(model,'dna5mtc','reactionFormula','amet[c] + dna[c] -> dna5mtc[c] + h[c] + ahcys[c]');
model = addReaction(model,'lipid_droplet','reactionFormula','0.19 dag_hs[c] + 0.0014 m00511c[c] + 0.0024 ak2gchol_hs[c] + 0.0006 lpchol_hs[c] + 0.005 chsterol[c] + 0.34 xolest2_hs[c] + 0.0092 pchol_hs[c] + 0.0034 pe_hs[c] + 0.0016 pail_hs[c] + 0.0002 ps_hs[c] + 0.0004 sphmyln_hs[c] + 0.44 m02958c[c] + 0.005 Rtotal[c] -> m02392c[c]'); 
model = addReaction(model,'acyl_PE_Pool','reactionFormula','pe_hs[c] + h2o[c] -> m00511c[c] + h[c] + Rtotal[c]');
model = addReaction(model,'ak2gchol_hs','reactionFormula','ak2lgchol_hs[c] + acoa[c] -> ak2gchol_hs[c] + coa[c]');
model = addReaction(model,'acoa','reactionFormula','crm_hs[c] + coa[c] + h[c] -> sphings[c] + acoa[c]');
model = addReaction(model,'TAG_LD_pool','reactionFormula','acoa[c] + dag_hs[c] -> m02958c[c] + coa[c]');

% Set bounds to avoid having an inconsistent core model
model.rev=zeros(numel(model.rxns),1);
model.rev(model.lb <0 & model.ub> 0)=1;
Irr=(model.lb >=0 & model.ub>0| model.ub<=0 & model.lb<0);
model.rev(Irr)=0;

FakeIrr= model.ub<=0 & model.lb<0;
model.S(:, FakeIrr)= -model.S(:,FakeIrr);
model.ub(FakeIrr)= -model.lb(FakeIrr);
model.lb(FakeIrr)= zeros(sum(FakeIrr),1);

for a = 1:size(model.subSystems,1)
     model.subSystems{a}=char(model.subSystems{a}); %I did this becaus otherwise it was not a cell of character arrays
end

A_final = zeros(numel(model.rxns),numel(colnames));
col = cell2table(colnames);
discretized = discretize_FPKM(fpkmIntestine, colnames);

load ('RPMI_HMR_Recon2.mat')
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
optional_settings.func = {'biomass_components'};

optional_settings.not_medium_constrained = not_medium_constrained;
optional_settings.medium = medium;

upper_bounds_data = readtable('UpperBounds_Zielinski_withoutBiomass.xlsx');
upper_bounds = table2array(upper_bounds_data(:,2:end));

add_core=upper_bounds_data.Var1;
for i=1:numel(col)   % need to define medium for the cells used here   
%% Single models
[~, A] = fastcormics_RNAseq_core(model, discretized(:,i), rownames, dico_RECON, ...
            already_mapped_tag, consensus_proportion, epsilon, optional_settings,add_core);
A_final(A,i)= 1;
    
delete *.log
poolobj = gcp('nocreate');
delete(poolobj);

end
models_keep_single = A_final;

save('models_CRC_CCLE_Growth_rate_HMR_MP','models_keep_single')
%% Set the new constraints according to Zielinski et al 2017
load models_CRC_CCLE_Growth_rate_HMR_MP.mat
optional_settings.func = {'biomass_components'};
Cmodel = model;
upper_bounds_data = readtable('UpperBounds_Zielinski_withoutBiomass.xlsx');
upper_bounds = table2array(upper_bounds_data(:,2:end));
reactions = table2array(upper_bounds_data(:,1));
lower_bounds_data = readtable('LowerBounds_Zielinski_withoutBiomass.xlsx');
lower_bounds = table2array(lower_bounds_data(:,2:end));

Cell_lines = erase(colnames,'_LARGE_INTESTINE');
changeCobraSolver('ibm_cplex')

production_rxns = cell(1,size(models_keep_single,2));
modified_bounds = cell(1,size(models_keep_single,2));
problem = cell(1,size(models_keep_single,2));

for i = 1:size(models_keep_single,2)
    cell_line = find(contains(upper_bounds_data.Properties.VariableNames,Cell_lines(i)));
    ind = find(~cellfun(@isempty, regexp(Cmodel.rxns,optional_settings.func{1})));
    model_out = removeRxns(Cmodel,Cmodel.rxns(setdiff(1:numel(Cmodel.rxns),find(models_keep_single(:,i))))); % create model based on active reactions
    model_out = changeObjective(model_out,Cmodel.rxns(ind));
    
    % Reset the model bounds to avoid problems identifying exchange
    % reactions
    exRxnsInd=find(sum(abs(model_out.S),1)==1); %find exchange rxns
    flip_EX_rxns = exRxnsInd(sum(model_out.S(:,exRxnsInd),1)==1); %exchange rxns written as: -> A
    model_out.S(:,flip_EX_rxns) = -model_out.S(:,flip_EX_rxns);
    tmp = model_out.lb(flip_EX_rxns);
    model_out.lb(flip_EX_rxns) = -model_out.ub(flip_EX_rxns);
    model_out.ub(flip_EX_rxns) = tmp; 
   
    for a = 1:size(reactions,1)
        if ismember(reactions(a,1),model_out.rxns)
            if (lower_bounds(a,cell_line-1)&&upper_bounds(a,cell_line-1)>0) || (upper_bounds(a,cell_line-1)>0&&lower_bounds(a,cell_line-1)==0)
                models_constraint = changeRxnBounds(model_out,reactions(a,1),0,'l');
                models_constraint = changeRxnBounds(models_constraint,reactions(a,1),upper_bounds(a,cell_line-1),'u');
                production_rxns{i} = [production_rxns{i},reactions(a,1)];
            else
                models_constraint = changeRxnBounds(model_out,reactions(a,1),0,'u');
                models_constraint = changeRxnBounds(models_constraint,reactions(a,1),lower_bounds(a,cell_line-1),'l');
            end
        
        model_out = models_constraint;
        modified_bounds{i} = [modified_bounds{i},reactions(a,1)];
        test = optimizeCbModel(models_constraint);
            if test.f<=0.005 | test.stat==0
                met = find(ismember(model_out.rxns,reactions(a,1)));
                metC = find(ismember(Cmodel.rxns,reactions(a,1)));
                model_out.ub(met) = Cmodel.ub(metC);
                model_out.lb(met) = Cmodel.lb(metC);
                problem{i} = [problem{i},a];
            end
        end
    end
    
    %Find other active EX reactions
    [up_exRxns_hmr, up_ex_mets_carbon_hmr,uptake_exRxnsInd_hmr]=find_uptakeCarbon_EX_Rxns(model_out,optional_settings.func);
    
    %Sum the uptake reactions
    uptake_rxns = find(contains(upper_bounds_data.Var1,setdiff(upper_bounds_data.Var1,production_rxns{i})));
    uptake_rxns(ismember(uptake_rxns,(find(ismember(reactions,'EX_o2(e)')))))=[]; %Remove O2 reaction, since it is not a carbon source
    value = sum(lower_bounds(uptake_rxns(1:end-1),cell_line-1));
    
    
    % 1% of the sum will be used as lb for other uptake exchange reactions,
    % with no experimental data. 
    test = addMetabolite(model_out,'CarbonAdd');
    test = addReaction(test,'EX_CarbonAdd','reactionFormula',' -> CarbonAdd');
    Ex_rxns_test = find(ismember(test.rxns,'EX_CarbonAdd'));
    test.lb(Ex_rxns_test) = 0;
    test.ub(Ex_rxns_test) = -0.01*value;
    
    reactions_ID=find(ismember(test.rxns, reactions));
    uptake_exRxnsInd_hmr=setdiff(uptake_exRxnsInd_hmr, reactions_ID);
  
    for z = 1:numel(uptake_exRxnsInd_hmr)
        if test.rev(uptake_exRxnsInd_hmr(z))==1
            %Uptake reaction
            test.rev(uptake_exRxnsInd_hmr(z))=0;
            test.lb(uptake_exRxnsInd_hmr(z)) = -1000;
            test.ub(uptake_exRxnsInd_hmr(z)) = 0;
            %Secretion reaction
            as = find(test.S(:,uptake_exRxnsInd_hmr(z))== -1);
            test.S(:,end+1)=0;
            test.S(as,end)=-1;
            test.rxns(end+1)=strcat('EX_sec_',test.mets(as));
            test.lb(end+1)=0;
            test.ub(end+1)=1000;
            test.rev(end+1)=0;
            test.c(end+1)=0;
        end    
    end
    CarbonAdd_met = find(ismember(test.mets,'CarbonAdd'));
    test.S(CarbonAdd_met,uptake_exRxnsInd_hmr) = 1;
    
    sol = optimizeCbModel(test);
    optimization_results_hmr{i} = sol.f;
    
    
    %Check why it is not working
    exRxnsInd=find(sum(abs(test.S),1)==1);
    ex_reactions = [];
    opt_values = [];
    for c = 1:size(exRxnsInd,2)
        if sol.x(exRxnsInd(c))~=0
            ex_reactions = [ex_reactions;test.rxns(exRxnsInd(c))];
            opt_values = [opt_values;num2cell(sol.x(exRxnsInd(c)))];
        end
    end
    table_exrxns = [ex_reactions,opt_values];
    flux_ex_rxns_hmt{i} = table_exrxns;
end

save('Growth_rate_Recon2_hmr_FinalScript','optimization_results_hmr')

%% Calculate growht rate based on doubling time
% growht rate = ln(2)/doubling time
annotation = readtable('Doubling_times_OConnor_1997.xlsx');

growth_rate = [];
for i = 1:size(Cell_lines)
    cell_line = find(contains(annotation.CellLine,Cell_lines(i)));
    growth = log(2)/annotation.DoublingTime(cell_line);
    growth_rate = [growth_rate,growth];
end

predicted_rates = [];
for a = 1:size(optimization_results_hmr,2)
    predicted_rates = [predicted_rates,optimization_results_hmr{a}];
end

figure
sz = 25;
s=scatter(categorical(colnames),growth_rate,sz,'filled')
hold on
s= scatter(categorical(colnames),predicted_rates,sz,'filled')
legend('Measured data','FBA prediction')
