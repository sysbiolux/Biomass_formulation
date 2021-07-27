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

discretized = discretize_FPKM(fpkmIntestine, colnames);

%% Coefficient test
load models_keep_single_CCLE_generic_biomass_medium

load Recon204.mat
Cmodel = CmodelR204;

Cmodel.S = full(Cmodel.S);
idx = find(~cellfun(@isempty,regexpi(Cmodel.rxns, 'biomass_reaction')));
values = find(Cmodel.S(:,idx)<0); % In order to consider only those components that are to the left of the reaction
coefficients = Cmodel.S(values,idx);
%Cmodel.mets(values(1));
model_keep = Cmodel;
change = [0,0.5,1,5,10];

optional_settings.func = {'biomass_reaction'};
load enrichment_results_genericR2_medium4val_L[c].mat
for i = 1:size(values,1)
     %Enrichment_results = cell(1,1);
    for m = 1:size(change,2) %We want to test how the change in coefficients affect the enrichment. To do that we will apply values from 10x more to 10x less
        model = model_keep;
        
        for n = 1:size(models_keep_single,2)
            ind = find(~cellfun(@isempty, regexp(model.rxns,optional_settings.func{1})));
            model_out = removeRxns(model,model.rxns(setdiff(1:numel(model.rxns),find(models_keep_single(:,n))))); % create model based on active reactions
            model_out = changeObjective(model_out,model.rxns(ind)); % set objective function
            idx = find(~cellfun(@isempty,regexpi(model_out.rxns, 'biomass_reaction')));
            values = find(model_out.S(:,idx)<0);
            model_out.S(values(i),idx) = (model_out.S(values(i),idx)*change(m));
            % Perform single gene deletion
            [grRatio, grRateKO, grRateWT, hasEffect, delRxns, fluxSolution, genelist] = singleGeneDeletionTamara(model_out,'FBA',[],0,1);
            grRatio_biomass(:,1) = grRatio;
        end
        
        for a = 1:size(grRatio_biomass,2)
            index = find(grRatio_biomass(:,a)<= 0.5);

            EG_names = genelist(index);
            [~, ia,ib] = intersect(EG_names,dico_RECON.ENTREZ); 
            EG_names_sym = dico_RECON.SYMBOL(ib);

            %Check enrichment 
            [enrichment] = True_False_positives_Coeff(EG_names_sym,colnames(n),'Recon');
            Enrichment_results{a}(5,:) = enrichment;           
        end
        idx = find(~cellfun(@isempty,regexpi(Cmodel.rxns, 'biomass_reaction')));
        values = find(Cmodel.S(:,idx)<0);
        name= strcat('enrichment_results_genericR2_medium', char(model.mets(values(i))));
        name2= 'Enrichment_results';
        save(name, name2)
    end 
        
end


%% Plot the results
% Here I manually checked which changes were leading to different EG
% predictions. Just arg, were idenitifed.
% Thus we pool them in a common variable
Enrichment_results_arg = Enrichment_results;
Enrichment_results_asn = Enrichment_results;
Enrichment_results_atp = Enrichment_results;
Enrichment_results_chsterol = Enrichment_results;
Enrichment_results_clpn = Enrichment_results;
Enrichment_results_pe_hs = Enrichment_results;
Enrichment_results_sphmyln = Enrichment_results;

enrichment_c2bbe1 = array2table(zeros(5,7));
enrichment_c2bbe1.Properties.VariableNames = {'arg_L','asn_L','atp','chsterol','clpn','pe_hs','sphmyln'};
enrichment_c2bbe1.Properties.RowNames = {'0x','0.5x','1','5x','10x'};
enrichment_c2bbe1(:,1) = table((cell2mat(Enrichment_results_arg{1,1}.Known_EG)./cell2mat(Enrichment_results_arg{1,1}.Predicted_EG)).*100);   
enrichment_c2bbe1(:,2) = table((cell2mat(Enrichment_results_asn{1,1}.Known_EG)./cell2mat(Enrichment_results_asn{1,1}.Predicted_EG)).*100);   
enrichment_c2bbe1(:,3) = table((cell2mat(Enrichment_results_atp{1,1}.Known_EG)./cell2mat(Enrichment_results_atp{1,1}.Predicted_EG)).*100);   
enrichment_c2bbe1(:,4) = table((cell2mat(Enrichment_results_chsterol{1,1}.Known_EG)./cell2mat(Enrichment_results_chsterol{1,1}.Predicted_EG)).*100);   
enrichment_c2bbe1(:,5) = table((cell2mat(Enrichment_results_clpn{1,1}.Known_EG)./cell2mat(Enrichment_results_clpn{1,1}.Predicted_EG)).*100);   
enrichment_c2bbe1(:,6) = table((cell2mat(Enrichment_results_pe_hs{1,1}.Known_EG)./cell2mat(Enrichment_results_pe_hs{1,1}.Predicted_EG)).*100);   
enrichment_c2bbe1(:,7) = table((cell2mat(Enrichment_results_sphmyln{1,1}.Known_EG)./cell2mat(Enrichment_results_sphmyln{1,1}.Predicted_EG)).*100);   

genes_metabolic_c2bbe1=(cell2mat(Enrichment_results_asn{1,1}.Metabolic_cancer_genes)./cell2mat(Enrichment_results_asn{1,1}.Metabolic_genes)).*100;;

z = [genes_metabolic_c2bbe1,table2array(enrichment_c2bbe1)];
x = {'0x','0.5x','1','5x','10x'};
cmap = [0 0 0;0 0.4470 0.7410;0.8500 0.3250 0.0980;0.9290 0.6940 0.1250;0.4940 0.1840 0.5560;0.4660 0.6740 0.1880;0.3010 0.7450 0.9330;0.6350 0.0780 0.1840;.2 .6 .5;0, 0.5, 0]%;1, 0, 0;0.75, 0, 0.75;0.75, 0.75, 0;0.95 0.61 0.73;0.77 0.38 0.06;1 0.75 0;0.8 0.58 0.46;0.7 0.75 0.71;0.44 0.31 0.22];    
figure
h = bar(z,'Facecolor','flat')
for k = 1:size(z,2)
    h(k).FaceColor = [cmap(k,:)];
end
set(gca,'XTickLabel',x,'XTick',1:numel(x))
ylabel('Fraction of genes known to be essential (%)')
l = cell(1,size(enrichment_c2bbe1,2));
names = enrichment_c2bbe1.Properties.VariableNames;
l{1}='Metabolic genes (1729 genes)'; l{2}= string(names(1)); l{3}=string(names(2));l{4}=string(names(3)); l{5}=string(names(4));l{6}=string(names(5)); l{7}=string(names(6)),l{8}=string(names(7));%l{8}='CRC patients h2o'; l{9}='CRC patients his_L';l{10}='CRC patients sphmyln';
legend(h,l);
ylim([0,55])

%These plot was repeated for all the cell lines, although differences in
%the numerical values, the plot of this cell line is representative for the
%plots obtained for the remaining cell lines. 
