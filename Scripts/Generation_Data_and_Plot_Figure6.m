%This script aims to analyse how changes in the stoichiometric coefficients 
%affect the growth rate prediction. Since only 5 CRC cell lines were tested
%for the experimental comparison, the same 5 cell lines will be considered here. 
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

%% Model reconstruction
load models_CRC_CCLE_Growth_rate_MP.mat %The order is 16,17,19,56,22 corresponding to the order in the NCI60 


%% Find the components of the biomass reaction
% We will analyse values from 10x less to 10x more the original
% coefficient. 
decrease = [0:0.1:0.9];
increase = [1:1:10];
change = [decrease,increase];

upper_bounds_data = readtable('UpperBounds_Zielinski_withoutBiomass.xlsx');
upper_bounds = table2array(upper_bounds_data(:,2:end));
reactions = table2array(upper_bounds_data(:,1));
lower_bounds_data = readtable('LowerBounds_Zielinski_withoutBiomass.xlsx');
lower_bounds = table2array(lower_bounds_data(:,2:end));
load Recon204.mat
Cmodel = CmodelR204;
load dico_recon.mat
optional_settings.func = {'biomass_reaction'};

%Change the coefficient of each component in each CRC model 
changeCobraSolver('ibm_cplex')

production_rxns = cell(1,size(models_keep_single,2));
modified_bounds = cell(1,size(models_keep_single,2));
problem = cell(1,size(models_keep_single,2));

for x = 1:size(models_keep_single,2)
    cell_line = find(contains(upper_bounds_data.Properties.VariableNames,erase(colnames(x),'_LARGE_INTESTINE')));
    ind = find(~cellfun(@isempty, regexp(Cmodel.rxns,optional_settings.func{1})));
    model_out = removeRxns(Cmodel,Cmodel.rxns(setdiff(1:numel(Cmodel.rxns),find(models_keep_single(:,x))))); % create model based on active reactions
    model_out = changeObjective(model_out,Cmodel.rxns(ind)); % set objective function
    
    % Find biomass coefficients in each model
    model_out.S = full(model_out.S);
    idx = find(ismember(model_out.rxns,'biomass_reaction'));
    values = find(model_out.S(:,idx)<0); % In order to consider only those components that are to the left of the reaction
    coefficients = model_out.S(values,idx);
    model_out=addReaction(model_out, 'DM_gudac_c_', 'gudac[c] <=> '); % Add this reaction, present in the experimental data
    model_out.subSystems(end)=model_out.subSystems{end};

    model_out.rev=zeros(numel(model_out.rxns,1))
    model_out.rev(model_out.lb<0 & model_out.ub>0)=1;
    sens=find(model_out.lb<0 & model_out.ub==0); %Find uptake reactions
    model_out.S(:,sens)=-model_out.S(:,sens);
    model_out.ub(sens)= -model_out.lb(sens);
    model_out.lb(sens)=0;
    
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
                production_rxns{x} = [production_rxns{x},reactions(a,1)];
            else
                models_constraint = changeRxnBounds(model_out,reactions(a,1),0,'u');
                models_constraint = changeRxnBounds(models_constraint,reactions(a,1),lower_bounds(a,cell_line-1),'l');
            end
        
        model_out = models_constraint;
        modified_bounds{x} = [modified_bounds{x},reactions(a,1)];
        test = optimizeCbModel(models_constraint);
            if test.f==0 & test.stat~=0
                met = find(ismember(model_out.rxns,reactions(a,1)));
                metC = find(ismember(Cmodel.rxns,reactions(a,1)));
                model_out.ub(met) = Cmodel.ub(metC);
                model_out.lb(met) = Cmodel.lb(metC);
                problem{x} = [problem{x},a];
            end
        end
    end
    
     %Find other active EX reactions
    [up_exRxns_Recon, up_ex_mets_carbon_Recon,uptake_exRxnsInd_Recon]=find_uptakeCarbon_EX_Rxns(model_out,optional_settings.func);
   
    %Sum the uptake reactions
    uptake_rxns = find(contains(upper_bounds_data.Var1,setdiff(upper_bounds_data.Var1,production_rxns{x})));
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
    
    %Optimise each time you change the coefficient value
    Final_table_growth = [];
    for i = 1:size(values,1)
        growth_values = [];
        for m = 1:size(change,2) %We want to test how the change in coefficients affect the enrichment. To do that we will apply values from 10x more to 10x less
        model = test;
        
        model.S(values(i),idx) = (model.S(values(i),idx)*change(m));
         
        % Optimise for the biomass
        sol = optimizeCbModel(model);
        growth_values = [growth_values,sol.f];
        end
        Table_growth = array2table(growth_values);
        Final_table_growth = [Final_table_growth;Table_growth];
    end
    Final_table_growth.Properties.RowNames = model.mets(values);
    
    name= strcat('GrowthRate_results_FinalScript', char(colnames(x)));
    name2= 'Final_table_growth';
    save(name, name2)
end

load GrowthRate_results_FinalScriptHCT116_LARGE_INTESTINE.mat
Final_table_growth_HCT116 = Final_table_growth;
% plot results HCT116
annotation = readtable('Doubling_times_OConnor_1997.xlsx');
Cell_lines = erase(colnames,'_LARGE_INTESTINE');

growth_rate = [];
for i = 1:size(Cell_lines)
    cell_line = find(contains(annotation.CellLine,Cell_lines(i)));
    growth = log(2)/annotation.DoublingTime(cell_line);
    growth_rate = [growth_rate,growth];
end

idx_table = [2,5,8,12,13,15,16,17,18,19,21,22,23,24,25,26,27,28,29,30,31,32,34,37,38]; %Components having an impact
cols = [2,4,7,9,11,13,16,20];%Select only a few columns covering the whole range
predicted_rates = [];
for a = 1:size(idx_table,2)
    predicted_rates = [predicted_rates;Final_table_growth_HCT116(idx_table(a),cols)];
end
%predicted_rates([6,9,14],:) = [];

colors = [0, 0.4470, 0.7410;0.8500, 0.3250, 0.0980;0.9290, 0.6940, 0.1250;0.4940, 0.1840, 0.5560;0.4660, 0.6740, 0.1880;0.3010, 0.7450, 0.9330;0.6350, 0.0780, 0.1840;1, 0, 0;0, 0.75, 0.75];
figure
sz = 40;
c='b';
s= scatter(categorical({'Experimental data'}),growth_rate(1),sz,'filled')
s.MarkerFaceColor=[0 0 0];
hold on
sz = 25;
for i = 1:size(cols,2)
c = colors(i,:);
s=scatter(categorical(predicted_rates.Properties.RowNames),table2array(predicted_rates(:,i)),sz,c,'filled')
hold on
end
%s.MarkerFaceColor =[0 0 0;0 0.4470 0.7410;0.8500 0.3250 0.0980;0.9290 0.6940 0.1250;0.4940 0.1840 0.5560;0.4660 0.6740 0.1880;0.3010 0.7450 0.9330;0.6350 0.0780 0.1840;.2 .6 .5;0, 0.5, 0];
legend('Experimental data','0.1xCoefficient','0.3xCoefficient','0.6xCoefficient','0.9xCoefficient','Original Coefficient','3xCofficient','6xCoefficient','10xCoefficient')
xlabel('Metabolite')
ylabel('Growth rate')

% plot results HCT15
load GrowthRate_results_FinalScriptHCT15_LARGE_INTESTINE.mat
Final_table_growth_HCT15 = Final_table_growth;

annotation = readtable('Doubling_times_OConnor_1997.xlsx');
Cell_lines = erase(colnames,'_LARGE_INTESTINE');

growth_rate = [];
for i = 1:size(Cell_lines)
    cell_line = find(contains(annotation.CellLine,Cell_lines(i)));
    growth = log(2)/annotation.DoublingTime(cell_line);
    growth_rate = [growth_rate,growth];
end

idx_table = [2,5,8,12,13,15,16,17,18,19,21,22,23,24,25,26,27,28,29,30,31,32,34,37,38]; %Components having an impact
cols = [2,4,7,9,11,13,16,20];%Select only a few columns covering the whole range
predicted_rates = [];
for a = 1:size(idx_table,2)
    predicted_rates = [predicted_rates;Final_table_growth_HCT15(idx_table(a),cols)];
end
%predicted_rates([6,9,14],:) = [];

colors = [0, 0.4470, 0.7410;0.8500, 0.3250, 0.0980;0.9290, 0.6940, 0.1250;0.4940, 0.1840, 0.5560;0.4660, 0.6740, 0.1880;0.3010, 0.7450, 0.9330;0.6350, 0.0780, 0.1840;1, 0, 0;0, 0.75, 0.75];
figure
sz = 40;
c='b';
s= scatter(categorical({'Experimental data'}),growth_rate(2),sz,'filled')
s.MarkerFaceColor=[0 0 0];
hold on
sz = 25;
for i = 1:size(cols,2)
c = colors(i,:);
s=scatter(categorical(predicted_rates.Properties.RowNames),table2array(predicted_rates(:,i)),sz,c,'filled')
hold on
end
%s.MarkerFaceColor =[0 0 0;0 0.4470 0.7410;0.8500 0.3250 0.0980;0.9290 0.6940 0.1250;0.4940 0.1840 0.5560;0.4660 0.6740 0.1880;0.3010 0.7450 0.9330;0.6350 0.0780 0.1840;.2 .6 .5;0, 0.5, 0];
legend('Experimental data','0.1xCoefficient','0.3xCoefficient','0.6xCoefficient','0.9xCoefficient','Original Coefficient','3xCofficient','6xCoefficient','10xCoefficient')
xlabel('Metabolite')
ylabel('Growth rate')

% plot results HT29
load GrowthRate_results_FinalScriptHT29_LARGE_INTESTINE.mat
Final_table_growth_HT29 = Final_table_growth;

annotation = readtable('Doubling_times_OConnor_1997.xlsx');
Cell_lines = erase(colnames,'_LARGE_INTESTINE');

growth_rate = [];
for i = 1:size(Cell_lines)
    cell_line = find(contains(annotation.CellLine,Cell_lines(i)));
    growth = log(2)/annotation.DoublingTime(cell_line);
    growth_rate = [growth_rate,growth];
end

idx_table = [2,5,8,12,13,15,16,17,18,19,21,22,23,24,25,26,27,28,29,30,31,32,34,37,38]; %Components having an impact
cols = [2,4,7,9,11,13,16,20];%Select only a few columns covering the whole range
predicted_rates = [];
for a = 1:size(idx_table,2)
    predicted_rates = [predicted_rates;Final_table_growth_HT29(idx_table(a),cols)];
end
%predicted_rates([6,9,14],:) = [];

colors = [0, 0.4470, 0.7410;0.8500, 0.3250, 0.0980;0.9290, 0.6940, 0.1250;0.4940, 0.1840, 0.5560;0.4660, 0.6740, 0.1880;0.3010, 0.7450, 0.9330;0.6350, 0.0780, 0.1840;1, 0, 0;0, 0.75, 0.75];
figure
sz = 40;
c='b';
s= scatter(categorical({'Experimental data'}),growth_rate(3),sz,'filled')
s.MarkerFaceColor=[0 0 0];
hold on
sz = 25;
for i = 1:size(cols,2)
c = colors(i,:);
s=scatter(categorical(predicted_rates.Properties.RowNames),table2array(predicted_rates(:,i)),sz,c,'filled')
hold on
end
%s.MarkerFaceColor =[0 0 0;0 0.4470 0.7410;0.8500 0.3250 0.0980;0.9290 0.6940 0.1250;0.4940 0.1840 0.5560;0.4660 0.6740 0.1880;0.3010 0.7450 0.9330;0.6350 0.0780 0.1840;.2 .6 .5;0, 0.5, 0];
legend('Experimental data','0.1xCoefficient','0.3xCoefficient','0.6xCoefficient','0.9xCoefficient','Original Coefficient','3xCofficient','6xCoefficient','10xCoefficient')
xlabel('Metabolite')
ylabel('Growth rate')

% plot results SW620
load GrowthRate_results_FinalScriptSW620_LARGE_INTESTINE.mat
Final_table_growth_SW620 = Final_table_growth;

idx_table = [2,5,8,12,13,15,16,17,18,19,21,22,23,24,25,26,27,28,29,30,31,32,34,37,38]; %Components having an impact
cols = [2,4,7,9,11,13,16,20];%Select only a few columns covering the whole range
predicted_rates = [];
for a = 1:size(idx_table,2)
    predicted_rates = [predicted_rates;Final_table_growth_SW620(idx_table(a),cols)];
end
%predicted_rates([6,9,14],:) = [];

colors = [0, 0.4470, 0.7410;0.8500, 0.3250, 0.0980;0.9290, 0.6940, 0.1250;0.4940, 0.1840, 0.5560;0.4660, 0.6740, 0.1880;0.3010, 0.7450, 0.9330;0.6350, 0.0780, 0.1840;1, 0, 0;0, 0.75, 0.75];
figure
sz = 40;
c='b';
s= scatter(categorical({'Experimental data'}),growth_rate(5),sz,'filled')
s.MarkerFaceColor=[0 0 0];
hold on
sz = 25;
for i = 1:size(cols,2)
c = colors(i,:);
s=scatter(categorical(predicted_rates.Properties.RowNames),table2array(predicted_rates(:,i)),sz,c,'filled')
hold on
end
%s.MarkerFaceColor =[0 0 0;0 0.4470 0.7410;0.8500 0.3250 0.0980;0.9290 0.6940 0.1250;0.4940 0.1840 0.5560;0.4660 0.6740 0.1880;0.3010 0.7450 0.9330;0.6350 0.0780 0.1840;.2 .6 .5;0, 0.5, 0];
legend('Experimental data','0.1xCoefficient','0.3xCoefficient','0.6xCoefficient','0.9xCoefficient','Original Coefficient','3xCofficient','6xCoefficient','10xCoefficient')
xlabel('Metabolite')
ylabel('Growth rate')

% plot results KM12
load GrowthRate_results_FinalScriptKM12_LARGE_INTESTINE.mat
Final_table_growth_KM12 = Final_table_growth;

idx_table = [2,5,8,12,13,15,16,17,18,19,21,22,23,24,25,26,27,28,29,30,31,32,34,37,38]; %Components having an impact
cols = [2,4,7,9,11,13,16,20];%Select only a few columns covering the whole range
predicted_rates = [];
for a = 1:size(idx_table,2)
    predicted_rates = [predicted_rates;Final_table_growth_KM12(idx_table(a),cols)];
end
predicted_rates([6,9,14],:) = [];

colors = [0, 0.4470, 0.7410;0.8500, 0.3250, 0.0980;0.9290, 0.6940, 0.1250;0.4940, 0.1840, 0.5560;0.4660, 0.6740, 0.1880;0.3010, 0.7450, 0.9330;0.6350, 0.0780, 0.1840;1, 0, 0;0, 0.75, 0.75];
figure
sz = 40;
c='b';
s= scatter(categorical({'Experimental data'}),growth_rate(4),sz,'filled')
s.MarkerFaceColor=[0 0 0];
hold on
sz = 25;
for i = 1:size(cols,2)
c = colors(i,:);
s=scatter(categorical(predicted_rates.Properties.RowNames),table2array(predicted_rates(:,i)),sz,c,'filled')
hold on
end
%s.MarkerFaceColor =[0 0 0;0 0.4470 0.7410;0.8500 0.3250 0.0980;0.9290 0.6940 0.1250;0.4940 0.1840 0.5560;0.4660 0.6740 0.1880;0.3010 0.7450 0.9330;0.6350 0.0780 0.1840;.2 .6 .5;0, 0.5, 0];
legend('Experimental data','0.1xCoefficient','0.3xCoefficient','0.6xCoefficient','0.9xCoefficient','Original Coefficient','3xCofficient','6xCoefficient','10xCoefficient')
xlabel('Metabolite')
ylabel('Growth rate')

%% Plot this results also using Recon3 as input model and generic biomass reaction as objective function. 

load GrowthRate_results_FinalScript_R3HCT116_LARGE_INTESTINE.mat
Final_table_growth_HCT116 = Final_table_growth;
% plot results HCT116
annotation = readtable('Doubling_times_OConnor_1997.xlsx');
Cell_lines = erase(colnames,'_LARGE_INTESTINE');

growth_rate = [];
for i = 1:size(Cell_lines)
    cell_line = find(contains(annotation.CellLine,Cell_lines(i)));
    growth = log(2)/annotation.DoublingTime(cell_line);
    growth_rate = [growth_rate,growth];
end

idx_table = [2,7,8,12,14,15,16,17,18,19,20,21,23,24,25,26,27,28,30,31,32,34,35,36,38]; %Components having an impact
cols = [1:8];%Select only a few columns covering the whole range
predicted_rates = [];
for a = 1:size(idx_table,2)
    predicted_rates = [predicted_rates;Final_table_growth_HCT116(idx_table(a),cols)];
end
%predicted_rates([6,9,14],:) = [];

colors = [0, 0.4470, 0.7410;0.8500, 0.3250, 0.0980;0.9290, 0.6940, 0.1250;0.4940, 0.1840, 0.5560;0.4660, 0.6740, 0.1880;0.3010, 0.7450, 0.9330;0.6350, 0.0780, 0.1840;1, 0, 0;0, 0.75, 0.75];
figure
sz = 25;
for i = 1:size(cols,2)
c = colors(i,:);
s=scatter(categorical(predicted_rates.Properties.RowNames),table2array(predicted_rates(:,i)),sz,c,'filled')
hold on
end
%s.MarkerFaceColor =[0 0 0;0 0.4470 0.7410;0.8500 0.3250 0.0980;0.9290 0.6940 0.1250;0.4940 0.1840 0.5560;0.4660 0.6740 0.1880;0.3010 0.7450 0.9330;0.6350 0.0780 0.1840;.2 .6 .5;0, 0.5, 0];
legend('0.1x','0.3x','0.6x','0.9x','1x','3x','6x','10x')
xlabel('Metabolite')
ylabel('Growth rate')

% plot results HCT15
load GrowthRate_results_FinalScript_R3HCT15_LARGE_INTESTINE.mat
Final_table_growth_HCT15 = Final_table_growth;

annotation = readtable('Doubling_times_OConnor_1997.xlsx');
Cell_lines = erase(colnames,'_LARGE_INTESTINE');

growth_rate = [];
for i = 1:size(Cell_lines)
    cell_line = find(contains(annotation.CellLine,Cell_lines(i)));
    growth = log(2)/annotation.DoublingTime(cell_line);
    growth_rate = [growth_rate,growth];
end

idx_table = [2,7,8,12,13,14,15,16,17,18,19,20,21,23,24,25,26,27,28,29,30,31,32,34,35,36,38]; %Components having an impact
cols = [1:8];%Select only a few columns covering the whole range
predicted_rates = [];
for a = 1:size(idx_table,2)
    predicted_rates = [predicted_rates;Final_table_growth_HCT15(idx_table(a),cols)];
end
%predicted_rates([6,9,14],:) = [];

colors = [0, 0.4470, 0.7410;0.8500, 0.3250, 0.0980;0.9290, 0.6940, 0.1250;0.4940, 0.1840, 0.5560;0.4660, 0.6740, 0.1880;0.3010, 0.7450, 0.9330;0.6350, 0.0780, 0.1840;1, 0, 0;0, 0.75, 0.75];
figure
sz = 25;
for i = 1:size(cols,2)
c = colors(i,:);
s=scatter(categorical(predicted_rates.Properties.RowNames),table2array(predicted_rates(:,i)),sz,c,'filled')
hold on
end
%s.MarkerFaceColor =[0 0 0;0 0.4470 0.7410;0.8500 0.3250 0.0980;0.9290 0.6940 0.1250;0.4940 0.1840 0.5560;0.4660 0.6740 0.1880;0.3010 0.7450 0.9330;0.6350 0.0780 0.1840;.2 .6 .5;0, 0.5, 0];
legend('0.1x','0.3x','0.6x','0.9x','1x','3x','6x','10x')
xlabel('Metabolite')
ylabel('Growth rate')

% plot results HT29
load GrowthRate_results_FinalScript_R3HT29_LARGE_INTESTINE.mat
Final_table_growth_HT29 = Final_table_growth;

annotation = readtable('Doubling_times_OConnor_1997.xlsx');
Cell_lines = erase(colnames,'_LARGE_INTESTINE');

growth_rate = [];
for i = 1:size(Cell_lines)
    cell_line = find(contains(annotation.CellLine,Cell_lines(i)));
    growth = log(2)/annotation.DoublingTime(cell_line);
    growth_rate = [growth_rate,growth];
end

idx_table = [2,7,8,12,13,14,15,16,17,18,19,20,21,23,24,25,26,27,28,29,30,31,32,34,35,36,38]; %Components having an impact
cols = [1:8];%Select only a few columns covering the whole range
predicted_rates = [];
for a = 1:size(idx_table,2)
    predicted_rates = [predicted_rates;Final_table_growth_HT29(idx_table(a),cols)];
end
%predicted_rates([6,9,14],:) = [];

colors = [0, 0.4470, 0.7410;0.8500, 0.3250, 0.0980;0.9290, 0.6940, 0.1250;0.4940, 0.1840, 0.5560;0.4660, 0.6740, 0.1880;0.3010, 0.7450, 0.9330;0.6350, 0.0780, 0.1840;1, 0, 0;0, 0.75, 0.75];
figure
sz = 25;
for i = 1:size(cols,2)
c = colors(i,:);
s=scatter(categorical(predicted_rates.Properties.RowNames),table2array(predicted_rates(:,i)),sz,c,'filled')
hold on
end
%s.MarkerFaceColor =[0 0 0;0 0.4470 0.7410;0.8500 0.3250 0.0980;0.9290 0.6940 0.1250;0.4940 0.1840 0.5560;0.4660 0.6740 0.1880;0.3010 0.7450 0.9330;0.6350 0.0780 0.1840;.2 .6 .5;0, 0.5, 0];
legend('0.1x','0.3x','0.6x','0.9x','1x','3x','6x','10x')
xlabel('Metabolite')
ylabel('Growth rate')

% plot results SW620
load GrowthRate_results_FinalScript_R3SW620_LARGE_INTESTINE.mat
Final_table_growth_SW620 = Final_table_growth;

idx_table = [2,7,8,12,13,14,15,16,17,18,19,20,21,23,24,25,26,27,28,29,30,31,32,33,35,36,37,38]; %Components having an impact
cols = [1:8];%Select only a few columns covering the whole range
predicted_rates = [];
for a = 1:size(idx_table,2)
    predicted_rates = [predicted_rates;Final_table_growth_SW620(idx_table(a),cols)];
end
%predicted_rates([6,9,14],:) = [];

colors = [0, 0.4470, 0.7410;0.8500, 0.3250, 0.0980;0.9290, 0.6940, 0.1250;0.4940, 0.1840, 0.5560;0.4660, 0.6740, 0.1880;0.3010, 0.7450, 0.9330;0.6350, 0.0780, 0.1840;1, 0, 0;0, 0.75, 0.75];
figure
sz = 25;
for i = 1:size(cols,2)
c = colors(i,:);
s=scatter(categorical(predicted_rates.Properties.RowNames),table2array(predicted_rates(:,i)),sz,c,'filled')
hold on
end
%s.MarkerFaceColor =[0 0 0;0 0.4470 0.7410;0.8500 0.3250 0.0980;0.9290 0.6940 0.1250;0.4940 0.1840 0.5560;0.4660 0.6740 0.1880;0.3010 0.7450 0.9330;0.6350 0.0780 0.1840;.2 .6 .5;0, 0.5, 0];
legend('0.1x','0.3x','0.6x','0.9x','1x','3x','6x','10x')
xlabel('Metabolite')
ylabel('Growth rate')

% plot results KM12
load GrowthRate_results_FinalScript_R3KM12_LARGE_INTESTINE.mat
Final_table_growth_KM12 = Final_table_growth;

idx_table = [1,2,7,8,9,12,13,14,16,17,18,19,20,21,23,24,25,26,27,28,29,30,31,32,34,35,36,37,38]; %Components having an impact
cols = [1:8];%Select only a few columns covering the whole range
predicted_rates = [];
for a = 1:size(idx_table,2)
    predicted_rates = [predicted_rates;Final_table_growth_KM12(idx_table(a),cols)];
end

colors = [0, 0.4470, 0.7410;0.8500, 0.3250, 0.0980;0.9290, 0.6940, 0.1250;0.4940, 0.1840, 0.5560;0.4660, 0.6740, 0.1880;0.3010, 0.7450, 0.9330;0.6350, 0.0780, 0.1840;1, 0, 0;0, 0.75, 0.75];
figure
sz = 25;
for i = 1:size(cols,2)
c = colors(i,:);
s=scatter(categorical(predicted_rates.Properties.RowNames),table2array(predicted_rates(:,i)),sz,c,'filled')
hold on
end
%s.MarkerFaceColor =[0 0 0;0 0.4470 0.7410;0.8500 0.3250 0.0980;0.9290 0.6940 0.1250;0.4940 0.1840 0.5560;0.4660 0.6740 0.1880;0.3010 0.7450 0.9330;0.6350 0.0780 0.1840;.2 .6 .5;0, 0.5, 0];
legend('0.1x','0.3x','0.6x','0.9x','1x','3x','6x','10x')
xlabel('Metabolite')
ylabel('Growth rate')
