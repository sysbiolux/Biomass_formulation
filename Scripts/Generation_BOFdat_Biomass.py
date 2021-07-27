#!/usr/bin/env python
# coding: utf-8

# In[1]:


import os
os.chdir(r"C:\Users\maria.moscardo\Desktop\Internship")
import cobra
import pandas as pd


# In[15]:


# Step1
from BOFdat import step1

# This module adds the macromolecules and calculates their stoichiometric
# coefficients.

#Set parameters based on dry weight composition based on (https://books.google.lu/books?id=pnbmCAAAQBAJ&pg=PA214&lpg=PA214&dq=HeLa+dry+weight&source=bl&ots=eXp9OYI9NY&sig=ACfU3U0nPZsQmvXL50gEr6HsAlVjHiyPoQ&hl=es&sa=X&ved=2ahUKEwi9kOmI2aPoAhUPCewKHX8TDr4Q6AEwBHoECAcQAQ#v=onepage&q=HeLa%20dry%20weight&f=false
dna_weight_fraction = 0.016
rna_weight_fraction = 0.066
protein_weight_fraction = 0.66
lipid_weight_fraction = 0.101

#Give the path to each file as function parameters
#Genome file in BioPython supported format (.faa, .fna) and GenBank file
#also in BioPython supported format (.gb, .gbff)
genome = r"C:\Users\maria.moscardo\Desktop\Internship\BOFdat_data\GCF_000001405.39_GRCh38.p13_genomic.fna"
genbank = r"C:\Users\maria.moscardo\Desktop\Internship\BOFdat_data\GCF_000001405_39_GRCh38_p13_genomic.gbff"

#OMICs data as a 2 column csv file, gene/protein and abundance
transcriptomic = r"C:\Users\maria.moscardo\Desktop\Internship\BOFdat_data\test.csv"
proteomic = r"C:\Users\maria.moscardo\Desktop\Internship\BOFdat_data\Proteomics_Data_BOFdat_IDswithdots.csv"

#Lipidomic abundances and conversion to model identifier
## The lipidomics data was obtained from: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6407004/ and the lipidomic conversion file was built using MetExplore (https://link.springer.com/article/10.1007/s11306-020-01663-5#Fig2)
lipidomic_abundances = r"C:\Users\maria.moscardo\Desktop\Internship\BOFdat_data\Lipid_abundance_BOFdat.csv"
lipidomic_conversion = r"C:\Users\maria.moscardo\Desktop\Internship\BOFdat_data\Lipid_annotation_BOFdat.csv"


model = r"C:\Users\maria.moscardo\Desktop\Internship\BOFdat_data\Recon204_final.xml"


# In[3]:


dna_coefficients = step1.generate_dna_coefficients(genome,model,DNA_WEIGHT_FRACTION=dna_weight_fraction)


# In[4]:


os.chdir(r"C:\Users\maria.moscardo\Desktop\Internship\BOFdat")
from core import rna
rna_coefficients = rna.generate_coefficients(genbank,model,transcriptomic,
                         RNA_WEIGHT_FRACTION=lipid_weight_fraction,
                         rRNA_WEIGHT_FRACTION=0.9,
                         tRNA_WEIGHT_FRACTION=0.05,
                         mRNA_WEIGHT_FRACTION=0.05,
                         identifier='geneID')


# In[17]:


os.chdir(r"C:\Users\maria.moscardo\Desktop\Internship\BOFdat")
from core import protein
protein_coefficients = protein.generate_coefficients(genbank,model,proteomic,PROTEIN_WEIGHT_FRACTION=protein_weight_fraction)


# In[11]:


os.chdir(r"C:\Users\maria.moscardo\Desktop\Internship\BOFdat\core")
import pandas as pd
import warnings
import cobra  
import lipid
R_WEIGHT = 284.486

lipid_coefficients = lipid.generate_coefficients(lipidomic_abundances,lipidomic_conversion,
                     model,
                     lipid_weight_fraction,
                     R_WEIGHT)


# In[28]:


# Due to the lack of data to obtain the growth-assocaited cost as suggested in the BOFdat paper, data found on the literature was used (references: PMID 18629932; 20933095; 5817088; 15903248)
maintenance_cost = {
    'GAM': 45,
    'NGAM': 0
}


# In[32]:


os.chdir(r"C:\Users\maria.moscardo\Desktop\Internship\BOFdat")
from util import update

# Create a biomass reaction with the metabolites and the coefficients obtained so far
xml_model = cobra.io.read_sbml_model(model)
# I had to modify the rxn id for atp hydrolisis to account for hard-coded rxn id required by BOFdat
try: xml_model.reactions.ATPM
except:
    atpm_rxn = xml_model.reactions.ATPasel.copy()
    atpm_rxn.id = "ATPM"
    xml_model.remove_reactions(["ATPasel"])
    xml_model.add_reaction(atpm_rxn)
    
    
bofdat_step1 = update.make_new_BOF(xml_model,False,True,dna_coefficients,rna_coefficients,protein_coefficients,
                    lipid_coefficients,maintenance=maintenance_cost)


# In[33]:


#Save the step1 objective function for use in step2
bofdat_step1.to_csv('data/bofdat_step1.csv')


# In[34]:


# Step 2
# First we obtain the metabolites
import os
os.chdir(r"C:\Users\maria.moscardo\Desktop\Internship\BOFdat")
from core import coenzymes_and_ions
path_to_model = r"C:\Users\maria.moscardo\Desktop\Internship\BOFdat_data\Recon204_final.xml"
selected_metabolites = coenzymes_and_ions.find_coenzymes_and_ions(path_to_model)


# In[35]:


# Then we calculate the coefficients
from BOFdat.util.update import determine_coefficients
WEIGHT_FRACTION = 0.1
model = cobra.io.read_sbml_model(path_to_model)
bd_step2 = determine_coefficients(selected_metabolites,model,WEIGHT_FRACTION)


# In[3]:


import os
os.chdir(r"C:\Users\maria.moscardo\Desktop\Internship\BOFdat")
from util import update

#Get the input from the previous step
path_to_biomass = 'data/bofdat_step1.csv'
bd_step1 = update.convert_to_dictionary(path_to_biomass)
bd_step2.update(bd_step1)
update.save_biomass(bd_step2,'data/bofdat_step2.csv')


# In[4]:


# Step 3 BOFdat
# Generation of the initial population
population_name = r"C:\Users\maria.moscardo\Desktop\Internship\BOFdat_data\Populations\test_pop"
path_to_model = r"C:\Users\maria.moscardo\Desktop\Internship\BOFdat_data\Recon204_final.xml"
base_biomass_path = r"C:\Users\maria.moscardo\Desktop\Internship\BOFdat\data\bofdat_step2.csv"
exp_essentiality_path = r"C:\Users\maria.moscardo\Desktop\Internship\BOFdat_data\Hela_EG.csv"


# In[1]:


# Set gurobi solver
import gurobipy
gurobipy.Model()
import optlang.gurobi_interface as grb
grb.Model()


# In[ ]:


from BOFdat import step3
step3.generate_initial_population(population_name,
                                  path_to_model,
                                  base_biomass_path,
                                  exp_essentiality_path,
                                  number_of_populations=5)


# In[ ]:


# Second part step 3
path_to_model = r"C:\Users\maria.moscardo\Desktop\Internship\BOFdat_data\Recon204_final.xml"
initial_population_path = r"C:\Users\maria.moscardo\Desktop\Internship\BOFdat_data\Populations\test_pop_1.csv"
exp_essentiality_path = r"C:\Users\maria.moscardo\Desktop\Internship\BOFdat_data\Hela_EG.csv"
base_biomass = r"C:\Users\maria.moscardo\Desktop\Internship\BOFdat\data\bofdat_step2.csv"
logbook_name = r"C:\Users\maria.moscardo\Desktop\Internship\BOFdat_data\Outputs\logbook_test_pop_1.csv"
hall_of_fame_name = r"C:\Users\maria.moscardo\Desktop\Internship\BOFdat_data\Outputs\hof_test_pop_1.csv"


# In[ ]:


import os
os.chdir(r"C:\Users\maria.moscardo\Desktop\Internship")
import cobra
import pandas as pd

from BOFdat import step3
step3.find_metabolites(path_to_model,
                       initial_population_path,
                       exp_essentiality_path,
                       base_biomass=True,logbook=True,history=False,processes=None,
                       logbook_name=logbook_name,hall_of_fame_name=hall_of_fame_name)


# In[ ]:


#Clustering Step 3
outpath = r"C:\Users\maria.moscardo\Desktop\Internship\BOFdat_data\Outputs"
path_to_model = r"C:\Users\maria.moscardo\Desktop\Internship\BOFdat_data\Recon204_final2.xml"
eps = 10
CONNECTIVITY_THRESHOLD = 15
BASELINE=0.2
show_frequency = True
show_matrix = True
THRESHOLD = 0.2


# In[ ]:


import cobra
import pandas as pd

from BOFdat import step3

selected_metabolites = step3.cluster_metabolites(data_path,
                              path_to_model,
                              CONNECTIVITY_THRESHOLD,
                              BASELINE,
                              eps,
                              show_frequency,
                              show_matrix,
                              frequency_fig_name=os.path.join(outpath,'frequency_fig.svg'),
                              matrix_fig_name=os.path.join(outpath,'matrix_fig.svg'))


# In[ ]:


#from BOFdat.util.update import determine_coefficients
WEIGHT_FRACTION = 0.05
model = cobra.io.read_sbml_model(path_to_model)
bd_step3 = determine_coefficients(list_of_metab,model,WEIGHT_FRACTION)


# In[ ]:


from BOFdat.util import update
path_to_biomass = r"C:\Users\maria.moscardo\Desktop\Internship\BOFdat\data\bofdat_step2.csv"
bd_step2 = update.convert_to_dictionary(path_to_biomass)
bd_step3.update(bd_step2)
update.save_biomass(bd_step3,r"C:\Users\maria.moscardo\Desktop\Internship\BOFdat\data\bofdat_step3.csv")

