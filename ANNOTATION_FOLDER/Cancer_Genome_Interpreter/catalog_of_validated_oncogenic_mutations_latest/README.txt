#
# Catalog of Validated Oncogenic Mutations
# 
# https://www.cancergenomeinterpreter.org/mutations
#

Compilation of mutations in cancer genes that are demonstrated to drive tumor growth or predispose to cancer.
 
This Catalog was compiled by combining the data contained in the DoCM⁠ (PMID:27684579), ClinVar (PMID:26582918)⁠ and OncoKB⁠ (PMID:28890946)
databases, as well as the results of several published experimental assays and additional manual curation⁠ efforts.

We also considered as oncogenic the mutations reported to be biomarkers of tumors response to targeted drugs, which are included in 
the Cancer Biomarkers Database of the Cancer Genome Interpreter. Germline variants found to predispose to cancer, which we retrieved 
from the ClinVar (PMID:26582918) and IARC (PMID:17311302) resources, were also included. The aggregation of the data includes 
(among others) the harmonization of the syntax of variants and the cancer type taxonomy (generically, "cancer" when the specific 
tumor type of the observation is not available) across the different data sources to guarantee the interoperability of all the 
resources that form the Cancer Genome Interpreter.

Contradictory data (i.e. a variant stated as oncogenic and neutral by different resources) was filtered out.

The contents of these resources are licensed under the following terms: 

DoCM (Database of Curated Mutations) license (http://docm.genome.wustl.edu/about),
ClinVar license (https://www.ncbi.nlm.nih.gov/home/about/policies/), 
OncoKB (Precision Oncology Knowledge Base) license (http://oncokb.org/#/terms), 
IARC (International Agency for Research on Cancer) license (http://p53.iarc.fr/Disclaimer.aspx) 
Cancer Biomarkers database license (https://www.cancergenomeinterpreter.org/biomarkers).


FILES:

cancer_acronyms.tsv
----------------------------
This file contains the equivalencies between cancer type acronyms used by the Cancer Genome Interpreter and their common names 
(description column). Both, acronyms and names correspond to the hierarchical organization of cancer types presented in the 
Cancer Genome Interpreter analysis page (https://www.cancergenomeinterpreter.org/analysis). 

catalog_of_validated_oncogenic_mutations.csv
---------------------------------------------------------------
List of mutations known to be involved in tumorigenesis (somatic) or in increasing the predisposition to develop cancer (germline). 
The oncogenic role of each of these mutations has been validated either experimentally or via clinical studies. Each row in the 
file corresponds to a separate oncogenic mutation, with the following columns:

gene: gene symbol
gdna: variant encoded in genomic coordinates and nucleotide change
protein: variant encoded in protein coordinates and aminoacid change (for coding variants)
transcript: Transcript used to map the variant to protein coordinates
info: Consequence of the variant according to ensembl
context: germline or somatic variant
cancer_acronym: acronym of the cancer type (see above description of cancer type acronyms) in which the variant has been observed as oncogenic
source: source of the oncogenic variant (see above, the sources detailed, with their corresponding license)
reference: References (PMIDs) supporting the oncogenic observed variant in the literature.


