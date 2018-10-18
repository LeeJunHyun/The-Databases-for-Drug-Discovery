**This version is just a draft, it is not ordered, categorized, completed yet.**
**If you know other databases, plz notice me. PR or Issues are always available.**
# Databases

* ZINC database [J. J. Irwin et al., 2012] [Sterling and Irwin, 2015] [Kusner et al., 2017]
    * molecule dataset
    * 250,000 drug like commercially abailable molecules
    * 35 million commercially-available compounds
    * maximum atom number 38
    * paper: Zinc: a free tool to discover chemistry for biology [J. J. Irwin et al., 2012]
    * paper: ZINC 15 – Ligand Discovery for Everyone [Sterling and Irwin, 2015]
    * http://zinc.docking.org/
    * http://zinc15.docking.org/

* Connectivity Map [Subramanian A, et al., 2017] [Lamb J, et al., 2006]
    * A Next Generation Connectivity Map: L1000 Platform And The First 1,000,000 Profiles [Subramanian A, et al.]
    * The Connectivity Map: using gene-expression signatures to connect small molecules, genes, and disease. [Lamb J, et al.]
    * https://clue.io/cmap

* PubChem

* Protein Data Bank [Helen Berman et al., 2003] [Kleywegt GJ 2018]
    * information about the 3D structures of proteins, nucleic acids, and complex assemblies.
    * http://www.wwpdb.org/
    * https://www.ebi.ac.uk/pdbe/
    * https://www.rcsb.org/
    * https://pdbj.org/

* GEO [Tanya Barrett, 2013]
    * Gene Expression Omnibus
    * international public repository that archives and freely distributes microarray, next-generation sequencing, and other forms of high-throughput functional genomics data submitted by the research community.
    * https://www.ncbi.nlm.nih.gov/geo/

* PharmGKB [M. Whirl-Carrillo et al., 2012]
    * PHARMACOGENOMICS. KNOWLEDGE BASE.
    * knowledge about the impact of genetic variation on drug response
    * relationship between genetic variations and how our body responds to medications.
    * Drugs, Pathways, Dosing Guidelines, Drug Labels
    * https://www.pharmgkb.org/

* STITCH
    * Chemical-Protein Interaction Networks
    * ORGANISMS 2031, CHEMICALS 0.5 mio, PROTEINS 9.6 mio, INTERACTIONS 1.6 bn
    * 68,000 different chemicals, including 2200 drugs, and connects them to 1.5 million genes across 373 genomes.
    * http://stitch.embl.de/

* RDKit [G. Landrum, 2006] [Landrum, 2016]
    * Rdkit: Open-source cheminformatics
    * SMILES -> chemical structure graph tool
    * https://www.rdkit.org/


* RNA-seq based expression profiles of genes extracted from TCGA bresast cancer level3 data[Prat Aparicio, 2012]

* Cancer hall-mark gene sets [Liberzon et al., 2015]

* PAM50 molecular subtype<br>
breast cancer subtype scheme. Luminal A, Luminal B, Basal-like, HER2.

* DrugBank
    * bioinformatics and cheminformatics resource that combines detailed drug data with comprehensive drug target information.
    * `version 5.1.1`, released 2018-07-03, contains 11,877 drug entries including 2,474 approved small molecule drugs, 1,180 approved biotech (protein/peptide) drugs, 129 nutraceuticals and over 5,748 experimental drugs. Additionally, 5,131 non-redundant protein (i.e. drug target/enzyme/transporter/carrier) sequences are linked to these drug entries. Each DrugCard entry contains more than 200 data fields with half of the information being devoted to drug/chemical data and the other half devoted to drug target or protein data.
    * https://www.drugbank.ca/

* STRING [Szklarczyk et al., 2014]
    * database of putatively associating genes from multiple pieces of evidence like biological experiments, text-mined literature information, computational prediction, etc.
    * ex) protein-protein interaction network topology

* SMILES [Weininger, 1988]
    * Simplified molecular-input line-entry system
    * Drug chemical structure
    * specification in form of a line notation for describing the structure of chemical species using short ASCII strings

* DTIP [Kyle Yingkai Gao et al., 2018]
    * IBM research dataset from BindingDB
    * paper : Interpretable Drug Target Prediction Using Deep Neural Representation
    * 39,747 positive examples and 31,218 negative examples
    * https://github.com/IBM/InterpretableDTIP

* BindingDB [Gilson etal., 2016]
    * Public, web-accessible database
    * binding affinities, focusing chiefly on the interactions of small molecules (drugs/drug candidates) and proteins (targets/target candidates)

* SIDER database [Kuhn et al., 2015]
    * drug side effect
    * 996 drugs and 4192 side effects

* OFFSIDES
    * drugs confounder-controlled side effects
    * 1,332 drugs and 10,093 side effects

* ECFP6
    * extended-connectivity fingerprints with diameter 6
    * drug structural features
    * the hashed 1,024-bit length vector encoding the presence or absence of substructure in a drug molecule

* Gene ontology (GO) annotation [Ashburner et al., 2009]

* Twosides database [Tatonetti et al., 2012]
    * 645 drugs and 1618 DDI(drug-drug interaction), in total 63,473 DDI pairs

* CPI database [Wishart et al., 2008]
    * chemical protein interactome
    * about how much power a drug needs to bind with its protein target

* ChEMBL (StARlite)
    * Chemical European Molecular Biology Laboratory
    * chemical database of bioactive molecules with drug-like properties.
    * 1.8M compounds, 1.1M assays, 69k documents, 12k targets, 11k drugs, 1.7k cells
    * https://www.ebi.ac.uk/chembl/
    * https://www.ebi.ac.uk/chembl/beta/
    * https://chembl.gitbook.io/chembl-interface-documentation/downloads


* TTD database [Chen et al., 2002] [Y. H. Li et al., 2018]
    * therapeutic target database
    * database to provide information about the known and explored therapeutic protein and nucleic acid targets, the targeted disease, pathway information and the corresponding drugs directed at each of these targets. Also included in this database are links to relevant databases containing information about target function, sequence, 3D structure, ligand binding properties, enzyme nomenclature and drug structure, therapeutic class, clinical development status. All information provided are fully referenced.
    * protein and nucleic acid target
    * https://db.idrblab.org/ttd/

* PDBBind dataset [Liu et al., 2017]
    * binding affinities for the protein-ligand complexes in the Protein Data Bank (PDB).
    * `version 2017` released by Jan 1st, 2017. This release provides binding data of a total of 17,900 biomolecular complexes, including protein-ligand (14,761), nucleic acid-ligand (121), protein-nucleic acid (837), and protein-protein complexes (2,181), which is currently the largest collection of this kind.
    * Liu et al., Acc. Chem. Res. 2017, 50, 302-309
    * http://www.pdbbind.org.cn/

* Tox21 Data Challenge 2014
    * for prediction of compounds' interference in biochemical pathways using only chemical structure data(SMILES).
    * https://tripod.nih.gov/tox21/challenge/data.jsp

* GDB Databases
    * GDB-11 [Fink, T et al., 2005, 2007]
        * small organic molecules up to 11 atoms of C, N, O and F following simple chemical stability and synthetic feasibility rules.
    * GDB-13 [Blum L. C. et al., 2009]
        * small organic molecules up to 13 atoms of C, N, O, S and Cl following simple chemical stability and synthetic feasibility rules.
    * GDB-17 [Ruddigkeit Lars et al., 2012]
        * 166.4 billion molecules of up to 17 atoms of C, N, O, S, and halogens.
        * Compared to known molecules in PubChem, GDB-17 molecules are much richer in nonaromatic heterocycles, quaternary centers, and stereoisomers, densely populate the third dimension in shape space, and represent many more scaffold types.
    * http://gdb.unibe.ch/downloads/

* Organic photovoltaics [Hachmann, J. et al., 2014]
    * candidate structures for organic electronic materials in particular photovoltaics
    * Harvard Clean Energy Project.
    * promising compounds that have emerged after studying 2.3 million molecular motifs by means of 150 million density functional theory calculations.

* Protein graph [Dobson, P. and Doig, 2003]
    * 1178 high-resolution proteins in a structurally non-redundant subset of the Protein Data Bank using simple features such as secondary-structure content, amino acid propensities, surface properties and ligands.
    * two functional groupings, enzymes and non-enzymes.
    * nodes are amino acids and two nodes are connected if they are less than 6 Angstroms apart.

* HMDD [Lu M et al., 2008] [Li Y et al., 2014]
    * the Human microRNA Disease Database
    * database that curated experiment-supported evidence for human microRNA (miRNA) and disease associations.
    * `HMDD v3.0`, released June 28 2018, 32281 miRNA-disease association entries which include 1102 miRNA genes, 850 diseases from 17412 papers.
    * http://www.cuilab.cn/hmdd

* DrugTargetCommons [Jing Tang et al., 2018]
    * Drug Target Commons (DTC) is a crowd-sourcing platform to improve the consensus and use of drug-target interactions.
    * https://drugtargetcommons.fimm.fi/
    
[2018.10.18]

|  | Assay annotation | Total |
| :--------: | :--------: | :--------: |
| Compounds | 4,276 | 1,746,997 |
| Targets | 1,007 | 13,023 |
| Publications | 346 | 69,955 |
| Bioactivities | 204,901 | 14,820,874 |


*  IDG Pharos
    * compound and target data resources on public domain
    * Ligand, disease, target
    * https://druggablegenome.net/
    * https://pharos.nih.gov/idg/index


* DDIExtraction2013 [BioNLP Challenge]
    * Extraction of Drug-Drug Interactions from BioMedical Texts
    * Task 1: Recognition and classification of drug names.
    * Task 2: Extraction of drug-drug interactions.
    * https://www.cs.york.ac.uk/semeval-2013/task9/

* Biocreative PPI [BioNLP Kaggle]
    * BioCreAtIvE (Critical Assessment of Information Extraction Systems in Biology)
    * text mining and information extraction systems applied to the biological domain.
    * Gene mention tagging [GM] 
    * Gene normalization [GN] 
    * Extraction of protein-protein interactions from text
    * http://biocreative.sourceforge.net/index.html


* Polysearch2 [Liu Y et al., 2015]
    * online text-mining system for identifying relationships between human diseases, genes, proteins, drugs, metabolites, toxins, metabolic pathways, organs, tissues, subcellular organelles, positive health effects, negative health effects, drug actions, Gene Ontology terms, MeSH terms, ICD-10 medical codes, biological taxonomies and chemical taxonomies.
    * http://polysearch.cs.ualberta.ca/

* SuperTarget [BioNLP]
    * database developed in the first place to collect informations about drug-target relations. It consist mainly of three different types of entities: DRUGS, PROTEINS, SIDE-EFFECTS.
    * database that contains a core dataset of about 7300 drug-target relations of which 4900 interactions have been subjected to a more extensive manual annotation effort. SuperTarget provides tools for 2D drug screening and sequence comparison of the targets. The database contains more than 2500 target proteins, which are annotated with about 7300 relations to 1500 drugs
    * data from DrugBank, BindingDB and SuperCyp
    * http://insilico.charite.de/supertarget/index.php?site=about

* ConsensusPathDB [Kamburov, A. et al., 2013]
    * [2018.10.04] unique physical entities: 170,276, unique interactions: 603,543, gene regulations: 17,410, protein interactions: 397,088, genetic interactions: 1,738, biochemical reactions: 23,482, drug-target interactions: 163,825, pathways: 5,359
    * Data originate from currently 32 public resources for interactions and interactions that we have curated from the literature.
    * http://cpdb.molgen.mpg.de/


* Others
    * http://polysearch.cs.ualberta.ca/otherdatabases

---

# Datasets in Papers

### Hybrid Approach of Relation Network and Localized Graph Convolutional Filtering for Breast Cancer Subtype Classification
Sungmin Rhee, Seokjun Seo, Sun Kim<br>
IJCAI 2018<br>
https://arxiv.org/abs/1711.05859
![dataset-01](/img/dataset-01.png)


### Interpretable Drug Target Prediction Using Deep Neural Representation
Kyle Yingkai Gao, Achille Fokoue, Heng Luo, Arun Iyengar, Sanjoy Dey, Ping Zhang<br>
IJCAI 2018<br>
https://astro.temple.edu/~tua87106/ijcai_dti.pdf
![dataset-02](/img/dataset-02.png)

### Drug Similarity Integration Through Attentive Multi-view Graph Auto-Encoders
Tengfei Ma, Cao Xiao, Jiayu Zhou, FeiWang<br>
IJCAI 2018<br>
https://arxiv.org/abs/1804.10850
![dataset-03](/img/dataset-03.png)



### Graph Convolutional Policy Network for Goal-Directed Molecular Graph Generation
Jiaxuan You, Bowen Liu, Rex Ying, Vijay Pande, Jure Leskovec<br>
NIPS 2018<br>
https://arxiv.org/pdf/1806.02473.pdf
![dataset-04](/img/dataset-04.png)


### Optimizing distributions over molecular space. An Objective-Reinforced Generative Adversarial Network for Inverse-design Chemistry (ORGANIC)
Benjamin Sanchez-Lengeling, Carlos Outeiral, Gabriel L. Guimaraes, Alan Aspuru-Guzik<br>
ChemRxiv e-prints, 8 2017.<br>
https://chemrxiv.org/articles/ORGANIC_1_pdf/5309668
![dataset-05](/img/dataset-05.png)


