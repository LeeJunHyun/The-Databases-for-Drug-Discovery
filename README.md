# Databases
**This version is just a draft, it is not sorted, orderd.**

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

* RDKit [G. Landrum, 2006] [Landrum, 2016]
    * Rdkit: Open-source cheminformatics
    * SMILES -> chemical structure graph tool
    * https://www.rdkit.org/


* RNA-seq based expression profiles of genes extracted from TCGA bresast cancer level3 data[Prat Aparicio, 2012]

* Cancer hall-mark gene sets [Liberzon et al., 2015]

* PAM50 molecular subtype<br>
breast cancer subtype scheme. Luminal A, Luminal B, Basal-like, HER2.

* DRUGBANK
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


* TTD database [Chen et al., 2002]
    * therapeutic target database
    * protein and nucleic acid target

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


