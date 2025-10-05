PlantCellChat: An R Toolkit for Predicting Plant Cell–Cell Communication

<img width="865" height="784" alt="image" src="https://github.com/user-attachments/assets/91e80131-d15e-4cd5-8705-1dceadcc0c2b" />

Graphical abstract illustrating the PlantCellChat workflow, including database construction, receptor prediction (PCC-GCN), and communication analysis.

🌱 Introduction

PlantCellChat is an R-based computational framework for inferring and visualizing cell–cell communication (CCC) in plants using single-cell RNA sequencing (scRNA-seq) and spatial transcriptomics data.
It integrates a high-confidence PlantCellChatDB of ligand–receptor pairs and implements a mass action–based model to compute communication strength between cell types.
PlantCellChat additionally incorporates a graph convolutional network (GCN) module, PCC-GCN, to predict hormone-binding receptor proteins that have not yet been experimentally characterized — a feature absent in all existing CCC tools.
<img width="870" height="695" alt="image" src="https://github.com/user-attachments/assets/c2c0cde2-4d65-451d-8a89-a90262597e63" />



⚙️ Installation
# install.packages("devtools");
library(devtools);
install_github("Jdlutt/PlantCellChat");
library(PlantCellChat)


All dependencies will be installed automatically.

🚀 Quick Start Example

# Step 1: Create analysis object
pccob <- CreatePlantCellChat(object = seob,
                      meta = seob@meta.data,
                      input.ident = "celltype",
                      assay = "SCT")

# Step 2: Extract signaling genes and identify overexpressed genes
pccob <- ExtractSignalingData(pccob)
pccob <- IdentifyOverExpressedGenes(pccob)

# Step 3: Detect overexpressed ligand–receptor interactions
pccob <- ExtractOverExpressedInteractions(pccob)
pccob <- CalculateAvgExp(pccob,methods = "average")

# Step 4: Compute communication strength
pccob <- CalculateCommunStrength(pccob,Kh = 0.5,n = 1,num.permutations = 100,seed = 123)
pccob <- CalculateSignalingStrength(pccob,Kh = 0.5,n = 1,num.permutations = 100,seed = 123)
# Step 5: Visualize cell–cell communication network
PlottingLRpairStats(pccob)
PlottingCommunNetwork(pccob,ligand.type = "lrpairs",comm.pattern = "paracrine",plot.type = "chord",input.color = color_list)
CompareSignalCommunStrength(pcc_obj = pccob,
                            pcc_obj_list = pccob_list,
                            ligand.type = "lrpairs",
                            key.source = "Leaf guard cell",
                            key.target = "Mesophyll",
                            key.signal = "JA",
                            top.n = 10,
                            input.color = c("#f68d8d","#88b2f2"))

For a full tutorial, see vignettes/PlantCellChat_Tutorial.Rmd.

🧠 PCC-GCN: Hormone Receptor Prediction Module
🔍 Overview

Many hormone receptors in plants remain uncharacterized. To address this, PlantCellChat integrates a deep learning module — PCC-GCN (PlantCellChat Graph Convolutional Network) — to predict hormone-binding receptor proteins from protein–protein interaction (PPI) networks.

PCC-GCN models the topological structure of hormone-related protein networks and the biochemical properties of each protein to classify receptor functional types (e.g., ABA, BR, GA receptors).

⚙️ Model Architecture

Input layer:
A graph of 684 hormone-related proteins (nodes) and 2,580 interactions (edges) derived from Arabidopsis thaliana.
Each protein node embeds 31 biochemical and structural features, including amino acid composition, isoelectric point, aliphatic index, and secondary structure ratios.

Hidden layer:
A single GCN layer (16 neurons) with LeakyReLU activation is used to aggregate neighboring node features, weighted by TM-score and STRING-score.

Output layer:
A multi-label sigmoid classifier returns the probability that a given protein belongs to one or more hormone receptor classes (ABA, BR, CTK, ET, GA, IAA, JA, SA).

🧩 Input Data

You can provide your own PPI network and node features, or use the built-in Arabidopsis model:

# Load built-in PCC-GCN model
data("PCC_GCN_model")

# Optional: visualize graph structure
PlotPCCNetwork(PCC_GCN_model)

🧪 Predicting Hormone Receptors
# Predict receptor classes for a new protein set
pred_results <- PredictHormoneReceptor(PCC_GCN_model, protein_features)

# View top predicted receptors
head(pred_results)


Output Example:

Protein_ID	ABA	BR	CTK	ET	GA	IAA	JA	SA
AT1G12340	0.91	0.04	0.02	0.01	0.09	0.03	0.02	0.05
AT3G56780	0.05	0.88	0.11	0.03	0.02	0.09	0.07	0.01

Each column represents the probability of binding to a specific hormone signal.

🧾 Interpretation

Probability > 0.8 indicates strong evidence of functional association with that hormone receptor class.

Predictions can guide receptor candidate validation in non-model plants.

Model performance (10-fold CV on Arabidopsis):

AUC: 0.96 AUCPR: 0.91 MCC: 0.82

📊 Key Functions Summary
Function	Description
CreatePlantCellChat()	Create analysis object
ExtractSignalingData()	Extract signaling-related genes
IdentifyOverExpressedGenes()	Identify differentially expressed signaling genes
ExtractOverExpressedInteractions()	Identify active ligand–receptor pairs
CalculateCommunStrength()	Compute CCC strength
PlottingCommunNetwork()	Visualize cell–cell communication
PredictHormoneReceptor()	Predict hormone receptor type via PCC-GCN
PlottingSignalContribution()	Visualize receptor contributions in signaling pathways
🌍 Repository Structure
PlantCellChat/
├── R/                      # Core functions
├── data/                   # Example model (PCC-GCN)
├── vignettes/              # Tutorials (.Rmd)
├── figures/                # Workflow figures
└── README.md               # Documentation

🧠 Citation

Liu W. et al. (2025). PlantCellChat: A computational framework for predicting plant cell–cell communication and hormone receptor networks in plants. Manuscript in preparation.

📨 Contact

📧 weiliu.bioinfo@outlook.com

GitHub: Jdlutt/PlantCellChat
