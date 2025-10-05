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

📧 weilau@fafu.edu.cn

GitHub: Jdlutt/PlantCellChat
