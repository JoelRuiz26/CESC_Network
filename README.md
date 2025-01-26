# NW_CESC
Networks for cervical cancer by HPV clade from TCGA-CESC project

1. Get data from TCGA

   0_Get_HPV_Clade.R

   1_Get_data_RNA.R

2. Preprocess data 

   2_Prepro_Data_RNA_HPV.R

3. MI_networks

   3_https://github.com/CSB-IG/ARACNE-multicore

   3_run.sh (modification to run 4 cores of ARACNE-multicore/launch)

4. Global_network_Analisis

   4_GiantComponent (threshold for minimum GiantComponent with highest MI)

   5_Network_Topology (Network global metrics)

   6_Louvain_Modularity
