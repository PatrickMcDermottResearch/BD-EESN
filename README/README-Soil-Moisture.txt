1. Start with the ReadInSM.R file in the Read-in directory to read in the SST and soil moisture data, along with reading in all of the functions used throughout.

2. At the very first line of the file ReadInDeepESN.R you will need to change the source function to point to the file  masterBRNNFunctions_001.R  located in the Functions directory. Note, the masterBRNNFunctions_001.R also reads in the ensembleESNCPPFuncts_001.cpp file. This file is also located in the Functions directory. Thus, the sourceCpp function at the top of the masterBRNNFunctions_001.R file must also point to the directory of the ensembleESNCPPFuncts_001.cpp file.

3. The setwd function located in the ReadInSM.R file needs to be pointed to the Datasets directory to read in both the soil moisture and SST data located in this directory.

4. Once these functions are pointing to the correct directory the ReadInSM.R file can be run to read in the data.

5. With the data generated, the D-EESN and BD-EESN model can be run using the Main_Soil_Moisture.R file located in the soil moisture directory.  Under the section “RUN DB-EESN” the sourceCpp function needs to be pointed to the DB-EESNNCPPFunctions.cpp file located in the Functions directory. Note, the parameters used to generate the reservoirs for both the BD-EESN and D-EESN model are selected using the genetic algorithm, using the file named GA_Soil_Moisture.R.
