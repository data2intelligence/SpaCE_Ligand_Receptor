The C++ part of network score computation in the SpaCE framework  
  
Test usage:  
Step 1, Go to the Release folder, and type make.  
Step 2, mkdir data and then copy st.matrix.txt and LRPairFile.txt into the data folder.  
Step 3, Type the following command  
./SpaCE_Ligand_Receptor ../data/st.matrix.txt ../data/LRPairFile.txt ../data/output 1000  

Output will contain the following columns  
Sample: sample ID  
Ratio: (network sum) / (network sum by random) - 1  
Zscore: (network sum - network sum by random) / standard deviation  
Pvalue: Prob(network sum by random >= network sum)  
