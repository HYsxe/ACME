# ACME
Supplementary data and codes for the ACME algorithm.

1. Making predictions for peptides: 

    1.1 Description
  
    Reads input peptide sequences and MHC alleles provided by the user and then makes binding predictions. The results are saved to a file.
    
    1.2 Requirements
  
    Python2, Linux system
    
    Python packages: re, keras (version 2.1.4), tensorflow (version 1.5.0), sklearn, math, scipy, numpy, distance
    
    We recommend that you run the software on a server instead of a PC due to the relatively large memory usage.
    
    1.3 Protocol
    
    (1) Download the ACME_codes/ folder. (About 350MB, might take a few minutes)
        
    (2) Change the working directory to this folder. 
        
         Example command: cd /home/user/ACME_codes
        
    (3) Change the path in ACME_codes/binding_prediction.py to the current path of this folder.

         Example: main_dir = "/home/user/ACME_codes/"
         
    (4) Paste the peptide sequences and the corresponding MHC alleles in /binding_prediction/prediction_input.txt
    
         Some examples are shown in ACME_codes/binding_prediction/prediction_input_example.txt
          
    (5) Run binding_prediction.py
    
         Example command: python binding_prediction.py
         
    (6) The results will be saved to ACME_codes/results/binding_prediction.txt
    
	Prediction scores range from 0 to 1. Higher scores indicate higher binding affinities.
	
	Peptides with scores above 0.42 can be considered to be strong binders. 
	
	However, the specific threshold for classification might vary in different experimental setups.
    
2. Repeating the experiments in the ACME paper.

    Please download the ACME_codes/ folder and follow the instructions in the Readme_ACME_repeat.txt file.
 

If you have problems using ACME, please contact hysxe97@gmail.com
  
Yan Hu

School of Life Sciences, Tsinghua University 
     
