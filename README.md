# ACME
Supplementary data and codes for the ACME algorithm.

1. Making predictions for peptides: 
  1.1 Description
    Reads input peptide sequences and MHC alleles provided by the user and then make binding predictions. The results are saved to a file.
    
  1.2 Requirements
    Python2, Linux system
    Python packages: re, keras, sklearn, math, scipy, numpy
    We recommend that you run the software on a server instead of a PC due to the relatively large memory usage.
    
  1.1 Protocol
    (1) Download the ACME_codes/ folder.
    (2) Change the working directory to this folder. 
          Example command: cd /home/user/ACME_codes
    (3) Change the path in ACME_codes/binding_prediction.py to the current path of this folder.
          Example: main_dir = "/home/user/ACME_codes/"
    (4) Paste the peptide sequences and the corresponding MHC alleles in /binding_prediction/prediction_input.txt
          Some examples are shown in ACME_codes/binding_prediction/prediction_input_example.txt
    (5) Run binding_prediction.py
          Example command: python binding_prediction.py
    (6) The results will be saved to ACME_codes/results/binding_prediction.txt
    
2. Repeating the experiments in the ACME paper.
  Please download the ACME_codes/ folder and follow the instructions in the Readme_ACME_repeat.txt file.
  
Yan Hu
School of Life Sciences, Tsinghua University 
     
