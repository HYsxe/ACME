ACME: Attention-basd Convlutional neural network for MHC Epitope prediction

1. Description

    The codes in this folder can be used to reproduce our work.

2. Requirements

    Python2, Linux system

    Python packages: re, keras, sklearn, math, scipy, numpy

3. Installation

    Before running the codes, you need to change the 'main_dir' in main.py
    
    to the path of this folder. For example
    
        main_dir = "/home/user_1/ACME-master/ACME_codes/"
    
    Then unzip the proteome data file using the following command.
    
        unzip proteome.zip

    You can invoke different functions using the main() function. Just change the argument
    
    passed to the main() function in main.py. For example, if you want to use the first function,
    
    which is to train new models and save them to a folder, you can change the last line of 
    
    main.py to main(0) where the '0' stands for the first function. ('i' stands for the i-1th function)
    
    Next, save the file and type the following command in the terminal.
    
        python main.py
	
    Then new models will be trained and saved.

4. Input and output

    (1)The IEDB MHC class I binding dataset used in 5-fold cross validation is in the /binding_data folder.

    (2)The IEDB weekly benchmark datasets are in the /IEDB_benchmarking_datasets folder.

    (3)Newly trained models will be saved to the /models folder.

    (4)The output of all the other functions will be saved to the /results folder.

5. Functions

    To run the k th function, change the last line of main.py to main(i). 
    
    You can directly run the k th fuction. (For example if you want to run the third function, just run main(2) without running main(0) and main(1)).

    0: Train new models and save them to /model.
    
	If you download the ACME_codes/ folder, we have some trained model in the /model folder that you can use directly,
	
	so normally you can just use these models and skip running main(0). If you want to train new models, 
	
	you can run main(0) which will train new models to replace the existing ones.

    1: 5-fold cross-validation. 
    
	To choose peptides of length L for testing, open
	
	main_cross_validation.py and change test_len in line 24 to L. 
	
	For example, when testing on 10-mers. Change line 24 to
	
	training_data, test_dicts = preparing_data(data_dict, n_splits, test_len = 10)
	
	The default value of test_len is 9.
	
	The output file is /results/cross_validation.txt.

    2: Testing the trained models on external datasets.

	The output file is /results/testing_on_IEDB_benchmark_datasets.txt	

    3: Test whether or not the attention module can detect important residues.

	For details, see section 1.9 in Supplementary Notes. 
	
	The output file is /results/attention_test.txt
	
	The output contains the average D1 and D2 values for each allele and the Kruskal test results.	
	
	It also contains the average D1, D2 and test results for all alleles (last line).

    4: Generate binding motifs (using peptides with high binding affinities).
    
	The output is a dictionary. The keys are alleles and their corresponding
	
	values are the attenttion matrix A described Supplementary Notes Section 1.9.
	
	The output file is /results/binder_motif.txt

    5: Generate binding motifs for non-binders (i.e. using peptides with low binding affinities).
    
	The output file is /results/non_binder_motif.txt

    6: Make predictions for alleles with no training data.
    
	See Section 3.2 for details. You can specify the alleles to be tested in ACME_codes/main_leave_one_out.py
	
	For example, to test A*02:01, change line 28 to: 
	
	    alleles = ["A*02:01"]
	
	The output file is /results/ACME_leave_one_out.txt

    7: Test whether the anti-anchor residues can impair binding. 
    
	The output file is /results/anti_anchor.txt	
	
	The output contains the allele, the amino acid type of the anti-anchor residue its position in the peptide,
	
	averaged prediction score after replacing peptide residues with the anti-anchor,
	
	averaged prediction score after replacing the same residues with a random residue
	
	and the Kruskal test results.

    8: Cross validation using a model without the attention module.
    
	The output file is results/cross_validation_without_attention_9mer.txt
	
	You can also run the test on peptides with other lengths. See the instructions for function 1.
	
    9: Cross validation using a model without the convolutional module.
    
        The output file is results/cross_validation_without_CNN_9mer.txt
	
	You can also run the test on peptides with other lengths. See the instructions for function 1.
	
    10: Calculating the reference information (RI) of different MHC alleles.
    
        Calculate the reference information of each allele when their data are not present in the training set.
	
	The output file is results/MHC_RI.txt 
	
	For the definition of Reference information, see Supplementary Notes Section 1.6. 

6. Problems

If you have problems using ACME, please contact yanhu@g.harvard.edu.

Hope you enjoy using the software!

Yan Hu, 
School of Life Sciences, Tsinghua University
Institute for Interdisciplinary Information Sciences, Tsinghua University
Oct, 2018




