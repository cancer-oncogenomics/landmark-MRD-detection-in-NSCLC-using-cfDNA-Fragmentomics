# landmark-MRD-detection-in-NSCLC-using-cfDNA-Fragmentomics-

For model generation: cox_model_construction.py
usage: cox_model_construction.py [-h] [-b] feature outfile

For Figure/Table generation: Figure_generation.R

60 40 split

python3 cox_model_60_40_random_split.py 12345 50 P1_feature.csv P1

python3 cox_model_60_40_random_split.py 12345 50 P2_feature.csv P2

K-fold CV

python3 cox_model_K_fold_cv.py 10 10 12345 P1_feature.csv P1

python3 cox_model_K_fold_cv.py 10 10 12345 P2_feature.csv P2

python3 cox_model_K_fold_cv.py 5 5 12345 P1_feature.csv P1

python3 cox_model_K_fold_cv.py 5 5 12345 P2_feature.csv P2


Coxnnet LOOCV

python3 nnet_model_validation.py P1_feature.csv P1

AUC: 0.6752717391304348

python3 nnet_model_validation.py P2_feature.csv P2

AUC: 0.4506001846722068
