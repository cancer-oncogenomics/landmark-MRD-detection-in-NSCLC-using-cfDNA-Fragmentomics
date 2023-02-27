import glob
import pandas as pd
import numpy as np
import sklearn
import argparse
from pathlib import Path
from sklearn.metrics import roc_auc_score
from cox_nnet3 import * 
import sys
import os
pd.options.mode.chained_assignment = None

parser = argparse.ArgumentParser(description="Just an example",
                                 formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument("feature", help="Input feature file location")
parser.add_argument("outfile", help="output file location")
args = parser.parse_args()
config = vars(args)
print("Input: "+config["feature"]+" Output: "+config["outfile"])

df = pd.read_csv(config["feature"])
df = df.dropna(axis='columns') #drop NA columns
output = pd.DataFrame(columns = ["SampleID","Fold_assign","DFS","DFS_adjusted","PD_event","Coxnnet_pred"])

old_stdout = sys.stdout
sys.stdout = open(os.devnull, "w")

for test_index in df.Fold_assign: 
    tmp_df = df.dropna(axis='columns') 
    train_var = tmp_df.loc[tmp_df['Fold_assign'] != test_index].iloc[:,9:]
    train_time = tmp_df.loc[tmp_df['Fold_assign'] != test_index].iloc[:,3]
    train_stat = tmp_df.loc[tmp_df['Fold_assign'] != test_index].iloc[:,4]

    train_time = np.array(train_time)
    train_stat = np.array(train_stat)

    test_var = tmp_df.loc[tmp_df['Fold_assign'] == test_index].iloc[:,9:]
    test_time = tmp_df.loc[tmp_df['Fold_assign'] == test_index].iloc[:,3]
    test_stat = tmp_df.loc[tmp_df['Fold_assign'] == test_index].iloc[:,4]

    test_pred = tmp_df.loc[tmp_df['Fold_assign'] == test_index].iloc[:,[0,1,2,3,4]] #to add above

    nnet_model, cost_iter= trainCoxMlp(train_var, 
                                       train_time, 
                                       train_stat, 
                                       verbose=False)
    cox_hazards = nnet_model.predictNewData(test_var)

    test_pred['Coxnnet_pred'] = cox_hazards.tolist()
    output = pd.concat([output,test_pred])

sys.stdout = old_stdout
output['PD_event'] = pd.to_numeric(output['PD_event'])

print("AUC: "+str(roc_auc_score(output.PD_event, output.Coxnnet_pred)))
output.to_csv("NNET_"+config["outfile"]+".csv", index=False)  
