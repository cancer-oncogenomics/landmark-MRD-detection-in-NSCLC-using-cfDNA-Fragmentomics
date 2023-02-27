import glob
import pandas as pd
import numpy as np
import argparse
from pathlib import Path
from sksurv.linear_model import CoxPHSurvivalAnalysis, CoxnetSurvivalAnalysis
from sklearn.model_selection import RepeatedKFold, KFold, cross_val_predict, cross_val_score
from sklearn.metrics import roc_auc_score
import random

parser = argparse.ArgumentParser(description="Just an example",
                                 formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument("Kfold", help="K value for K-fold cross-validation")
parser.add_argument("Repeat", help="N for repeated K-fold cross-validation")
parser.add_argument("Seeds", help="Random seeds")
parser.add_argument("feature", help="Input feature file location")
parser.add_argument("outfile", help="output file location")
args = parser.parse_args()
config = vars(args)
print("Input: "+config["feature"]+" Output: "+config["outfile"]+" Kfold: "+config["Kfold"]+" Repeats: "+config["Repeat"]+" Seed:"+config["Seeds"])

config["Repeat"] = int(config["Repeat"])
config["Kfold"] = int(config["Kfold"])
config["Seeds"] = int(config["Seeds"])

random.seed(config["Seeds"])
random_list = random.sample(range(1, 2**20), config["Repeat"])

df = pd.read_csv(config["feature"])
df = df.dropna(axis='columns') #coxnet can't handle NA values
output_pred = df.iloc[:,0:8]
output_auc = pd.DataFrame(columns = ["Run","Seed","AUC"])

df.PD_event = df.PD_event.astype(bool)
df_obs = df.iloc[:,[4,3]]
df_obs = df_obs.to_records(index=False)
df_var = df.iloc[:,9:]

cox_model = CoxnetSurvivalAnalysis(
    l1_ratio=1, 
    alpha_min_ratio= 0.01, 
)

for run_number in range(0,config["Repeat"]):
    cv = KFold(n_splits= config["Kfold"], shuffle = True, random_state= random_list[run_number])
    scores = cross_val_predict(cox_model, df_var, df_obs, cv=cv)
    output_pred["Run_"+str(run_number)+"_pred"] = scores.tolist()

    tmp_out = pd.DataFrame([[run_number,random_list[run_number],roc_auc_score(output_pred.PD_event, scores)]])
    tmp_out.columns = ['Run',"Seed","AUC"]
    output_auc = pd.concat([output_auc,tmp_out])

tmp_col = output_pred.iloc[:,8:]
output_pred['Pred_average'] = tmp_col.mean(axis = 1)
output_pred.to_csv(config["outfile"]+"_"+str(config["Repeat"])+"_Repeat_"+str(config["Kfold"])+"_fold_score.txt", index = False)

average_out = pd.DataFrame([["Average","NA",roc_auc_score(output_pred.PD_event, output_pred.Pred_average)]])
average_out.columns = ['Run',"Seed","AUC"]
output_auc = pd.concat([output_auc,average_out])

output_auc.to_csv(config["outfile"]+"_"+str(config["Repeat"])+"_Repeat_"+str(config["Kfold"])+"_fold_AUC.txt", index = False)