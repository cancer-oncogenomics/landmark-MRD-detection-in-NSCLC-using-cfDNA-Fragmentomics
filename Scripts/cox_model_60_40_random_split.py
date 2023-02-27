import glob
import random
import pandas as pd
import numpy as np
import argparse
from pathlib import Path
from sksurv.linear_model import CoxPHSurvivalAnalysis, CoxnetSurvivalAnalysis
from sklearn.model_selection import GridSearchCV, GroupKFold, train_test_split
from sklearn.metrics import roc_auc_score
pd.options.mode.chained_assignment = None

parser = argparse.ArgumentParser(description="Just an example",
                                 formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument("seeds", help="Input feature file location")
parser.add_argument("repeats", help="Input feature file location")
parser.add_argument("feature", help="Input feature file location")
parser.add_argument("outfile", help="output file location")
args = parser.parse_args()
config = vars(args)
print("Input: "+config["feature"]+" Output: "+config["outfile"]+" Seed: "+config["seeds"]+" Reps: "+config["repeats"])

config["repeats"] = int(config["repeats"])
config["seeds"] = int(config["seeds"])

random.seed(config["seeds"])
random_list = random.sample(range(1, 2**20), config["repeats"])

df_orig = pd.read_csv(config["feature"])

df = df_orig
df = df.dropna(axis='columns') #coxnet can't handle NA values
df.PD_event = df.PD_event.astype(bool)
    
output = pd.DataFrame(columns = ["Run","Seed","AUC"])

for runs in range(0,config["repeats"]):
    Train, Test = train_test_split(df, test_size=0.4, random_state=random_list[runs])
    Train_var = Train.iloc[:,9:]
    Train_obs = Train.iloc[:,[4,3]]
    Train_obs = Train_obs.to_records(index=False)

    Test_var = Test.iloc[:,9:]
    Test_pred = Test.iloc[:,[0,1,2,3,4]]

    cox_model = CoxnetSurvivalAnalysis(
        l1_ratio=1, 
        alpha_min_ratio= 0.01, 
    )
    cox_model = cox_model.fit(Train_var,Train_obs)
    cox_hazards = cox_model.predict(Test_var)
    Test_pred['Cox_pred'] = cox_hazards.tolist()
    Test_pred.PD_event = Test_pred.PD_event.astype(int)

    Test_pred.to_csv("60_40_random_split_output/Run_"+str(runs)+"_"+config["outfile"]+"_score.csv",index =False)
    
    tmp_out = pd.DataFrame([[runs,random_list[runs],roc_auc_score(Test_pred.PD_event, Test_pred.Cox_pred)]])
    tmp_out.columns = ['Run',"Seed","AUC"]
    output = pd.concat([output,tmp_out])

output.to_csv(config["outfile"]+"_AUC.txt", index = False)