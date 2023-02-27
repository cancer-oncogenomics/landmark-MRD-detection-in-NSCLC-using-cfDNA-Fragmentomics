import glob
import pandas as pd
import numpy as np
import argparse
from pathlib import Path
from sksurv.linear_model import CoxPHSurvivalAnalysis, CoxnetSurvivalAnalysis
from sklearn.model_selection import GridSearchCV, GroupKFold
import sksurv
print(sksurv.__version__)


parser = argparse.ArgumentParser(description="Just an example",
                                 formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument("-b", "--baseline", action="store_true", help="adding baseline to features")
parser.add_argument("feature", help="Input feature file location")
parser.add_argument("outfile", help="output file location")
args = parser.parse_args()
config = vars(args)
print("Input: "+config["feature"]+" Output: "+config["outfile"])

df = pd.read_csv("P1_feature.csv") #### needs to be deleted
df = pd.read_csv(config["feature"])
df = df.dropna(axis='columns') #drop NA values
df.PD_event = df.PD_event.astype(bool)

output = pd.DataFrame(columns = ["SampleID","Fold_assign","DFS","DFS_adjusted","PD_event","Cox_pred"])

test_index = 1 #### needs to be deleted

for test_index in df.Fold_assign: 
    tmp_df = df.dropna(axis='columns') 
    train_var = tmp_df.loc[tmp_df['Fold_assign'] != test_index].iloc[:,9:]
    train_obs = tmp_df.loc[tmp_df['Fold_assign'] != test_index].iloc[:,[4,3]]
    train_obs = train_obs.to_records(index=False)
    test_var = tmp_df.loc[tmp_df['Fold_assign'] == test_index].iloc[:,9:]
    test_obs = tmp_df.loc[tmp_df['Fold_assign'] == test_index].iloc[:,[4,3]]
    test_pred = tmp_df.loc[tmp_df['Fold_assign'] == test_index].iloc[:,[0,1,2,3,4]] #to add above

    cox_model = CoxnetSurvivalAnalysis(
        l1_ratio=1, 
        alpha_min_ratio= 0.01, 
    )
    cox_model = cox_model.fit(train_var,train_obs)
    cox_lp = cox_model.predict(test_var)
    test_pred['Cox_pred'] = cox_lp.tolist()
    output = pd.concat([output,test_pred])

    coefficients_lasso = pd.DataFrame(
        cox_model.coef_,
        index=train_var.columns,
        columns=np.round(cox_model.alphas_, 5)
        )
    non_zero = np.sum(coefficients_lasso.iloc[:, 0] != 0)





output.to_csv(config["outfile"], index=False)  
