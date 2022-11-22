import glob
import pandas as pd
import argparse
from pathlib import Path
from sksurv.linear_model import CoxPHSurvivalAnalysis, CoxnetSurvivalAnalysis
from sklearn.model_selection import GridSearchCV, GroupKFold

parser = argparse.ArgumentParser(description="Just an example",
                                 formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument("feature", help="Input feature file")
parser.add_argument("output", help="output file")
args = parser.parse_args()
config = vars(args)
print("Input: "+config["feature"]+" Output: "+config["output"])


df = features[j]
df = df.dropna(axis='columns') #only necessary for coxnet as it can't handle NA values
df.PD_event = df.PD_event.astype(bool)

output = pd.DataFrame(columns = ["SampleID","Fold_assign","DFS_orig","DFS","PD_event","Cox_pred","GBM_pred","RSF_pred","SVM_pred"])
    
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
        cox_hazards = cox_model.predict(test_var)
        test_pred['Cox_pred'] = cox_hazards.tolist()
        
        gbm_model = GradientBoostingSurvivalAnalysis(
            random_state=1234
        )
        gbm_model = gbm_model.fit(train_var,train_obs)
        gbm_hazards = gbm_model.predict(test_var)
        test_pred['GBM_pred'] = gbm_hazards.tolist()

        rsf_model = RandomSurvivalForest(n_jobs = 8,
                           random_state= 1234)
        rsf_model = rsf_model.fit(train_var,train_obs)
        rsf_hazards = rsf_model.predict(test_var)
        test_pred['RSF_pred'] = rsf_hazards.tolist()

        SVM_model = FastSurvivalSVM(
            max_iter = 1000,
            random_state=1234
        )
        SVM_model = SVM_model.fit(train_var,train_obs)
        SVM_hazards = SVM_model.predict(test_var)
        
        test_pred['SVM_pred'] = SVM_hazards.tolist()

        output = pd.concat([output,test_pred])

    tmp_str = "/work/dtang/MRD_P1_P3_seperate_models/Output/KY249_N87/P1_without_baseline/" + str(feature_names[j]) + ".csv"
    filepath = Path(tmp_str)  
    filepath.parent.mkdir(parents=True, exist_ok=True)  
    output.to_csv(filepath, index=False)  



# In[ ]:




