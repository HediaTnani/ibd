import json
import cobra
import os
import numpy as np
import pandas as pd


import matplotlib.pyplot as plt
import seaborn as sns
from seaborn import clustermap

from sklearn.ensemble import RandomForestClassifier
from sklearn.ensemble import RandomForestRegressor
from sklearn import metrics

# Matplotlib defaults
import matplotlib.pyplot as plt
import matplotlib
import matplotlib.ticker as mticker
matplotlib.rcParams["savefig.dpi"] = 300
matplotlib.rcParams['figure.dpi'] = 300
#Rivanna hasn't updated fonts yet to include sans-serif
# matplotlib.rcParams['font.sans-serif'] = "Arial"
# matplotlib.rcParams['font.family'] = "sans-serif"

SMALLER_SIZE = 8
SMALL_SIZE = 14
MEDIUM_SIZE = 16
BIGGER_SIZE = 16

matplotlib.rc('font', size=SMALL_SIZE)          # controls default text sizes
matplotlib.rc('axes', titlesize=MEDIUM_SIZE)     # fontsize of the axes title
matplotlib.rc('axes', labelsize=SMALL_SIZE)    # fontsize of the x and y labels
matplotlib.rc('xtick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
matplotlib.rc('ytick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
matplotlib.rc('legend', fontsize=SMALLER_SIZE)    # legend fontsize
matplotlib.rc('figure', titlesize=BIGGER_SIZE)  # fontsize of the figure title

# read in the base model
# load the metabolic model. We'll use iHsa with liver-specific function removed.
ihsa = cobra.io.load_json_model('../data/ihsa/ihsa_liver_removed_open_exchanges.json')

# read in the kallisto output, which we've converted the IDs for using gene_symbols_to_EntrezID.py
data = pd.read_csv(
    '../data/transcriptomics/13464__coding_52EED_25Ctl_17Celiac_ENTREZ.txt',
    index_col=0, sep = '\t')
# make sure the index is string, not int, for integration with GIMME
new_index = ['gene'+str(i) for i in data.index.tolist()]
data.index = new_index

# read in the metadata/patient IDs
metadata = pd.read_csv('../data/transcriptomics/52EED_25Ctl_17celiac_metadata_sendout_Sana_Aug6_2019.txt',
                      index_col=0, sep = '\t')

# rename the erroneously-labelled 'gene' column to 'patient_rnaseq_id'
metadata.index.name = 'patient_rnaseq_id'

# read in the anthropometric + other metada
anthro = pd.read_excel('../data/yael_anthro_20190801 no names_LE.xlsx')
anthro = anthro.rename({'studyid_':'PID'},axis=1)
control_celiac_metadata = pd.read_excel('../data/WUPAX_CCHMC_metadata_23Sep2019_LM.xlsx')

# add sex for AKU biopsies
all_sex = pd.concat([anthro[['PID','gender']].rename({'gender':'Gender'},axis=1),control_celiac_metadata[['PID','Gender']]])
metadata = metadata.reset_index().merge(all_sex,on='PID',how='left').set_index('patient_rnaseq_id')

# Load the reaction fluxes from the gimme solutions
gimme_fluxes = pd.read_csv('../results/compiled_gimme_solutions.tsv', sep='\t',index_col=0)

# remove the celiac samples
gimme_fluxes = gimme_fluxes.loc[metadata.index[(metadata['Diagnosis'] != 'CELIAC')]]
#remove all 0 columns
gimme_fluxes = gimme_fluxes[gimme_fluxes.columns[(abs(gimme_fluxes).sum() > 1E-15)]]


# repeat using the internal reactions only
gimme_internal_diagnosis = RandomForestClassifier(n_estimators = 50000,
                                                  oob_score = True, 
                                                  class_weight='balanced_subsample')

x = gimme_fluxes.copy()

#Add sex
x['Gender'] = metadata.loc[x.index,'Gender'].values == 'Female'

y = metadata.loc[x.index,'Diagnosis'].values
feature_names = x.columns

fit_gimme = gimme_internal_diagnosis.fit(x,y)

from sklearn.metrics import confusion_matrix
y_pred = [a<b for a,b in fit_gimme.oob_decision_function_]
y_true = [b == 'FTT_enteropathy' for b in y]
print(confusion_matrix(y_true, y_pred)) #TODO: convert to more readily-saved output

# ROC curve and precision-recall curve
from sklearn import metrics
# scores are the probability assigned to the positive class, which we defined as EED
scores = fit_gimme.oob_decision_function_[:,1]
fpr, tpr, ss_thresholds = metrics.roc_curve(y_true, scores, pos_label=True)
precision, recall, pr_thresholds = metrics.precision_recall_curve(y_true, scores, pos_label=True)

# plot the ROC curve
fig,ax = plt.subplots()
ax.plot(fpr,tpr,lw=2, color = 'darkorange', label='ROC curve (area = %0.2f)' % metrics.roc_auc_score(y_true,scores))
ax.plot([0,1],[0,1], color = 'navy', linestyle='--')
ax.set_xlabel('False Positive Rate (1 - Sensitivity)')
ax.set_ylabel('True Positive Rate (Specificity)')
ax.legend()
ax.set_title('Classifier receiver operating characteristic\nusing internal reaction fluxes')

plt.savefig('../results/classifier_performance/all_rxns_roc.svg')
plt.savefig('../results/classifier_performance/all_rxns_roc.png')

# Plot the precision-recall curve
from inspect import signature
average_precision = metrics.average_precision_score(y_true, scores)
fig,ax = plt.subplots()
step_kwargs = ({'step': 'post'}
               if 'step' in signature(plt.fill_between).parameters
               else {})
ax.step(recall, precision, color='b', alpha=0.2,
         where='post')
ax.fill_between(recall, precision, alpha=0.2, color='b', **step_kwargs)

ax.set_xlabel('Recall')
ax.set_ylabel('Precision')
ax.set_ylim([0.0, 1.05])
ax.set_xlim([0.0, 1.0])
ax.set_title('Classifier precision-recall using internal\nreaction fluxes (Avg. precision={0:0.2f})'.format(
          average_precision))
plt.savefig('../results/classifier_performance/all_rxns_precision_recall.svg')
plt.savefig('../results/classifier_performance/all_rxns_precision_recall.png')

# extract feature importances
# TODO; add saving step
print(fit_gimme.oob_score_)
sorted_importance = pd.DataFrame(fit_gimme.feature_importances_,x.columns, columns = ['importance']).sort_values(by='importance', ascending = False)
top10_imp = sorted_importance.head(30)
ordered_rxns = []
for rxn in top10_imp.index:
    rxn_obj = ihsa.reactions.get_by_id(rxn)
    print(rxn_obj.name)
    ordered_rxns.append(rxn_obj.name)
    print(rxn_obj)
    print(rxn_obj.build_reaction_string(use_metabolite_names=True))
    print()
    
# For the top 10 important features, plot the gimme fluxes
rxns = top10_imp.index
scatter_data = x[top10_imp.index]
# normalize the fluxes by column max
scatter_data = scatter_data.divide(abs(scatter_data).mean())

scatter_data['Diagnosis'] = metadata.loc[scatter_data.index]['Diagnosis']
scatter_data.loc[scatter_data['Diagnosis'] == 'FTT_enteropathy','Diagnosis'] = "Enteropathy"
scatter_data.loc[scatter_data['Diagnosis'] == 'CONTROL','Diagnosis'] = "Control"
scatter_data['sample'] = scatter_data.index
scatter_data = pd.melt(scatter_data,id_vars=['sample','Diagnosis'],value_vars=rxns, var_name = 'reaction', value_name = 'flux')

f, ax = plt.subplots(figsize=(30, 20))
sns_plot = sns.violinplot(y='reaction',x='flux',hue='Diagnosis',data=scatter_data.iloc[::-1], order = rxns,orient='h')
sns_plot.set_yticklabels(labels = ordered_rxns,size=20, rotation=0, ha='right')
sns_plot.legend(prop={'size': 18})
sns_plot.axvline(color='grey')
for i in range(0,len(ordered_rxns)):
    sns_plot.axhline(i+0.5, color='grey', alpha = 0.5)
sns_plot.set_xlabel('Reaction flux',size=28)
sns_plot.set_ylabel('Reactions',size=28)
plt.tight_layout()
plt.savefig('../results/gimme_internal_top10.svg')
plt.savefig('../results/gimme_internal_top10.png')