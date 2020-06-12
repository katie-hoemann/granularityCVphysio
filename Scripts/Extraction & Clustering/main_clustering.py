from __future__ import division
from tkinter import filedialog

import warnings
warnings.filterwarnings("ignore")

from sklearn.mixture import BayesianGaussianMixture
from sklearn.preprocessing import StandardScaler
import pandas as pd
import pickle
from collections import OrderedDict
from utils.helpers import *

clustering_method = 'DPGMM'
just_means = True
THRESHOLD = 0.7
MAX_CHANGE_THRESHOLD = 400
MAX_SV_THRESHOLD = 200
COVARIANCE_TYPE = 'full'
RANDOM_RESTARTS = 200



feature_list = ['RSA Mean', 'IBI Mean', 'RR Mean', 'PEP Mean']
def get_immediate_subdirectories(a_dir):
    return [name for name in os.listdir(a_dir)
            if os.path.isdir(os.path.join(a_dir, name))]

print ('Please browse to the subject folder:')

data_directory = filedialog.askdirectory(initialdir =
                                           r'C:\Users\zulqa\Documents\PhD Lab\ARI_Data\Sitting_Segments')
participant_directory = data_directory + '/'


s_id = data_directory[-3:]
original_summary_file = participant_directory + 'sitting_features_summary_' + s_id + '.xlsx'
original_feature_summary_pd = pd.read_excel(original_summary_file)
feature_summary_pd = original_feature_summary_pd[original_feature_summary_pd['Bad Flag'] == 0]
feature_matrix = feature_summary_pd[feature_list].values

original_feature_matrix = feature_matrix.copy()
standard_transform =StandardScaler().fit(original_feature_matrix)

feature_matrix = standard_transform.transform(feature_matrix)

save_summary_name = participant_directory + 'sitting_features_summary_results_' + s_id + '.xlsx'
model_file = participant_directory + 'model_' + s_id + '.pk'


print ('TOTAL USED:')
print (len(feature_matrix))


avg_ll = np.zeros(shape=RANDOM_RESTARTS)
if clustering_method == 'DPGMM':
    for r in range(RANDOM_RESTARTS):
        N_COMPONENTS = len(feature_matrix)
        model = BayesianGaussianMixture(n_components=N_COMPONENTS,max_iter=300,verbose=False,random_state=r,covariance_type=COVARIANCE_TYPE).fit(feature_matrix)
        avg_ll[r] = np.sum(model.score(feature_matrix))
        print (avg_ll[r])
    model = BayesianGaussianMixture(n_components=N_COMPONENTS,max_iter=300,verbose=True,random_state=np.argmax(avg_ll),covariance_type=COVARIANCE_TYPE).fit(feature_matrix)
    weights = model.weights_
    predict_prob = model.predict_proba(feature_matrix)
    labels = model.predict(feature_matrix)
    clusters = np.unique(labels)
    cluster_weights = weights[clusters]
    print (np.sum(model.score(feature_matrix)))

model_file = open(model_file, 'wb')
pickle.dump(model, model_file)
model_file.close()

for c in clusters:
    visualize = original_feature_matrix[np.where(labels == c)[0],:]
    print (len(np.where(labels == c)[0]))

posterior_prob = predict_prob[:,clusters]
cluster_means = model.means_
cluster_means = cluster_means[clusters,:]

new_columns = ['C_ID']
for i in range(len(clusters)):
    new_columns.append('p_'+str(i))

new_labels = np.array([int(np.where(clusters == l)[0][0]) for l in labels])
save_pd = original_feature_summary_pd.copy()
save_pd = save_pd.reindex(save_pd.columns.tolist() + new_columns, axis=1)
use_index = save_pd.index[original_feature_summary_pd['Bad Flag'] == 0]
save_pd.loc[use_index,'C_ID'] = new_labels

for i in range(len(clusters)):
    save_pd.loc[use_index, 'p_'+str(i)] = posterior_prob[:,i]


transformed_cluster_means = standard_transform.inverse_transform(cluster_means)
column_names = feature_list + ['C_ID','Cluster_Weights']
dictionary = OrderedDict([(column_names[0],transformed_cluster_means[:,0]),
                          (column_names[1],transformed_cluster_means[:,1]),
                          (column_names[2],transformed_cluster_means[:,2]),
                          (column_names[3],transformed_cluster_means[:,3]),
                          (column_names[4],np.arange(len(clusters))),
                          (column_names[5],cluster_weights),])
writer = pd.ExcelWriter(save_summary_name)
save_pd.to_excel(writer, 'Sheet1')
pd.DataFrame(dictionary).to_excel(writer,'Sheet2')
writer.save()


