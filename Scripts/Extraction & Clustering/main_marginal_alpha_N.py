from __future__ import division
from tkinter import filedialog
import matplotlib.pyplot as plt

import warnings
warnings.filterwarnings("ignore")

from sklearn.mixture import BayesianGaussianMixture, GaussianMixture
from sklearn.preprocessing import StandardScaler
import pandas as pd
import pickle
from collections import OrderedDict
from utils.helpers import *

clustering_method = 'DPGMM'
# clustering_method = 'BICGMM'
just_means = True
COVARIANCE_TYPE = 'full'
RANDOM_RESTARTS = 200


column_width = 3.3
auto_save = True

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

if clustering_method == 'DPGMM':
    save_summary_name = participant_directory + 'sitting_features_summary_results_' + s_id + '.xlsx'
    model_file = participant_directory + 'model_' + s_id + '.pk'

print ('TOTAL USED:')
print (len(feature_matrix))


avg_ll = np.zeros(shape=RANDOM_RESTARTS)
if clustering_method == 'DPGMM':
    for r in range(RANDOM_RESTARTS):
        N_COMPONENTS = len(feature_matrix)
        model = BayesianGaussianMixture(n_components=N_COMPONENTS,max_iter=300,
                                        verbose=False,random_state=r,
                                        covariance_type=COVARIANCE_TYPE).fit(feature_matrix)
        avg_ll[r] = np.sum(model.score(feature_matrix))
        print (avg_ll[r])
    N = [n for n in range(1, N_COMPONENTS, 5)]
    N.append(N_COMPONENTS)
    elbo = []
    n_clusters = []
    n_clusters_below = []
    for n in N:
        model = BayesianGaussianMixture(n_components=n,max_iter=300,
                                        verbose=True,random_state=np.argmax(avg_ll),
                                        covariance_type=COVARIANCE_TYPE,
                                        weight_concentration_prior=1 / N_COMPONENTS).fit(feature_matrix)
        elbo.append(model.lower_bound_)
        weights = model.weights_
        predict_prob = model.predict_proba(feature_matrix)
        labels = model.predict(feature_matrix)
        clusters = np.unique(labels)
        cluster_weights = weights[clusters]
        n_clusters.append(np.sum(cluster_weights >= 0.05))
        n_clusters_below.append(np.sum(cluster_weights < 0.05))

fig,ax = plt.subplots()
plt.title("Marginal likelihood for prior N")
plt.xlabel("prior N_components")
plt.ylabel("marginal likelihood")
plt.plot(elbo,'-x')
ax.set_xticks(range(len(N)))
ax.set_xticklabels(N, fontdict={'fontsize': 4})
# ax.set_xticklabels(alpha_list)
for i in range(len(n_clusters)):
    plt.annotate(str(n_clusters[i]),(i,elbo[i] + 20),fontsize=6)
    plt.annotate(str(n_clusters_below[i]),(i,elbo[i]), (i,elbo[i] - 20),fontsize=6)
if auto_save:
    # fig.set_size_inches(column_width,column_width)
    plt.savefig('marginal_N_' + str(s_id) + '.png',dpi=300)
else:
    plt.show()
print (np.sum(model.score(feature_matrix)))


alpha = [0.00001, 0.0001, 0.001, 0.01, .1, 1, 10, 100, 1000, 10000]
alpha_list = ['10^-5', '10^-4', '10^-3', '10^-2', '10^-1', '1', '10', '10^2', '10^3', '10^4']
elbo = []
n_clusters = []
n_clusters_below = []
for a in alpha:
    model = BayesianGaussianMixture(n_components=N_COMPONENTS, max_iter=300,
                                    verbose=True,random_state=np.argmax(avg_ll),
                                    covariance_type=COVARIANCE_TYPE,
                                    weight_concentration_prior=a).fit(feature_matrix)
    elbo.append(model.lower_bound_)
    weights = model.weights_
    predict_prob = model.predict_proba(feature_matrix)
    labels = model.predict(feature_matrix)
    clusters = np.unique(labels)
    cluster_weights = weights[clusters]
    n_clusters.append(np.sum(cluster_weights >= 0.05))
    n_clusters_below.append(np.sum(cluster_weights < 0.05))

fig,ax = plt.subplots()
plt.title("Marginal likelihood for changing alpha")
plt.xlabel("alpha value")
plt.ylabel("marginal likelihood")
plt.plot(elbo,'-x')
ax.set_xticks(range(len(alpha)))
ax.set_xticklabels(alpha_list)
for i in range(len(n_clusters)):
    plt.annotate(str(n_clusters[i]),(i,elbo[i] + 80),fontsize=6)
    plt.annotate(str(n_clusters_below[i]),(i,elbo[i]), (i,elbo[i] - 80),fontsize=6)
if auto_save:
    # fig.set_size_inches(column_width,column_width)
    plt.savefig('marginal_a_' + str(s_id) + '.png',dpi=300)
else:
    plt.show()