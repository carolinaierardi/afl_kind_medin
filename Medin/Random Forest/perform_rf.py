import os
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import train_test_split
from sklearn.metrics import accuracy_score, roc_auc_score, classification_report
import numpy as np
import pandas as pd

import matplotlib.pyplot as plt
import seaborn as sns

from sklearn.decomposition import PCA
from sklearn.manifold import TSNE
import umap

plt.rcParams['font.family'] = 'Arial'
plt.rcParams['font.family'] = 'Arial'
plt.rcParams['xtick.labelsize'] = 12
plt.rcParams['ytick.labelsize'] = 12
plt.rcParams['axes.titlesize'] = 20
plt.rcParams['axes.labelsize'] = 15

def calculate_metrics(true_values, predicted_values, y_pred_proba,average='weighted'):
    metrics = {}
    metrics['classification_report'] = classification_report(true_values, predicted_values, output_dict=True)
    metrics['accuracy'] = accuracy_score(true_values, predicted_values)    
    if y_pred_proba is not None:
        if len(y_pred_proba.shape) == 1 or y_pred_proba.shape[1] == 2:
            metrics['roc_auc'] = roc_auc_score(true_values, y_pred_proba if len(y_pred_proba.shape) == 1 else y_pred_proba[:, 1])
        else:
            metrics['roc_auc'] = roc_auc_score(true_values, y_pred_proba, multi_class='ovr', average=average)
    else:
        metrics['roc_auc'] = 'N/A (y_pred_proba not provided)'
    return metrics

def perform_dim_red(fgp, aff, title, png_path, scaling_factor = 2):

    """
     Parameters
    ----------
    fgp : np.array
        Fingerprints; each row represents a molecule.
    aff : np.array
        Radius to use as parameters.
    title: str
        Title of graph
    png_path: str (.png)
        Path to store .png final image
    scaling_factor : int, optional. 
        Scaling factor for labels and titles. The default is 2.

    """

    pca = PCA(n_components=2)
    pca_result = pca.fit_transform(fgp)
    print("PCA done")

    # t-SNE
    tsne = TSNE(n_components=2, random_state=42)
    tsne_result = tsne.fit_transform(fgp)
    print("t-SNE done")

    # UMAP
    umap_model = umap.UMAP(n_components=2, random_state=42)
    umap_result = umap_model.fit_transform(fgp)
    print("UMAP done")


    cm = plt.get_cmap('rainbow')
    fig, axes = plt.subplots(nrows=1, ncols=3, figsize=(15*scaling_factor,5*scaling_factor))

    pc = axes[0].scatter(pca_result[:, 0], pca_result[:, 1], s=15, c=aff, cmap = cm, alpha=0.7)
    axes[0].set_title("PCA", size = 20*scaling_factor)
    axes[0].set_xlabel('Component 1', size = 15*scaling_factor)
    axes[0].set_ylabel('Component 2', size = 15*scaling_factor)
    cbar = fig.colorbar(pc, ax = axes[0])
    cbar.ax.tick_params(labelsize=10*scaling_factor)

    tn = axes[1].scatter(tsne_result[:, 0], tsne_result[:, 1], s=15,c=aff, cmap = cm,alpha=0.7)
    axes[1].set_title("t-SNE", size = 20*scaling_factor)
    axes[1].set_xlabel('Component 1',  size = 15*scaling_factor)
    axes[1].set_ylabel('Component 2', size = 15*scaling_factor)
    cbar = fig.colorbar(tn, ax = axes[1])
    cbar.ax.tick_params(labelsize=10*scaling_factor)

    um = axes[2].scatter(umap_result[:, 0], umap_result[:, 1], s=15, c=aff, cmap = cm,alpha=0.7)
    axes[2].set_title("UMAP", size = 20*scaling_factor)
    axes[2].set_xlabel('Component 1',  size = 15*scaling_factor)
    axes[2].set_ylabel('Component 2', size = 15*scaling_factor)
    cbar = fig.colorbar(um, ax = axes[2])
    cbar.ax.tick_params(labelsize=10*scaling_factor)


    fig.suptitle(title, size = 25*scaling_factor)
    fig.tight_layout(h_pad = 2)                                                          #tight layout so there is no overlay between plots
    plt.savefig(png_path)


#Import data
os.chdir("/home/cm2231/rds/project/rds-fYhPa3It0hY/carolina/inferences/analysis/figures/figures/hype")
feat = pd.read_csv("/home/cm2231/rds/project/rds-mphil-rkVIP0ZCj0k/carolina/test_Medin/PaDEL_class.csv")
aff = pd.read_csv("/home/cm2231/rds/project/rds-fYhPa3It0hY/carolina/inferences/data/processed/hype/rf_sample/top_bottom_sample.csv")

merged_df = pd.merge(feat, aff, how = "left", left_on = ["Name"], right_on = ["id"])
clean_df = merged_df.dropna()

train, test= train_test_split(clean_df, test_size=0.1, shuffle = True, random_state=42)

clf = RandomForestClassifier(max_depth=2, random_state=0)
res = clf.fit(train[train.columns[1:-4]], train["affinity class"])

#See how well the data predicts affinity
test_prediction = clf.predict(test[test.columns[1:-4]])
test_proba = clf.predict_proba(test[test.columns[1:-4]])

met = calculate_metrics(test["affinity class"], test_prediction, test_proba)
print(met)
# met = pd.DataFrame.from_dict(data = met, orient = "index").T
# met.to_csv("rf_class_prediction.csv", index = False)

feature_names = clean_df.columns[1:-4]
feature_names = feature_names[res.feature_importances_ != 0]
nonz = res.feature_importances_[res.feature_importances_ != 0]

df_features = pd.DataFrame({"features":feature_names, "importance": nonz})
df_features.sort_values(by = "importance", ascending = False, inplace = True)

selected_feat = clean_df[df_features["features"].head(15).to_list()]
selected_feat["affinity"] = clean_df[" pred_affinity"]
selected_feat["id"] = clean_df["id"]
df_features.to_csv("best_class_features_importance.csv", index = False)
selected_feat.to_csv("best_Padel_features_class.csv", index = False)

perform_dim_red(selected_feat[selected_feat.columns[:-2]], clean_df[" pred_affinity"], "Molecular descriptors", "moldescrip_class.png")

