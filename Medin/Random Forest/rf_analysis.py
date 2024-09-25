
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import os
import numpy as np

plt.rcParams['font.family'] = 'Arial'
plt.rcParams['font.family'] = 'Arial'
plt.rcParams['xtick.labelsize'] = 12
plt.rcParams['ytick.labelsize'] = 12
plt.rcParams['axes.titlesize'] = 20
plt.rcParams['axes.labelsize'] = 15


regression = pd.read_csv("/home/cm2231/rds/project/rds-fYhPa3It0hY/carolina/inferences/analysis/figures/figures/hype/best_Padel_features.csv")
regression_importance = pd.read_csv("/home/cm2231/rds/project/rds-fYhPa3It0hY/carolina/inferences/analysis/figures/figures/hype/best_features_importance.csv")
classification = pd.read_csv("/home/cm2231/rds/project/rds-fYhPa3It0hY/carolina/inferences/analysis/figures/figures/hype/best_Padel_features_class.csv")
classification_importance = pd.read_csv("/home/cm2231/rds/project/rds-fYhPa3It0hY/carolina/inferences/analysis/figures/figures/hype/best_class_features_importance.csv")

# List of variables to plot
os.chdir("/home/cm2231/rds/project/rds-fYhPa3It0hY/carolina/inferences/analysis/figures/figures/hype")
regression[['SpDiam_D', 'ECCEN', 'SpAD_D', 'EE_D', 'SpMax_D']] = regression[['SpDiam_D', 'ECCEN', 'SpAD_D', 'EE_D', 'SpMax_D']].apply(lambda x: np.log10(x))

# Create scatter plots
for i, var in enumerate(regression.columns[:5]):
    plt.figure(figsize = (10,6))
    ax = sns.jointplot(data=regression, x=var, y="affinity", kind = "reg", 
              color = '#7BF1A8')
    ax.fig.suptitle(f"r = {round(np.corrcoef(regression[var], regression['affinity'])[1,0],3)}", size=16)
    ax.set_axis_labels(f"log10("+var+")", 'affinity', fontsize = 20)
    plt.tight_layout()
    plt.savefig(f"regression{i}.png")




plt.figure(figsize = (8,8))
plt.title('Feature Importances')
plt.barh(range(len(regression_importance)), width = regression_importance["importance"], height = 0.2, color='#8f63f4', align='center')
plt.yticks(range(len(regression_importance)), regression_importance["features"])
plt.xlabel('Relative Importance')
plt.gca().invert_yaxis()
plt.tight_layout()
plt.savefig("regression_importance.png")

classification[['ECCEN', 'SpMax8_Bhp', 'SpMax7_Bhe', 'SpDiam_Dzi', 'ATS8p']] = classification[['ECCEN', 'SpMax8_Bhp', 'SpMax7_Bhe', 'SpDiam_Dzi', 'ATS8p']].apply(lambda x: np.log10(x))


for i, var in enumerate(classification.columns[:5]):
    plt.figure(figsize = (10,6))
    ax = sns.jointplot(data=classification, x=var, y="affinity", color = '#FE654F')
    ax.fig.suptitle(f"r = {round(np.corrcoef(classification['affinity'], classification[var])[1,0],3)}", size=16)
    ax.set_axis_labels(f"log10("+var+")",'affinity', fontsize = 20)
    plt.tight_layout()
    plt.savefig(f"classification{i}.png")

plt.figure(figsize = (8,8))
plt.title('Feature Importances')
plt.barh(range(len(classification_importance.iloc[:15])), width = classification_importance.iloc[:15]["importance"], height = 0.5, color='#8f63f4', align='center')
plt.yticks(range(len(classification_importance.iloc[:15])), classification_importance.iloc[:15]["features"])
plt.xlabel('Relative Importance')
plt.gca().invert_yaxis()
plt.tight_layout()
plt.savefig("classification_importance.png")










