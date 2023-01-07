import numpy as np
from sklearn.model_selection import LeaveOneOut
from sklearn.model_selection import GridSearchCV
from sklearn.neighbors import KNeighborsClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.ensemble import RandomForestClassifier
from sklearn.feature_selection import SelectKBest, mutual_info_classif
from sklearn.svm import SVC
from sklearn.utils.fixes import loguniform
from sklearn.metrics import confusion_matrix
from sklearn.utils.fixes import loguniform
from sklearn.neural_network import MLPClassifier
from tqdm import tqdm
import pandas as pd
from sklearn.naive_bayes import GaussianNB
import warnings
warnings.filterwarnings('ignore')
from sklearn.tree import DecisionTreeClassifier
import sys
import os
import seaborn as sns
import matplotlib


def ranking():
    loo = LeaveOneOut()
    sales_data = pd.read_csv("expressao_27686.csv")
    genes = open("genes.txt")
    vetor = str.split(genes.readline(),",")
    saida = open("out.csv","w")
    for w in tqdm(vetor):
        excluir = []
        excluir.append(str(w[:-1]))
        excluir.append(str(w))

        #com tudo
        cols = [col for col in sales_data.columns if col in
        excluir]

        data = sales_data[cols]
        target = sales_data['y']
        X = data.to_numpy()
        y = target.to_numpy()
        a_list = list(range(1, 100))


        param_grid = [{'n_estimators': [50, 100, 200],
                                'criterion': ['gini', 'entropy'],
                                'max_depth': np.arange(3, 10),
                                'min_samples_split' :[0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9]}]

        param_grid = [{'n_estimators': [50],
                                'criterion': ['gini', 'entropy'],
                                'max_depth': np.arange(3, 6),
                                'min_samples_split' :[0.1,0.2]}]

        param_grid = [{'C': [10, 100, 1000], 'gamma': [0.0001, 0.001],'kernel': ['rbf','linear']}]
        posi = 0
        X_new = X #SelectKBest(chi2, k=i).fit_transform(X, y)
        tn1=0
        fp1=0
        fn1=0
        tp1=0
        rfF={}
        for train_index, test_index in loo.split(X):
            X_train, X_test = X_new[train_index], X_new[test_index]
            y_train, y_test = y[train_index], y[test_index]
            model =  SVC()
            pred=model.fit(X_train, y_train).predict(X_test)
            tn1,fp1,fn1,tp1 = (confusion_matrix(y_test,pred,labels=[0,1]).ravel())
            rfF[posi]=[tn1,fp1,fn1,tp1]
            posi = posi+1
        TN = 0
        FP = 0
        FN = 0
        TP = 0
        for j in rfF.keys():
            TN = TN + rfF[j][0]
            FP = FP + rfF[j][1]
            FN = FN + rfF[j][2]
            TP = TP + rfF[j][3]
        saida.write(w + "," + str((TP+TN)/(TN+FP+FN+TP))+"\n")
    saida.close()



a = []

rankfile = open("out.csv")
for i in rankfile.readlines():
    data = (str.split(i,","))
    if float(data[1]) >= 1.0:
        a.append(data[0])
rankfile.close()



if len(a)>=10-1:
    matplotlib.use('Agg')
        # Read in the data with `read_csv()`
    sales_data = pd.read_csv("expressao_27686.csv")
    cols = [col for col in sales_data.columns if col in a]
    data = sales_data[cols]
    species = sales_data["y"]

    lut = dict(zip(species.unique(), "rgb"))
    row_colors = species.map(lut)
    h = sns.clustermap(data,  cmap="coolwarm", row_colors=row_colors)
    h.savefig(str(len(a))+".png")
saida = open("27K-ranking.txt","w")
for j in a:
    saida.write(j+"\n")
