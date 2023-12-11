# -*- coding: utf-8 -*-
"""
Created on Mon Nov  6 20:51:38 2023

@author: Dongdong Zhou
"""

import numpy as np
import pandas as pd
from sklearn import svm, metrics
from sklearn.metrics import confusion_matrix, roc_curve
from sklearn.preprocessing import Normalizer
from sklearn.model_selection import GridSearchCV, train_test_split, StratifiedKFold
from sklearn.pipeline import Pipeline
from sklearn.feature_selection import SelectKBest, chi2
from sklearn.metrics import make_scorer, precision_score


precision_scorer=make_scorer(precision_score, zero_division=0)
custom_scoring={"accuracy":"accuracy",
                "precision":precision_scorer,
                "recall":"recall",
                "f1":"f1",
                "roc_auc":"roc_auc" 
    }




filename='FCGR_features.txt'
data = pd.read_table(filename,sep=",",header=None)  # for FCGR features
#data = pd.read_table(filename,sep=",")    # for classical or oc features
alldata=np.array(data.iloc[:,1:])
labels=np.array(data.iloc[:,0])

sum_acc=[]        
sum_acvb=[]
sum_aauc=[]
sum_sen=[]
sum_spe=[]

for i in range(1,101):  
    
    X_train, X_test, y_train, y_test = train_test_split(alldata, labels, test_size=0.2, 
                                                        random_state=None,stratify=labels)  
    inner_cv = StratifiedKFold(n_splits=5, shuffle= True)
    acvb=[]
    acc = []
    aauc = []
    spe=[]
    sen=[]
    c_best=[]    
    score=[]
    aprecision=[]
    arecall=[]
    af1=[]
           
    pipe = Pipeline([
        ('nm', Normalizer() ),
         # ('reduce_dim', SelectKBest(chi2)),  # for FCGR features
            ('classify', svm.LinearSVC(
                            class_weight='balanced',
                max_iter=500000
                ))
    ]    )
    param_grid = [
        {
            'classify__C':np.logspace(-3,3,50)
            }
    ]
    clf = GridSearchCV(
            pipe, 
            param_grid, 
            scoring=custom_scoring,
            n_jobs=-1, 
            cv = inner_cv,
            refit="accuracy"
            )    
    
    clf.fit(X_train, y_train)
    clfb = clf.best_estimator_
    pred = clfb.predict(X_test)
    acvb.append(clf.best_score_)
    acc.append(metrics.accuracy_score(y_test, pred))
    scores = clfb.decision_function(X_test)
    score.append(scores)
    fpr, tpr, threshold = roc_curve(y_test, scores, pos_label=1)
    aauc.append(metrics.auc(fpr,tpr))
    aprecision.append(metrics.precision_score(y_test, pred))
    arecall.append(metrics.recall_score(y_test, pred))
    af1.append(metrics.f1_score(y_test, pred))
    cmat=confusion_matrix(y_test, pred)
    sen.append(cmat[1][1] / np.sum(cmat[1,:]) ) 
    spe.append(cmat[0][0] / np.sum(cmat[0,:]))
        
    sum_acc.append(acc)        
    sum_acvb.append(acvb)  
    sum_aauc.append(aauc)  
    sum_sen.append(sen) 
    sum_spe.append(spe) 
    print(i)


print('accuracy is...',  np.mean(np.array(sum_acc)))
print('acvb is ...',np.mean(np.array(sum_acvb)))
print('auc is ...',np.mean(np.array(sum_aauc)))
print('sensitivity is ...',np.mean(np.array(sum_sen)))
print('specificity is...',np.mean(np.array(sum_spe)))


# results=np.hstack((np.array(sum_acc),np.array(sum_aauc),np.array(sum_sen),np.array(sum_spe)))
# mlname='ML_'+filename
# np.savetxt(mlname, results, delimiter=',')

