import matplotlib
from collections import deque
from collections import Counter
import itertools
import pandas as pd
import math
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import random
from sympy import simplify
#from sympy import *
import warnings
from sympy import symbols
#from sympy import log, exp, I, E, pi
from sympy.parsing.sympy_parser import parse_expr
from sympy import expand
import scipy.stats
from statsmodels.tsa.stattools import acf,pacf
from statsmodels import graphics
import statsmodels.api as sm
from sklearn.metrics import roc_auc_score

from sklearn.model_selection import train_test_split
from sklearn.ensemble import RandomForestRegressor
from sklearn.metrics import mean_absolute_error


warnings.filterwarnings("ignore")


data=pd.read_excel('data_do_25MPa.xlsx')
y=data['Z']
colUse=['T','P','Molar_mass','CO2','N2']
#colUse=['P','T']
X=data[colUse]
train_X, val_X, train_y, val_y = train_test_split(X, y,random_state = 0)
val_y_ROC=val_y

print("n_test = ", len(val_y))
print('n_train = ',len(train_y))
mean=val_y.mean()
print(mean)

y_mean_y_val=list([])
val_y=np.array(val_y)
train_y=np.array(train_y)
for i in range(len(val_y)):
    y_mean_y_val.append(val_y[i]-val_y.mean())
S2_0_val=sum( i*i for i in y_mean_y_val)
y_mean_y_train=list([])
for i in range(len(train_y)):
    y_mean_y_train.append(train_y[i]-train_y.mean())
S2_0_train=sum( i*i for i in y_mean_y_train)
y_mean_y=list([])
for i in range(len(y)):
    y_mean_y.append(y[i]-y.mean())
S2_0=sum( i*i for i in y_mean_y)

WAPE_train=list([])
WAPE_test=list([])
num_tree=100
for i in range(1,num_tree):
    forest_model = RandomForestRegressor(n_estimators=i,min_samples_leaf=1,n_jobs=4,random_state=1)
    forest_model.fit(train_X, train_y)
    y_preds = forest_model.predict(val_X)
    s=0
    val_y=list(val_y)
    e=list([])
    for j in range(len(val_y)):
        s=s+(val_y[j]-y_preds[j])**2
        e.append(val_y[j]-y_preds[j])

    e=np.array(e)
    RMSE=(s/len(val_y))**(1/2)
    R2=1-s/S2_0_val
    MAE=mean_absolute_error(val_y, y_preds)
    WAPE=(MAE/mean)*100
    AARD = (sum(abs(np.array(e))/abs(np.array(val_y))))/len(val_y) * 100
    train_y=np.array(train_y)
    WAPE_train.append((mean_absolute_error(train_y,forest_model.predict(train_X))/train_y.mean())*100)
    WAPE_test.append(WAPE)
    #print(i)
    if i==27:
        
        ################# test data
        print("Метрики качества на тестовых данных: ")
        print('AARD: ', AARD, '%')
        print('MAE: ', MAE)
        print('RMSE: ', RMSE)
        print('WAPE: ', WAPE, '%')
        print('R2: ', R2)
        ################# train data
        train_y=list(train_y)
        e=list([])
        y_preds = forest_model.predict(train_X)
        for j in range(len(train_y)):
            s=s+(train_y[j]-y_preds[j])**2
            e.append(train_y[j]-y_preds[j])

        e=np.array(e)
        RMSE=(s/len(train_y))**(1/2)
        R2=1-s/S2_0_train
        MAE=mean_absolute_error(train_y, y_preds)
        WAPE=(MAE/mean)*100
        AARD = (sum(abs(np.array(e))/abs(np.array(train_y))))/len(train_y) * 100
        print("Метрики качества на обучающих данных")
        print('AARD: ', AARD, '%')
        print('MAE: ', MAE)
        print('RMSE: ', RMSE)
        print('WAPE: ', WAPE, '%')
        print('R2: ', R2)
        ################# full data
        print("Метрики качества на полном наборе данных")
        y=list(y)
        e=list([])
        y_preds = forest_model.predict(X)
        for j in range(len(y)):
            s=s+(y[j]-y_preds[j])**2
            e.append(y[j]-y_preds[j])

        e=np.array(e)
        RMSE=(s/len(y))**(1/2)
        R2=1-s/S2_0
        MAE=mean_absolute_error(y, y_preds)
        WAPE=(MAE/mean)*100
        AARD = (sum(abs(np.array(e))/abs(np.array(y))))/len(y) * 100
        print('AARD: ', AARD, '%')
        print('MAE: ', MAE)
        print('RMSE: ', RMSE)
        print('WAPE: ', WAPE, '%')
        print('R2: ', R2)
        fig2, ax2=plt.subplots()
        ax2.scatter(y,e)
        ax2.set_title('График остатков')
        ax2.set_xlabel('Коэффициент сжимаемости Z')
        ax2.set_ylabel('Абсолютная ошибка')


fig, ax=plt.subplots()
ax.set_title("График метрики WAPE относительно количества деревьев в лесу")
ax.set_xlabel("Количество деревьев")
ax.set_ylabel("WAPE, %")
sns.lineplot(range(1,num_tree),WAPE_train, label="WAPE на обучающей выборке")
sns.lineplot(range(1,num_tree),WAPE_test, label="WAPE на тестовой выборке")
plt.show()
    #print(MAE, RMSE, WAPE, R2)
