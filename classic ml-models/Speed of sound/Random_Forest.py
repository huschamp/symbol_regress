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

from sklearn.model_selection import train_test_split
from sklearn.ensemble import RandomForestRegressor
from sklearn.metrics import mean_absolute_error


warnings.filterwarnings("ignore")


data=pd.read_excel('data_W.xlsx')
y=data['W, m/s']
colUse=['T, K','PRES','Molar mass, g/mol']
X=data[colUse]
train_X, val_X, train_y, val_y = train_test_split(X, y,random_state = 0)

mean=val_y.mean()
print(mean)

y_mean_y=list([])
for i in range(len(y)):
    y_mean_y.append(y[i]-mean)
S2_0=sum( i*i for i in y_mean_y)


forest_model = RandomForestRegressor(random_state=1)
forest_model.fit(train_X, train_y)
y_preds = forest_model.predict(val_X)
s=0
val_y=list(val_y)
e=list([])
for i in range(len(val_y)):
    s=s+(val_y[i]-y_preds[i])**2
    e.append(val_y[i]-y_preds[i])

e=np.array(e)
RMSE=(s/len(val_y))**(1/2)
R2=1-s/S2_0
MAE=mean_absolute_error(val_y, y_preds)
WAPE=(MAE/mean)*100
AARD = (sum(abs(np.array(e))/abs(np.array(val_y))))/len(val_y) * 100


print('AARD: ', AARD, '%')
print('MAE: ', MAE)
print('RMSE: ', RMSE)
print('WAPE: ', WAPE, '%')
print('R2: ', R2)
#print(MAE, RMSE, WAPE, R2)
