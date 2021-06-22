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
from tensorflow import keras
from tensorflow.keras import layers
from sklearn.model_selection import train_test_split
from sklearn.metrics import mean_absolute_error

warnings.filterwarnings("ignore")


data=pd.read_excel('data_W.xlsx')
y=data['W, m/s']
colUse=['T, K','PRES','Molar mass, g/mol']
X=data[colUse]

######### data test
data_test=pd.read_excel('data_W.xlsx')
y_test=data_test['W, m/s']
col_test=['T, K','PRES','Molar mass, g/mol']
X_test=data_test[col_test]

###################
X_train, X_valid, y_train, y_valid = train_test_split(X, y,random_state = 0)

mean=y_valid.mean()
print(mean)



train_X=X_train
train_y=y_train
val_X=X_valid
val_y=y_valid
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


model = keras.Sequential([
    layers.Dense(64, activation='tanh', input_shape=[3]),
    layers.Dense(32, activation='elu'),
    layers.Dense(64, activation='relu'),
    layers.Dense(32, activation='elu'),
    layers.Dense(64, activation='relu'),
    layers.Dense(1),
])

model.compile(
    optimizer='adam',
    loss='mae',
)

history = model.fit(
    X_train, y_train,
    validation_data=(X_valid, y_valid),
    batch_size=256,
    epochs=200,
)

# convert the training history to a dataframe
history_df = pd.DataFrame(history.history)
# use Pandas native plot method
history_df['loss'].plot();

#####################
##X_valid=X_test
##y_valid=y_test
#####################

predict_data=model.predict(X_valid)

s=0
y_valid=list(y_valid)
e=list([])
##predict_data=predict_data
p=list([])
for i in range(len(y_valid)):
    d=predict_data[i]
    p.append(d[0])
predict_data=p
##print(predict_data)
for i in range(len(y_valid)):
    s=s+(y_valid[i]-predict_data[i])**2
    e.append(y_valid[i]-predict_data[i])
##print(predict_data)
RMSE=(s/len(y_valid))**(1/2)
R2=1-s/S2_0_val
MAE=mean_absolute_error(y_valid, predict_data)
WAPE=(MAE/mean)*100
AARD = (sum(abs(np.array(e))/abs(np.array(y_valid))))/len(y_valid) * 100

print("Метрики качества на тестовых данных: ")
print('AARD: ', AARD, '%')
print('MAE: ', MAE)
print('RMSE: ', RMSE)
print('WAPE: ', WAPE, '%')
print('R2: ', R2)

################# train data

train_y=list(train_y)
e=list([])
y_preds = model.predict(train_X)
p=list([])
for i in range(len(y_train)):
    d=y_preds[i]
    p.append(d[0])
y_preds=p
s=0
print(len(train_y))
print(len(y_preds))
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

print("Метрики качества на полном наборе данных")
y=list(y)
e=list([])
s=0
y_preds = model.predict(X)
p=list([])
for i in range(len(y)):
    d=y_preds[i]
    p.append(d[0])
y_preds=p
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
ax2.set_xlabel('Скорость звука, м/с')
ax2.set_ylabel('Абсолютная ошибка, м/с')

plt.show()
