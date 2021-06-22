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
warnings.filterwarnings("ignore")

unknown=['x1','x2','x3']
binary_op=['+','-','*','/','^']
unar_op=['U-','ln','exp','tanh']
a=random.randint(0,25)
b=random.randint(a+1,55)
start_length=15
max_length=400

num_eqn=3 #6 - norm
#num_generation=5
num_generation=30 #20-30 norm
#num_generation=3

data_full=pd.read_excel('data_W.xlsx')
data_full=pd.read_excel('data_train.xlsx')
#data=data_full[(data_full.index % 5) !=0]
#data_full=pd.read_excel('data_norm_alpha007.xlsx')
#data_full=pd.read_excel('data_delete.xlsx')
##data=data_full.sample(frac=0.8)
data=data_full

x1=data['T, K']
x2=data['PRES']
x3=data['Molar mass, g/mol']


y=data['W, m/s']

##data_test=data_full[data_full.index.map(lambda x: x not in data.index)]
##data_test.to_excel('data_test_P.xlsx')
##
##x1.index=np.arange(len(x1))
##x2.index=np.arange(len(x2))
##x3.index=np.arange(len(x3))
##x4.index=np.arange(len(x4))
##y.index=np.arange(len(y))
##x1=np.array(x1)
##x2=np.array(x2)
##x3=np.array(x3)
##x4=np.array(x4)
##y=np.array(y)

#X=np.array([(x1*0+1),x1,np.power(x1,2),np.sin(x1),x2,x3,x4]).T

##data.to_excel('data_train_P.xlsx')

mean=sum(abs(y))/len(y)
print('Среднее значение: ', mean)

y_mean_y=list([])
for i in range(len(y)):
    y_mean_y.append(y[i]-mean)
S2_0=sum( i*i for i in y_mean_y)
print(S2_0)
print(np.dot(np.array(y_mean_y).T,np.array(y_mean_y)))



def del_exp(string):
    string_new=str(string)
    index1=0
    num=[]
    for i in range(len(string)-1):
        if string[i]=='e' and string[i+1]!='x':
            index12=i-1
            if string[i+1] in ['-','+']:#=='-':
                for j in range(i+2,len(string)):
                    if is_digit(string[j]):
                        continue
                    else:
                        index2=j-1
                        break
            else:
                for j in range(i+1,len(string)):
                    if is_digit(string[j]):
                        continue
                    else:
                        index2=j-1
                        break
            for j in range(index12,index1,-1):
                if is_digit(string[j]) or string[j]=='.':
                    continue
                else:
                    index1=j+1
                    break
            num.append(string[index1:index2+1])
            print(num)
            #num_new=format(num, '.f')
    for i in range(len(num)):
        zamena=format(float(num[i]),'.15f')
        print('zamena: ',zamena)
        print(num[i])
        string_new=string_new.replace(num[i],zamena)
    return string_new
                    



def is_digit(string):
    if string.isdigit():
       return True
    else:
        try:
            float(string)
            return True
        except ValueError:
            return False


##def is_int(str):
##    try:
##        int(str)
##        return True
##    except ValueError:
##        return False

def is_int(str):
    num=float(str)
    if int(num)==num:
        return True
    else:
        return False
#print(int(2.01)==2.01)

def prioritet(operation):
    p=0
    if operation=='^':
        p=3
    elif operation=='*' or operation=='/':
        p=2
    elif operation=='-' or operation=='+':
        p=1
    return p


def eqn_to(eqn):
    eqn=str(eqn)
    eqn=eqn.replace('--','+')
    eqn=eqn.replace('**','^')
    eqn=eqn.replace('log','ln')
    eqn=eqn.replace(' ','')
    eqn=eqn.replace('E','2.718281828459')
    eqn=del_exp(eqn)
##    eqn=eqn.replace('e-','*10^-')
##    eqn=eqn.replace('e+','*10^')
    #print('eqn_to: ',eqn)
    stack=deque()
    exit_str=deque()
    num=''
    i=0
    while i<len(eqn):
        if is_digit(eqn[i]):
            num=''
            while (eqn[i] not in binary_op) and eqn[i]!=')':
                num=num+eqn[i]
                i=i+1
                if i==len(eqn):  #-1
                    break
            exit_str.append(str(num))
        elif eqn[i]=='x':
            exit_str.append((eqn[i]+eqn[i+1]))
            i=i+2
        elif eqn[i] in ['U', 'l', 'e','t']:
            if eqn[i]+eqn[i+1] in ['U-', 'ln']:
                stack.append((eqn[i]+eqn[i+1]))
                i=i+2
            elif eqn[i]+eqn[i+1]+eqn[i+2] in ['exp','tanh']:
                stack.append((eqn[i]+eqn[i+1]+eqn[i+2]))
                i=i+3
        elif eqn[i]=='(':
            stack.append(eqn[i])
            i=i+1
        elif eqn[i]==')':
            while stack[-1]!='(':
                exit_str.append(stack[-1])
                stack.pop()
            stack.pop()
            i=i+1
        elif eqn[i] in binary_op:
            if i==0 and eqn[i]=='-':
                stack.append('U-')
                i=i+1
            elif eqn[i]=='-' and eqn[i-1] in ['/','*','(','^']:
                stack.append('U-')
                i=i+1
            else:
                if len(stack)!=0:
                    while (stack[-1] in unar_op) or (prioritet(stack[-1])>prioritet(eqn[i])) or ((prioritet(stack[-1])==prioritet(eqn[i]) and (stack[-1]=='-' or stack[-1]=='/'))):
                        exit_str.append(stack[-1])
                        stack.pop()
                        if len(stack)==0:
                            break
                stack.append(eqn[i])
                i=i+1
    for i in range(len(stack)):
        exit_str.append(stack[-1])
        stack.pop()
    return exit_str




def rss(chromosoma):
    RSS=0
    for i in range(len(y)):
        stack=[]
        for gen in chromosoma:
            if gen in binary_op and gen!='^':
                try:
                    if stack[-2]=='inf':
                        print(stack,'1')
                        RSS=float('Inf')
                        return RSS
                    elem=eval(str(stack[-2])+gen+str(stack[-1]))
                    stack.pop()
                    stack.pop()
                    stack.append(elem)
                except:
                    RSS=float('Inf')
                    print(stack,'2')
                    return RSS
                    #print(gen, 'problem')
                    #break
            elif gen=='^':
                if stack[-2]<0 and is_int(str(stack[-1]))==False:
                    RSS=float('Inf')
                    print(stack,'3')
                    return RSS
                try:
                    elem=pow(stack[-2],stack[-1])
                    if elem=='inf':
                        print(stack,'4')
                        RSS=float('Inf')
                        return RSS
                    stack.pop()
                    stack.pop()
                    stack.append(elem)
                except:
                    print(stack,'5')
                    RSS=float('Inf')
                    return RSS
                    #print(gen, 'problem')
                    #break
            elif gen in unar_op:
                if gen=='U-':
                    stack[-1]=-stack[-1]
                if gen=='ln':
                    if stack[-1]<=0.0000000000001:
                        RSS=float('Inf')
                        print(stack,'6')
                        return RSS
                    try:
                        elem=math.log(stack[-1])
                        stack.pop()
                        stack.append(elem)
                    except:
                        RSS=float('Inf')
                        print(stack,'7')
                        return RSS
                        #print('ln problem')
                        #break
                if gen=='exp':
                    try:
                        elem=math.exp(stack[-1])
                        if elem=='inf':
                            RSS=float('Inf')
                            print(stack,'8')
                            return RSS 
                        stack.pop()
                        stack.append(elem)
                    except:
                        RSS=float('Inf')
                        print(stack,'9')
                        return RSS
                        #print('exp problem')
                        #break
                if gen=='tanh':
                    try:
                        elem=math.tanh(stack[-1])
                        stack.pop()
                        stack.append(elem)
                    except:
                        RSS=float('Inf')
                        print(stack,'10')
                        return RSS
                        #print('exp problem')
                        #break
            else:
                if is_digit(gen):
                    stack.append(float(gen))
                else:
                    stack.append(eval(gen+'[i]'))
        if type(stack[0])==complex:
            RSS=float('Inf')
            print(stack,'11')
            return RSS
        #print(type(stack[0]))
        if stack[0]<0:
            RSS=float('Inf')
            print(stack,'12')
            return RSS
        try:
            RSS=RSS+abs(stack[0]-y[i])
        except:
            RSS=float('Inf')
            print(stack,'13')
            return RSS
    if RSS<0:
        RSS=float('Inf')
        print(stack)
    try:
        if math.isnan(RSS):
            RSS=float('Inf')
            print(stack)
    except:
        RSS=float('Inf')
        print(stack,'14')
    return RSS

def predict(chromosoma,y):
    pred=list([])
    #print(rss(chromosoma))
    for i in range(len(y)):
        stack=[]
        for gen in chromosoma:
            
            if gen in binary_op and gen!='^':
                try:
                    elem=eval(str(stack[-2])+gen+str(stack[-1]))
                    stack.pop()
                    stack.pop()
                    stack.append(elem)
                except:
                    print(stack)
                    print(gen)
                    print(chromosoma)
            elif gen=='^':
                if stack[-2]<0 and is_int(str(stack[-1]))==False:
                    print('ОШИБКА в ^')
                    print(stack)
                elem=pow(stack[-2],stack[-1])
                stack.pop()
                stack.pop()
                stack.append(elem)
            elif gen in unar_op:
                if gen=='U-':
                    stack[-1]=-stack[-1]
                if gen=='ln':
                    try:
                        elem=math.log(stack[-1])
                        stack.pop()
                        stack.append(elem)
                    except:
                        print('ОШИБКА в ln')
                        print(stack[-1])
                if gen=='exp':
                    #print(stack)
                    try:
                        elem=math.exp(stack[-1])
                        stack.pop()
                        stack.append(elem)
                    except:
                        print('ОШИБКА в exp')
                if gen=='tanh':
                    print(stack[-1])
                    elem=math.tanh(stack[-1])
                    stack.pop()
                    stack.append(elem)
            else:
                if is_digit(gen):
                    stack.append(float(gen))
                else:
                    #print(gen)
                    #print(i)
                    stack.append(eval(gen+'[i]'))
        pred.append(stack[0])
    return np.array(pred)


#def predict_Z(EQN,x11,x22,x33,x44,x55,x66,x77,x88,x99):
    #predict=list([])
    #for j in range(len(y)):
    #    eqn=EQN
    #    digit=eqn.subs(x1,x11[j])
    #    digit=digit.subs(x2,x22[j])
    #    digit=digit.subs(x3,x33[j])
    #    digit=digit.subs(x4,x44[j])
    #    digit=digit.subs(x5,x55[j])
    #    digit=digit.subs(x6,x66[j])
    #    digit=digit.subs(x7,x77[j])
    #    digit=digit.subs(x8,x88[j])
    #    digit=digit.subs(x9,x99[j])
    #    predict.append(eval(str(simplify(digit))))
    #return predict

#def predict_P(eqn,x11,x22,x33,x44):
#   predict=list([])
#    for j in range(len(y)):
#        digit=eqn.subs(x1,x11[j])
#        digit=digit.subs(x2,x22[j])
#        digit=digit.subs(x3,x33[j])
#        digit=digit.subs(x4,x44[j])
#        predict.append(eval(str(simplify(digit))))
#    return predict

def MAE(predict,y):
    ae=list([])
    for i in range(len(y)):
        ae.append(y[i]-predict[i])
    mae=sum(abs(np.array(ae)))/len(ae)
    accurate_mae=(1-mae/mean)*100
    return ae,mae,accurate_mae

def MSE(predict,y):
    se=list([])
    for i in range(len(y)):
        se.append((y[i]-predict[i])**2)
    mse=(sum(abs(np.array(se)))/len(se))**(1/2)
    accurate_mse=(1-mse/mean)
    return se,mse,accurate_mse

def S2(ae):
    S2_sum=0
    for i in range(len(ae)):
        S2_sum=S2_sum+(ae[i])**2
    return S2_sum

def R2(S2_0,S2):
    #print(type(S2_0))
    R2_coef=1-(S2/S2_0)
    #print(type(R2_coef))
    return R2_coef


def split_eqn(eqn):
    eqn=str(eqn)
    eqn=parse_expr(eqn)
    eqn=str(expand(eqn))
    #print('IN SPLIT DEF: ')
    #print(eqn)
    #eqn=eqn.replace('**','^')
    eqn=eqn.replace('log','ln')
    eqn=eqn.replace(' ','')
    eqn=eqn.replace('E','2.718281828459')
    eqn=del_exp(eqn)
##
##    eqn=eqn.replace('e+16','*10^(16)')
##    eqn=eqn.replace('e+17','*10^(17)')
##    eqn=eqn.replace('e+18','*10^(18)')
##    eqn=eqn.replace('e+19','*10^(19)')
##    eqn=eqn.replace('e+20','*10^(20)')
##    eqn=eqn.replace('e+21','*10^(21)')
##    eqn=eqn.replace('e+22','*10^(22)')
##    eqn=eqn.replace('e+23','*10^(23)')
##    eqn=eqn.replace('e+24','*10^(24)')
##    eqn=eqn.replace('e+25','*10^(25)')
##    eqn=eqn.replace('e+26','*10^(26)')
##    eqn=eqn.replace('e+27','*10^(27)')
##    eqn=eqn.replace('e+28','*10^(28)')
##    eqn=eqn.replace('e+29','*10^(29)')
##    eqn=eqn.replace('e+30','*10^(30)')
##    eqn=eqn.replace('e+31','*10^(31)')
##    eqn=eqn.replace('e+32','*10^(32)')
##    eqn=eqn.replace('e+33','*10^(33)')
##    eqn=eqn.replace('e+34','*10^(34)')
##
##    eqn=eqn.replace('e-16','*10^(-16)')
##    eqn=eqn.replace('e-17','*10^(-17)')
##    eqn=eqn.replace('e-18','*10^(-18)')
##    eqn=eqn.replace('e-19','*10^(-19)')
##    eqn=eqn.replace('e-20','*10^(-20)')
##    eqn=eqn.replace('e-21','*10^(-21)')
##    eqn=eqn.replace('e-22','*10^(-22)')
##    eqn=eqn.replace('e-23','*10^(-23)')
##    eqn=eqn.replace('e-24','*10^(-24)')
##    eqn=eqn.replace('e-25','*10^(-25)')
##    eqn=eqn.replace('e-26','*10^(-26)')
##    eqn=eqn.replace('e-27','*10^(-27)')
##    eqn=eqn.replace('e-28','*10^(-28)')
##    eqn=eqn.replace('e-29','*10^(-29)')
##    eqn=eqn.replace('e-30','*10^(-30)')
##    eqn=eqn.replace('e-31','*10^(-31)')
##    eqn=eqn.replace('e-32','*10^(-32)')
##    eqn=eqn.replace('e-33','*10^(-33)')
##    eqn=eqn.replace('e-34','*10^(-34)')
##
##    eqn=eqn.replace('e-1','*10^(-1)')
##    eqn=eqn.replace('e+1','*10^(1)')
##    eqn=eqn.replace('e-2','*10^(-2)')
##    eqn=eqn.replace('e+2','*10^(2)')
##    eqn=eqn.replace('e-3','*10^(-3)')
##    eqn=eqn.replace('e+3','*10^(3)')
##    eqn=eqn.replace('e-4','*10^(-4)')
##    eqn=eqn.replace('e+4','*10^(4)')
##    eqn=eqn.replace('e-5','*10^(-5)')
##    eqn=eqn.replace('e+5','*10^(5)')
##    eqn=eqn.replace('e-6','*10^(-6)')
##    eqn=eqn.replace('e+6','*10^(6)')
##    eqn=eqn.replace('e-7','*10^(-7)')
##    eqn=eqn.replace('e+7','*10^(7)')
##    eqn=eqn.replace('e-8','*10^(-8)')
##    eqn=eqn.replace('e+8','*10^(8)')
##    eqn=eqn.replace('e-9','*10^(-9)')
##    eqn=eqn.replace('e+9','*10^(9)')
##    eqn=eqn.replace('e-10','*10^(-10)')
##    eqn=eqn.replace('e+10','*10^(10)')
##    eqn=eqn.replace('e-11','*10^(-11)')
##    eqn=eqn.replace('e+11','*10^(11)')
##    eqn=eqn.replace('e-12','*10^(-12)')
##    eqn=eqn.replace('e+12','*10^(12)')
##    eqn=eqn.replace('e-13','*10^(-13)')
##    eqn=eqn.replace('e+13','*10^(13)')
##    eqn=eqn.replace('e-14','*10^(-14)')
##    eqn=eqn.replace('e+14','*10^(14)')
##    eqn=eqn.replace('e-15','*10^(-15)')
##    eqn=eqn.replace('e+15','*10^(15)')
    print(eqn)
    count=0
    list_eqn=list([])
    j=0
    for i in range(1,len(eqn)):  #1
        if eqn[i]=='(':
            count=count+1
        elif eqn[i]==')':
            count=count-1
        if eqn[i]=='-' and eqn[i-1] in ['e', '^']:
            continue
        if eqn[i] in ['+','-'] and count==0:
            for q in range(j,i):
                s=''
                if eqn[q]=='*': #and (eqn[q+1]!='(' and eqn[q+2]!='1' and eqn[q+3]!='0'):
                    jold=j
                    j=q+1
                    if eqn[q-3]=='x':
                        j=jold
                    s=''
                    count=0
                    break
                elif eqn[q]=='/':
                    j=q
                    s='1'
                    count=0
                    break
##                elif eqn[q]=='^':
##                    
            #print(eqn[i])
            #print(eqn[j])
            #print('list: ')
            #print(list_eqn)
            list_eqn.append(s+eqn[j:i])
            j=i+1
##        elif eqn[i]=='^' and count==0:
##            for q in range(j,i):
##                if eqn[q]=='*' and (eqn[q+1]!='(' and eqn[q+2]!='1'):
##                    j=q+1
##                    s=''
    for q in range(j,len(eqn)):
        if eqn[q]=='*':
            j=q+1
            s=''
            break
        elif eqn[q]=='/':
            j=q
            
            s='1'
            break
    list_eqn.append(s+eqn[j:len(eqn)])
    #print(list_eqn)
    for i in range(len(list_eqn)):
        if is_digit(list_eqn[i]):
            list_eqn[i]='1'
    i=0
    vec=[]
    print('in split_eqn: ',list_eqn)
    while i<len(list_eqn):
        vec=(predict(eqn_to(list_eqn[i]),y))
        set_vec=set(vec)
        if len(set_vec)==1 and vec[0]==0:
            list_eqn.pop(i)
        i=i+1
    return list_eqn


##########################

class EX(object):
    def __init__(self,eps,x):
        self.eps=eps
        self.x=x

    def eps(self):
        return self.eps
    def x(self):
        return self.x
    
    def __str__(self):
        return f"{self.eps, self.x}"
    def __repr__(self):
        return f"{self.eps, self.x}"
        #return "" % (self.eps, self.x)


def byX_key(ex):
    return ex.x
##########################

#def sort_EX(EX):
    

def mnk(X,y):
    #print(X)
    b=np.dot(np.dot(np.linalg.inv(np.dot(X,X.T)),X),y)
    y1=np.dot(X.T,b)
    e=y-y1
    #print('ERROR: ',e)
    S2=np.dot(e.T,e)
    sigma2=S2/(np.shape(X)[0]-np.shape(X)[1])
    S0=np.dot((y-np.mean(y)).T,(y-np.mean(y)))
    r2=1-S2/S0
    return y1,b,S2,r2


def M_const(e,X,k,m):
    print(np.shape(X))
    count=0
    for q in range(np.shape(X)[0]):
        x=X[q]
        mas=list([])
        n=len(e)
        for i in range(len(e)):
            mas.append(EX(e[i],x1[i]))
        #repr(mas)
        mas=sorted(mas,key=byX_key)
        e_s=list([])
        for i in range(len(e)):
            e_s.append(mas[i].eps)
        mas=list([])
        kolvo=len(e)//k
        for i in range(k-1):
            mas.append(e_s[i*kolvo:(i+1)*kolvo])
            index=(i+1)*kolvo
        mas.append(e_s[index:len(e_s)])
        mean_e=sum(e)/len(e)
        S0=0
        for j in range(k):
            mean_group=sum(mas[j])/len(mas[j])
            S0=S0+(mean_group-mean_e)**2
        S=0
        for j in range(k):
            S_v=0
            mean_group=sum(mas[j])/len(mas[j])
            for i in range(len(mas[j])):
                S_v=S_v+(mas[j][i]-mean_group)**2
            S=S+S_v
        Criter=(S0/(k-m))/(S/(n-k))
        F=scipy.stats.f.ppf(q=1-0.05,dfn=k-m,dfd=n-k)
        #print(Criter)
        #F=scipy.stats.f.cdf(Criter,k-m,n-k)
##        if Criter>F:
##            print('Математическое ожидание остатков не постоянно')
##        else:
##            print('Математическое ожидание остатков постоянно')
        if Criter>F:
            count=count
        else:
            count=count+1
        #print(count)
    if count==np.shape(X)[0]:
        print('Математическое ожидание остатков постоянно')
        return 0
    elif count>1 and count<np.shape(X)[0]:
        print('Математическое ожидание остатков постоянно не по всем факторам')
        return 1
    else:
        print('Математическое ожидание остатков не постоянно')
        return 1




def Crit_White(e,X):
    #print(np.shape(X))
    ######
    e_new=list([])
    for i in range(len(e)):
        e_new.append(e[i]**2)
    e=e_new
    ######
    x1=X
    #for q in range(np.shape(X)[0]):
    #x=X[q]
    #print(x)
    alpha=0.05
    X_vsp=list([])
    for i in range(np.shape(X)[0]):
        if i==0:
            X_vsp.append(X[i]**0)
        X_vsp.append(X[i])
        #print(X_vsp)
        X_vsp.append(X[i]**2)
    for i in range(np.shape(X)[0]):
        for j in range(np.shape(X)[0]):
            if i!=j:
                X_vsp.append(X[i]*X[j])
    #for i in range(np.shape(X_vsp)[0]):
    #    for j in range(np.shape(X_vsp)[0]):
    #        if i!=j and X_vsp[i]==X_vsp[j]:
    #            X_vsp.pop(j)
    #print(np.shape(X_vsp))
    cf=[]
    for i in range(np.shape(X_vsp)[0]):
        cf.append(sum(X_vsp[i])/len(X_vsp[i]))
    X_vsp2=list([])
    cf2=list([])
    for i in range(len(cf)):
        if cf[i] not in cf2:
            X_vsp2.append(X_vsp[i])
            cf2.append(cf[i])
    #X_vsp=set(X_vsp)
    X_vsp=np.array(X_vsp2)
    #X_vsp=np.array([x1**0,x1,x2,x3,x4,x1**2,x2**2,x3**2,x4**2,x1*x2,x1*x3,x1*x4,x2*x3,x2*x4,x3*x4])
    #i=0
    #while i!=1:
        #X_vsp.append(X**0)
        #i=1
    #X_vsp=np.array([x**0,x,np.dot(np.dot(x,A),x)])
    #print(np.shape(X_vsp))
    n_vsp=np.shape(X_vsp)[1]
    m_vsp=np.shape(X_vsp)[0]
    #print(m_vsp)
    y1_vsp,b_vsp,S2_vsp,R2_vsp=mnk(X_vsp,e)
    X0=list([])
    X0.append(X[0]**0)
    #print(X0)
    y1_0,b_0,S2_0,R2_0=mnk(np.array(X0),e)
    #####
##    e=np.array(e)
##    me=sum(abs(e))/len(e)
##    delta=list([])
##    for i in range(len(e)):
##        delta.append(e[i]-me)
##    S2_0=0
##    for i in range(len(delta)):
##        S2_0=S2_0+delta[i]**2
    #####
    Statistica=((S2_0-S2_vsp)/(m_vsp-1))/(S2_vsp/(n_vsp-m_vsp))
    F=scipy.stats.f.ppf(q=1-0.05,dfn=m_vsp-1,dfd=n_vsp-m_vsp)
    LM=len(e)*(1-S2_vsp/S2_0)
    Chi2=scipy.stats.chi2.ppf(q=1-0.05,df=m_vsp-1)
    #print(LM)
    #print(Chi2)
    #print(Statistica)
    #print(F)
    if abs(Statistica)>F and abs(LM)>Chi2:
        print('Дисперсия остатков не постоянна по критерию Уайта')
    else:
        print('Дисперсия остатков постоянна по критерию Уайта')


def Crit_Goldfeld(e,X):
    count=0
    k=3
    e_new=list([])
    for i in range(len(e)):
        e_new.append(e[i]**2)
    e=e_new
    X_copy=list(X)
    for q in range(np.shape(X)[0]):
        x=X[q]
        mas=list([])
        n=len(e)
        for i in range(len(e)):
            mas.append(EX(e[i],x1[i]))
        #repr(mas)
        mas=sorted(mas,key=byX_key)
        e_s=list([])
        for i in range(len(e)):
            e_s.append(mas[i].eps)
        mas=list([])
        kolvo=len(e)//k
        for i in range(k-1):
            mas.append(e_s[i*kolvo:(i+1)*kolvo])
            index=(i+1)*kolvo
        mas.append(e_s[index:len(e_s)])
        n1=len(mas[0])
        X=list(X)
        for i in range(np.shape(X)[0]):
            X_X=list(X[i])
            #print(len(X_X))
            X[i]=X_X[0:n1]
        #print(np.shape(X))
        X=np.array(X)
        alpha=0.05
        X_vsp=list([])
        for i in range(np.shape(X)[0]):
            if i==0:
                X_vsp.append(X[i]**0)
            X_vsp.append(X[i])
            X_vsp.append(X[i]**2)
        for i in range(np.shape(X)[0]):
            for j in range(np.shape(X)[0]):
                if i!=j:
                    X_vsp.append(X[i]*X[j])

        cf=[]
        for i in range(np.shape(X_vsp)[0]):
            cf.append(sum(X_vsp[i])/len(X_vsp[i]))
        X_vsp2=list([])
        cf2=list([])
        for i in range(len(cf)):
            if cf[i] not in cf2:
                X_vsp2.append(X_vsp[i])
                cf2.append(cf[i])
        X_vsp=np.array(X_vsp2)

        n_vsp=np.shape(X_vsp)[1]
        m_vsp=np.shape(X_vsp)[0]
        y1_vsp,b_vsp,S2_vsp,R2_vsp=mnk(X_vsp,mas[0])


        X=X_copy
        n2=len(mas[1])
        X=list(X)
        for i in range(np.shape(X)[0]):
            X_X=list(X[i])
            #print(len(X_X))
            X[i]=X_X[n1:n1+n2]
        #print(np.shape(X))
        X=np.array(X)
        alpha=0.05
        X_vsp=list([])
        for i in range(np.shape(X)[0]):
            if i==0:
                X_vsp.append(X[i]**0)
            X_vsp.append(X[i])
            X_vsp.append(X[i]**2)
        for i in range(np.shape(X)[0]):
            for j in range(np.shape(X)[0]):
                if i!=j:
                    X_vsp.append(X[i]*X[j])

        cf=[]
        for i in range(np.shape(X_vsp)[0]):
            cf.append(sum(X_vsp[i])/len(X_vsp[i]))
        X_vsp2=list([])
        cf2=list([])
        for i in range(len(cf)):
            if cf[i] not in cf2:
                X_vsp2.append(X_vsp[i])
                cf2.append(cf[i])
        X_vsp=np.array(X_vsp2)

        n_vsp2=np.shape(X_vsp)[1]
        m_vsp2=np.shape(X_vsp)[0]
        y1_vsp2,b_vsp2,S2_vsp2,R2_vsp2=mnk(X_vsp,mas[1])

        Statistica=(S2_vsp/(n1-m_vsp))/(S2_vsp2/(n2-m_vsp2))
        F=scipy.stats.f.ppf(q=1-0.05,dfn=n1-m_vsp,dfd=n2-m_vsp2)
        #print(F)
        #print('F: ', Statistica)
##        if Statistica>F:
##            print('Дисперсия остатков не постоянна по критерию Голдфелда-Куандта')
##        else:
##            print('Дисперсия остатков постоянна по критерию Голдфелда-Куандта')
        if Statistica>F:
            count=count
        else:
            count=count+1
    if count==np.shape(X)[0]:
        print('Дисперсия остатков постоянна по критерию Голдфелда-Куандта')
        return 0
    elif count>1 and count<np.shape(X)[0]:
        print('Дисперсия остатков постоянна не по всем факторам согласно критерию Голдфелда-Куандта')
        return 1
    else:
        print('Дисперсия остатков не постоянна по критерию Голдфелда-Куандта')
        return 1 

        

def Q_Criter_Liunga(e,X,S):
    e_copy=list(e)
    for q in range(np.shape(X)[0]):
    #q=0
        x=X[q]
        mas=list([])
        n=len(e)
        m=np.shape(X)[0]
        for i in range(len(e)):
            mas.append(EX(e[i],x[i]))
        mas=sorted(mas,key=byX_key)
        e_s=list([])
        for i in range(len(e)):
            e_s.append(mas[i].eps)
        e=e_s
        ####
        autocov=list([])
        for k in range(len(e)):
            pr=0
            for i in range(k,n):
                pr=pr+e[i]*e[i-k]
            autocov.append(pr/(n-k))
        
        sigma2=S/(n-m)
        #autocorr=autocov/sigma2
        autocorr=acf(e)
        
        dover_int=list([])
        for i in range(len(autocorr)):
            dover_int.append(2/((n-i)**(1/2)))
        #print('acf: ', autocorr)
        #print('interval: ', dover_int)
        if q!=-1:
            e_frame=pd.DataFrame(e)
            fig = plt.figure(figsize=(12,8))
            ax1 = fig.add_subplot(211)
            fig = sm.graphics.tsa.plot_acf(e_frame.values.squeeze(), lags=15, ax=ax1)
            ax2 = fig.add_subplot(212)
            fig = sm.graphics.tsa.plot_pacf(e_frame, lags=15, ax=ax2)
            
            q_test = sm.tsa.stattools.acf(e_frame, qstat=True)
            #print(q_test[0])
            #print(pd.DataFrame({'Q-stat':q_test[1], 'p-value':q_test[2]}))  # если p-value ниже значения alpha=0.05, то гипотезу - отвергаем
            #plt.show()#  - good
        Q_dop=0
        k=15
        #k=len(autocorr)
        for i in range(1,k):
            Q_dop=Q_dop+(autocorr[i]**2)/(n-i)
        Q=n*(n+2)*Q_dop
        Chi2=scipy.stats.chi2.ppf(q=1-0.05,df=k-m) #- НАВЕРНОЕ ДОЛЖНО БЫТЬ ТАК
        #Chi2=scipy.stats.chi2.ppf(q=1-0.05,df=k-1)
        #print('Q: ',Q)
        #print('Chi2 for Q: ',Chi2)
        p_value=scipy.stats.chi2.sf(Q,df=k-m)
        #print(p_value)
        if p_value<0.05:
            print('Ошибки коррелированны согласно Q-критерию Льюнга-Бокса')
        else:
            print('Ошибки некоррелированны согласно Q-критерию Льюнга-Бокса')
##        if Q<Chi2:
##            print('Ошибки некоррелированны согласно Q-критерию Льюнга-Бокса')
##        else:
##            print('Ошибки коррелированны согласно Q-критерию Льюнга-Бокса')
        

def Criter_Pirsona(e,X,k):
    n=len(e)
    m=np.shape(X)[0]
    ######
##    e_new=list([])
##    for i in range(len(e)):
##        e_new.append(e[i]**3-e[i]**2)
##    e=e_new
    ######
    ####
    #k=7
    #m=2
    ####
    e2=list([])
    for i in range(n):
        e2.append(e[i]**2)
    sigma=(sum(e2)/(n-m))**(1/2)
    e=e/sigma
    e_s=sorted(e)
    mas=list([])
    min_e=min(e)
    max_e=max(e)
    razn=(max_e-min_e)/k
    intervals=list([min_e])
    for i in range(k):
        intervals.append(min_e+razn*(i+1))
        ee=list([])
        for j in range(len(e_s)):
            if i!=k-1:
                if e_s[j]<min_e+razn*(i+1) and e_s[j]>=min_e+razn*i:
                    ee.append(e_s[j])
            else:
                if e_s[j]<=min_e+razn*(i+1) and e_s[j]>=min_e+razn*i:
                    ee.append(e_s[j])
        mas.append(ee)
    #print(intervals)
    p=list([])
    for i in range(len(intervals)-1):
        p.append(scipy.stats.norm.cdf(intervals[i+1])-scipy.stats.norm.cdf(intervals[i]))
    #print(p)
    n_interval=list([])
    for i in range(len(mas)):
        e_interval=mas[i]
        n_interval.append(len(mas[i]))
    npn=list([])
    for i in range(len(n_interval)):
        npn.append(n_interval[i]/n)
    #print('Chi2test: ',scipy.stats.chisquare(npn,f_exp=p,ddof=k-m-2))
    Statistica=0
    for i in range(k):
        Statistica=Statistica+(((n_interval[i]-n*p[i])**2)/(n*p[i]))
    #print(Statistica)
    fig2 = sm.qqplot(pd.DataFrame(e_s), scipy.stats.norm,fit=True,line="45",label='Q-Q Plot')
    #ax2.set_title('Вероятностная бумага')
    fig3 = plt.figure()
    sns_plot=sns.distplot(pd.DataFrame(e_s),fit=scipy.stats.norm)
    Chi2=scipy.stats.chi2.ppf(q=0.975,df=k-m-2)
    p_value=scipy.stats.chi2.sf(Statistica,df=k-m-2)
    print(Statistica)
    print(scipy.stats.chi2.ppf(0.95,k-m-2))
    if p_value<0.05:
        print('Ошибки не согласуются с нормальным распределением согласно Хи^2-критерию Пирсона-Фишера')
    else:
        print('Ошибки имеют нормальное распределение согласно Хи^2-критерию Пирсона-Фишера')
##  Код ниже АНАЛОГИЧЕН
##    print(scipy.stats.chi2.sf(Statistica,df=k-m-2))
##    #print(Chi2)
##    if Statistica<Chi2:
##        print('Ошибки имеют нормальное распределение согласно Хи^2-критерию Пирсона-Фишера')
##    else:
##        print('Ошибки не согласуются с нормальным распределением согласно Хи^2-критерию Пирсона-Фишера')




def minimum(a):
    min = a[0]
    pos = 0
    for i in range(len(a)):
        if a[i]<min:
            min=a[i]
            pos=i
    return [min,pos]


def Crit_Studenta(X,y,b,S,znach):
    X=np.array(X)
    n=np.shape(X)[1]
    m=np.shape(X)[0]
    sigma2=S/(n-m)
    V=np.linalg.inv(np.dot(X,X.T))
    t=list([])
    for j in range(m):
        t.append(b[j]/((sigma2*V[j,j])**(1/2)))
    tnm=scipy.stats.t.ppf(0.95,df=n-m)
    t_index=list([])
    for j in range(m):
        if (t[j]<tnm) and (t[j]>-tnm):
            t_index.append(j)
            #print('Фактор ', j, ' незначим по Стьюденту')
        #else:
            #print('Фактор ', j, ' значим по Стьюденту')
    t_copy=list(t)
    min_t_index=list([])
    if len(znach)==2:
        t_index.pop()
        t_index.pop()
    for i in range(len(t_copy)):
        if i in znach:
            t_copy[i]=float('Inf')
    for j in range(len(t_index)):
        #print(min_t_index)
        min_t_mas=minimum(absl(t_copy))
        #print('znach: ',znach)
        min_t_index.append(min_t_mas[1])
        #t_copy.pop(min_t_mas[1])
        t_copy[min_t_mas[1]]=float('Inf')
    return min_t_index


def Crit_Fishera(X,y,min_t_index):
    X=list(X)
    znach=list([])
    n=np.shape(X)[1]
    m=np.shape(X)[0]
    if len(min_t_index)==0:
        return X,znach
    elif len(min_t_index)==1:
        X_copy=list(X)
        for i in range(m):
            if i in min_t_index:
                X.pop(i)
        return X,znach
    else:
        min_t_index=min_t_index[0:2]
        min_t_index=sorted(min_t_index)
        min_t_index.reverse()
        #print(min_t_index)
        X_2=list(X)
        #print(np.shape(X_2)[0])
        for i in range(m,-1,-1):
            if i in min_t_index:
                X_2.pop(i)
        X=np.array(X)
        X_2=np.array(X_2)
        y1,b,S2,R2=mnk(X,y)
        y1_2,b_2,S2_2,R2_2=mnk(X_2,y)
        X=list(X)
        F_stats=((S2_2-S2)/(2))/((S2)/(n-m))
        F=scipy.stats.f.ppf(0.95,dfn=2,dfd=n-m)
        #print(F_stats)
        #print(F)
        if F_stats<F:
            #print('Факторы ',min_t_index,' незначимы по Фишеру')
            for i in range(m,-1,-1):
                if i in min_t_index:
                    X.pop(i)
            return X,znach
        else:
            #print('Факторы ',min_t_index,' значимы по Фишеру')
            znach=min_t_index
            return X,znach

def absl(x):
    for i in range(len(x)):
        x[i]=abs(x[i])
    return x



def get_eqn(formula):
    #print(formula)
    stack=deque()
    count=0
    for gen in formula:
        #print(stack)
        if gen=='-' and count==0:
            stack.append(gen)
        elif is_digit(gen) or gen in unknown:
            stack.append(gen)
        elif gen in unar_op:
            elem=stack[-1]
            stack.pop()
            if gen!='U-':
                stack.append(gen+'('+elem+')')
            else:
                stack.append('(-'+elem+')')
        elif gen in binary_op:
            elem2=stack[-1]
            elem1=stack[-2]
            stack.pop()
            stack.pop()
            stack.append('('+elem1+gen+elem2+')')
        count=count+1
    eqn=stack[0].replace('--','+')
    eqn=eqn.replace('+-','-')
    eqn=eqn.replace('-+','-')
    eqn=eqn.replace('^','**')
    eqn=eqn.replace('++','+')
##    eqn=stack[0].replace('--','+')
##    eqn=stack[0].replace('+-','-')
##    eqn=stack[0].replace('-+','-')
##    eqn=stack[0].replace('^','**')
    return eqn



def make_formula(b,X,X_copy,split):
    index=list([])
    X=list(X)
    X_copy=list(X_copy)
    for i in range(np.shape(X)[0]):
        for j in range(np.shape(X_copy)[0]):
            x=X[i]
            x_copy=X_copy[j]
            if x[0]==x_copy[0]:
                index.append(j)
    index=list(set(index))
    #print(index)
    split_new=list([])
    for i in range(len(split)):
        if i in index:
            split_new.append(split[i])
    
    index=range(len(b))
    #print('SPLIT NEW:        <', split_new)
    string=''
    for i in index:
        #print(string)
        if i!=0:
            string=string+'+'+str(b[i])+'*'+'('+split_new[i]+')'
        else:
            string=string+str(b[i])+'*'+'('+split_new[i]+')'
    #print('do smipli:', string)
    string=string.replace('+-','-')
    string=string.replace('-+','-')
    string=string.replace('--','+')
    #eqn=simplify(string)
    eqn=string
    #print('posle simpli: ',eqn)
    return str(eqn)


#print('EQQNN: ', get_eqn(eqn_to('2.75*exp(x3)*10^-3*x1')))

sss='0**(1/x1)'
sss=parse_expr(sss)
print('!!!!!!!!!!!!!', expand(sss))


def analys(F,y):
    #print('RSS FOR F: ',rss(eqn_to(F)))
    #print('OBRATNOE: ',get_eqn(eqn_to(get_eqn(eqn_to(F)))))
    split=split_eqn(F)
    print('split: ',split)
    X=list([])
    for i in range(len(split)):     
        X.append(predict(eqn_to(split[i]),y))
    X=np.array(X)
    EXCEL_X=pd.DataFrame(X)
    EXCEL_X.to_excel('matrix_X.xlsx')
    #############################

    print('Матрица эксперимента: ',X)
    polsk=eqn_to(F)
    print(get_eqn(polsk))
    spl=[]
    for i in split:
        spl.append(eqn_to(i))
    length=1
    c=0
    mae_old=999999999
    while c!=1:
        X=[]
        print('len spl',len(spl))
        print(length)
        for i in range(length):
            try:
                X.append(predict(spl[i],y))
            except:
                spl.pop(i)
        #print(np.linalg.det(np.dot(X,X.T)))
        X=np.array(X)
        det=np.linalg.det(np.dot(X,X.T))
        if det==0:
           #length=length+1
           spl.pop(length-1)
           if length>len(spl):
            c=1
            print(mae)
            print(spl)
           continue
        y1,b,S2_2,R2_2=mnk(X,y)
        #mae=np.mean(np.array(y-np.dot(X.T,b)))
        pred_y=np.dot(X.T,b)
        ae,mae,accurate_mae=MAE(pred_y,y)
        #print('pred_y',pred_y)
        print(mae)
        if abs(mae)>abs(mae_old):
            spl.pop(length-1)
            mae=mae_old
        else:
            mae_old=mae
            length=length+1
        if length>len(spl):
            c=1
            print(mae)
            print(spl)
    X=[]
    for i in range(len(spl)):
        X.append(predict(spl[i],y))
    #print(np.linalg.det(np.dot(X,X.T)))
    X=np.array(X)
    y1,b,S2_2,R2_2=mnk(X,y)
    mae=np.mean(np.array(y-np.dot(X.T,b)))
    print(np.shape(X))
    string=''
    num_eqn=len(spl)-1
    for i in range(num_eqn+1):
        #print('')
        #print(best_individual[i])
        #print(get_eqn(best_individual[i]))
        if i!=0:
            string=string+'+'+str(b[i])+'*'+'('+get_eqn(spl[i])+')'
        else:
            string=string+str(b[i])+'*'+'('+get_eqn(spl[i])+')'
    string=string.replace('+-','-')
    string=string.replace('-+','-')
    string=string.replace('--','+')
    string=string.replace('++','+')
    eqn=simplify(string)
    FF=eqn
    F=eqn
    #split=list(spl)
    split=[]
    for i in range(len(spl)):
        split.append(get_eqn(spl[i]))
    print(FF)

    
    #pred_y=predict(eqn_to(FF),y)
    pred_y=np.dot(X.T,b)
    print('pred_y',pred_y)
    ae,mae,accurate_mae=MAE(pred_y,y)
    se,mse,accurate_mse=MSE(pred_y,y)
    S2_f=S2(ae)
    #print(type(S2_0))
    #print(type(S2_f))
    global S2_0
    S2_0=float(S2_0)
    S2_f=float(S2_f)
    print(type(S2_0))
    print(type(S2_f))
    r2=R2(S2_0,S2_f)
    #r2_adj=1-(1-r2)*((n-1)/(n-m))
    print('Прогноз: ', pred_y)
    print('Абсолютные ошибки: ', ae)
    print('Средняя абсолютная ошибка: ',mae)
    print('Cреднеквадратичная ошибка: ', mse)
    print('WAPE: ',100-accurate_mae,'%')
    print('Коэффициент детерминации: ', r2)
    #print('Скорректированный коэффициент детерминации: ', r2_adj)

    #############################
    y1,b,S2_2,R2_2=mnk(X,y)
    print('ПРОГНОЗЗЗЗЗЖ ', y1)
    X_copy=list(X)
    #print('Formula do znac:',F)
    print('Количество факторов до проверки на значимость по критериям Стьюдента и Фишера: ',np.shape(X)[0])
    znach=list([])
    for i in range(50):
        index=Crit_Studenta(X,y,b,S2_2,znach)
        X,znach=Crit_Fishera(X,y,index)
    X=np.array(X)
    y1,b,S2_2,R2_2=mnk(X,y)
    print('Количество факторов после проверки на значимость по критериям Стьюдента и Фишера: ',np.shape(X)[0])
    F=make_formula(b,X,X_copy,split)
    print('FFFFFF')
    print(F)
    #F_view=simplify(F)
    #print('Итоговая формула: ',F_view)

    split=split_eqn(F)
    print(split)
    #print('split posle znach: ',split)
    X=list([])
    for i in range(len(split)):     
        X.append(predict(eqn_to(split[i]),y))
    X=np.array(X)
    print(np.shape(X))
    y1,b,S2_2,R2_2=mnk(X,y)
    n=np.shape(X)[1]
    m=np.shape(X)[0]
    #X=X.T
    print('Матрица эксперимента: ',X)
    pred_y=predict(eqn_to(F),y)
    ae,mae,accurate_mae=MAE(pred_y,y)
    se,mse,accurate_mse=MSE(pred_y,y)
    for i in range(len(ae)):
        if ae[i]>80:
            print(i)
            print(y[i])
#############################
    fig9, ax9=plt.subplots()
    ax9.scatter(y,pred_y)
    ax9.set_title('График остатков')
    ax9.set_xlabel('Предсказанные значения')
    ax9.set_ylabel('Эталонные значения')
    plt.show()
##############################
    S2_f=S2(ae)
    r2=R2(S2_0,S2_f)
    r2_adj=1-(1-r2)*((n-1)/(n-m))
    print('Прогноз: ', pred_y)
    print('Абсолютные ошибки: ', ae)
    print('Средняя абсолютная ошибка: ',mae)
    print('Cреднеквадратичная ошибка: ', mse)
    print('WAPE: ',100-accurate_mae,'%')
    print('Коэффициент детерминации: ', r2)
    print('Скорректированный коэффициент детерминации: ', r2_adj)
    #fig3=plt.figure()
    fig3, ax=plt.subplots()
    ax.scatter(y,ae)
    ax.set_title('График остатков')
    ax.set_xlabel('Скорость звука, м/c')
    ax.set_ylabel('Абсолютная ошибка, м/c')
    ##########################################
    fig4, ax2=plt.subplots()
    ax2.scatter(x1,ae)
    ax2.set_title('Зависимость остатков от T')
    ax2.set_xlabel('Температура')
    ax2.set_ylabel('Абсолютная ошибка')
    fig5, ax3=plt.subplots()
    ax3.scatter(x2,ae)
    ax3.set_title('Зависимость остатков от P')
    ax3.set_xlabel('Давление')
    ax3.set_ylabel('Абсолютная ошибка')
    fig6, ax4=plt.subplots()
    ax4.scatter(x3,ae)
    ax4.set_title('Зависимость остатков от молярной массы')
    ax4.set_xlabel('Молярная масса')
    ax4.set_ylabel('Абсолютная ошибка')
    ##########################################
    #plt.show()
    #возможно, вместо X нужно x1,x2,...,x9
    M_const(ae,X,4,3)
    Crit_White(ae,X)
    Crit_Goldfeld(ae,X)
    Q_Criter_Liunga(ae,X,S2_f)
    #print(np.shape(X)[0])
    #sigma2=S2_f/(n-m)
    #normed_ae=(ae)/sigma2
    #print(scipy.stats.kstest(normed_ae,'norm'))
    #print(sm.stats.lilliefors(ae)[1])
    if sm.stats.lilliefors(ae)[1]>0.05: # тест Лиллиефорса
        print('Ошибки имеют нормальное распределение согласно критерию Лиллиефорса')
    else:
        print('Согласно критерию Лиллиефорса ошибки имеют не нормальное распределение')
    if sm.stats.normal_ad(np.array(ae))[1]>0.05: # тест Андерсона-Дарлинга
        print('Ошибки имеют нормальное распределение согласно критерию Андерсона-Дарлинга')
    else:
        print('Согласно критерию Андерсона-Дарлинга ошибки имеют не нормальное распределение')
    #print('Критерий Лилифорса: ', sm.stats.lilliefors(ae))
    #print('Критерий Андерсона-Дарлинга: ', sm.stats.normal_ad(np.array(ae)))
    Criter_Pirsona(ae,X,np.shape(X)[0]+3) #--- чтобы всегда норм работал, надо склеивать интервалы

    fig_hist=plt.figure()
    sns_plot_p=sns.distplot(pd.DataFrame(ae),fit=scipy.stats.norm, bins=9)


    stat, p = scipy.stats.shapiro(pd.DataFrame(ae)) # тест Шапиро-Уилк
    print('Statistics=%.3f, p-value=%.3f' % (stat, p))
    alpha = 0.05
    if p > alpha:
        print('Распределение нормальное согласно тесту Шапиро-Уилка')
    else:
        print('Распределены не нормальное согласно тесту Шапиро-Уилка')
    
    aepd=pd.DataFrame(ae)
    #aepd.to_excel('ae_P.xlsx') - for 2 formula

    RMSD = math.sqrt(sum(np.array(ae)**2)/n)
    print('RMSD: ', RMSD)

    ARD=list([])
    for i in range(len(ae)):
        ARD.append(abs((ae[i])/(y[i]))*100)
    per=list([])
    for i in range(len(ARD)):
        count=0
        for j in range(len(ARD)):
            if ARD[j]<ARD[i]:
                count=count+1
        per.append((count/len(ARD))*100)

    fig_ard, ax_ard=plt.subplots()
    plt.plot(sorted(ARD),sorted(per))
    ax_ard.set_title('Накопительная частота ошибки ARD')
    ax_ard.set_xlabel('Максимальное значение ARD, (%)')
    ax_ard.set_ylabel('Накопительная частота для скорости звука, (%)')
    #df = pd.DataFrame([(i, s) for i in per for s in ARD], columns=['per', 'ARD'])
    #df.plot(x='ARD', y='per', kind='line', 
    #    figsize=(10, 8), legend=False, style='mo-')
    #plt.show()
    
    AARD = (sum(abs(np.array(ae))/abs(np.array(y))))/n * 100
    print('AARD: ', AARD, '%')

    #standing=list([])
##    ae_standing=list([])
##    standing=(18.2*((x1/x2)**0.83 * ((10**(0.00091*x4))/(10**(0.0125*x3))) -1.4 ))
##    for i in range(n):
##        #standing.append(18.2*((x1/x2)**0.83 * ((10**(0.00091*x4))/(10**(0.0125*x3))) -1.4 ))
##        #print(standing)
##        ae_standing.append(standing[i]-y[i])
##    AARD_standing = (sum(abs(np.array(ae_standing))/abs(np.array(y))))/n * 100
##    print('AARD Standing: ', AARD_standing, '%')
    return F


def clear_unar_minus(chromosoma):
    length=len(chromosoma)
    i=1
    while i<length:
        if chromosoma[i-1]==chromosoma[i] and chromosoma[i]=='U-':
            chromosoma.pop(i)
            chromosoma.pop(i-1)
            length=length-2
        i=i+1
    return chromosoma

######################## процедура МНК для сгенерированных формул
######################## то что надо!!!
def model_of_best_individual(best_individ):
    X=list([])
    #print(np.shape(best_individ))
    for i in range(np.shape(best_individ)[0]): #мб 0 вместо 1
        X.append(predict(best_individ[i],y))
    return np.array(X)


def test(F,data):
    x1=data['T, K']
    x2=data['PRES']
    x3=data['Molar mass, g/mol']
##    x4=data['CO2']
##    x5=data['N2']
    print(len(x1))
    y=data['W, m/s']

    mean=sum(abs(y))/len(y)
    print(mean)
    y_mean_y=list([])
    for i in range(len(y)):
        y_mean_y.append(y[i]-mean)
    S2_0=sum( i*i for i in y_mean_y)
    n=len(y)
    m=15
    #pred_y=predict(eqn_to(F),y)
    

    pred_y=list([])
    x11=list(x1)
    x22=list(x2)
    x33=list(x3)
##    x44=list(x4)
##    x55=list(x5)
    x1,x2,x3=symbols('x1 x2 x3')
    eqn=parse_expr(F)
    for i in range(len(y)):
        eqn_predict=eqn
        eqn_predict=eqn_predict.subs(x1,x11[i])
        eqn_predict=eqn_predict.subs(x2,x22[i])
        eqn_predict=simplify(eqn_predict.subs(x3,x33[i]))
        eqn_predict=str(eqn_predict)
        eqn_predict=eval(eqn_predict)
        print(eqn_predict)
        pred_y.append(eqn_predict)
    ae,mae,accurate_mae=MAE(pred_y,y)
    se,mse,accurate_mse=MSE(pred_y,y)
    S2_f=S2(ae)
    r2=R2(S2_0,S2_f)
    r2_adj=1-(1-r2)*((n-1)/(n-m))
    print('Прогноз: ', pred_y)
    print('Абсолютные ошибки: ', ae)
    print('Средняя абсолютная ошибка: ',mae)
    print('Cреднеквадратичная ошибка: ', mse)
    print('WAPE: ',100-accurate_mae,'%')
    print('Коэффициент детерминации: ', r2)
    print('Скорректированный коэффициент детерминации: ', r2_adj)
    AARD = (sum(abs(np.array(ae))/abs(np.array(y))))/n * 100
    print('AARD: ', AARD, '%')
    fig3, ax=plt.subplots()
    ax.scatter(y,ae)
    ax.set_title('График остатков')
    ax.set_xlabel('Скорость звука, м/c')
    ax.set_ylabel('Абсолютная ошибка, м/c')

###
#best_individual=[['x1', '0', '*', '1', '+'], ['x3', 'ln', 'ln', 'x1', '^', 'ln', 'x3', '*', 'x3', 'x2', 'x3', '*', '-', '-'], ['x3', 'x2', 'ln', '52.52800399877734', '*', '55.600144949819', 'U-', 'exp', '19.869686257439803', '+', '+', '^', 'x3', '/', '80.87972005504302', '*', '36.932352341528784', 'ln', '/', 'ln'], ['x3', 'x3', 'x2', '+', '83.34434751518123', 'ln', '-', '*'], ['x2', 'U-', 'x3', '*', 'U-', 'x2', '23.227835634616298', 'x1', 'ln', '-', 'ln', '77.35883461525425', '*', 'U-', '+', '-'], ['x2', 'ln', 'exp', '74.21715491394288', 'ln', 'U-', 'U-', '/', 'x3', 'U-', '/', 'U-', 'exp', '54.90811787090503', 'x1', '+', 'x3', '+', '*']]
#best_individual=[['x1', '0', '*', '1', '+'], ['x1', '51.502826512557455', '-', 'x2', '-', 'x2', 'x3', '*', 'x3', '-', 'x2', '-', '+'], ['x3', 'exp', 'x3', '41.49419542959295', '/', '+', '24.46931169083946', '^', 'ln', 'x2', '+', '43.90148299209225', '19.56679294277791', '34.79279884021608', '/', 'ln', 'exp', 'x3', '*', '+', '+'], ['x2', '30.155025331759482', 'x3', '/', '-', 'x2', '^', 'ln', '21.73915889169582', '+', '23.78306820283264', 'x3', '*', '19.5236968544625', '49.47397847297518', 'exp', '/', '-', '+'], ['x1', 'x2', '20.62985581149512', 'ln', '-', 'x1', 'x3', '/', '51.52540780991984', '-', '*', '17.152603613037545', '36.7616060209905', 'U-', '/', '*', '+'], ['x1', 'x3', 'x2', '+', 'x1', 'x3', '36.29935888608503', 'x2', 'ln', '*', '*', '/', '/', '+']]
#best_individual=[['x1', '0', '*', '1', '+'], ['x2', 'ln', 'x3', 'x2', 'ln', '/', '/', 'exp', 'x1', '*', '37.00260157881856', '27.185734076097965', 'ln', '*', 'x1', '/', 'x3', '*', '+'], ['x1', 'U-', 'U-', 'x2', 'ln', '33.37314277771566', 'x2', '*', '51.28314504037124', '49.00883222438101', 'U-', 'U-', '/', 'U-', 'exp', '*', '+', '+'], ['x2', 'ln', '35.16470565461305', '50.737778462898014', '+', '*', 'ln', '35.16470565461305', '50.543466952067334', '+', '*', 'x2', '+'], ['x2', '38.27296644908639', '41.588458426051716', '24.788791486963884', '/', '^', 'x2', '+', '+']]
#best_individual=[['x1', '0', '*', '1', '+'], ['x2', 'ln', '16.496870590469655', '^', '25.806966178948077', '*', 'ln', '24.83870364443528', 'ln', 'exp', '*', '42.53343078839217', 'x2', 'U-', '20.534632963670266', '-', '+', '-'], ['x1', '45.730310596707554', 'x2', 'x3', '44.03059308906431', '47.780514486933725', '-', '^', '-', 'x3', '43.975362390553386', '47.780514486933725', '/', '^', '*', '-', '-'], ['x1', 'x2', '26.790516771686644', '38.661724633756506', 'ln', 'x3', '26.741471900338237', 'x3', '^', '/', 'exp', '/', '+', 'x1', '/', '/', '+']]
#best_individual=[['x1', '0', '*', '1', '+'], ['x2', '44.65120370241437', 'x1', '+', '34.59216915889572', '-', 'x2', '28.952408640121167', '+', '+', '34.594778614880774', '/', 'x2', '30.121562888695983', '+', '*', '+'], ['x2', 'x1', 'x1', '+', 'ln', '*', '40.16780822891951', '+', '27.92754861610171', 'x1', 'x3', 'x1', '/', '+', '+', '+']]
#best_individual=[['x1', '0', '*', '1', '+'], ['x1', '35.534221622469204', '+', 'x1', '-', 'ln', 'ln', 'ln', 'U-', 'x1', '+', 'x2', '18.5274502337934', '-', 'x1', 'ln', '*', '+', 'x2', '18.517247688104522', '+', 'x1', 'ln', '*', '+']]
#best_individual=[['x1', '0', '*', '1', '+'], ['x2', '42.142524905734575', 'x2', '-', 'x2', 'x3', '+', '/', '/', 'x2', 'exp', 'ln', '*', 'x3', '*', '57.18270006537522', '^', 'ln'], ['x2', '77.67920432700937', 'x2', '^', 'x2', '27.128674312593223', '/', '^', 'ln', '+', 'x1', '59.57799442232747', 'x2', 'x2', '+', 'x2', '+', 'x2', '+', '+', '+', '+'], ['x1', 'x2', '-', 'x2', '-', '61.35070038654061', 'x2', 'exp', 'ln', 'x3', '+', '-', '+', 'x2', '-', '60.797946684867675', 'x2', 'exp', 'ln', 'x3', '*', '-', '-'], ['x2', 'x3', '*', 'x2', '42.54418067584284', '51.39178492369522', '+', '-', '46.128811217159075', '-', '46.12308497620506', '50.09102621952748', '+', '-', '48.69675550805538', 'x2', '+', '+', '46.070918383381', '-', '-'], ['x2', '13.59492120576958', '-', '82.97371329243872', 'U-', 'x2', 'ln', '*', '-', '32.102022076593045', 'x2', '-', '-', '32.10370545105927', 'x1', '-', '-']]
##best_individual=[['x1', '0', '*', '1', '+'], ['x3', 'x2', '*', 'x3', '+', 'ln', '69.6213081851031', 'x2', '+', '*'], ['x2', 'U-', '16.532360763401147', 'ln', 'x2', '*', '78.56226877168623', 'ln', '*', '+', '16.532360763401147', 'ln', 'x2', '/', '78.54558914321603', 'ln', '*', '+', 'x1', '17.032446818332453', '29.097270172742284', '/', '+', '+'], ['x3', 'x3', '84.27272270151613', '25.878844319259475', '-', '76.41728331987666', '/', '15.298664440834823', '17.512798632200095', '+', '*', 'ln', 'U-', 'x2', '+', '/', '42.66944270658032', '31.685241515817165', 'x2', '-', '-', '+', '*'], ['x3', 'x2', 'x3', '+', '76.57257232607094', '16.43363228691422', '/', '-', '*']]
##best_individual=[['x1', '0', '*', '1', '+'], ['x1', '53.14606458901729', '54.87121182106304', '+', 'x1', '-', '54.87027473031292', '+', 'x2', '-', '-', 'x2', '-', 'x2', '58.33095743341281', '39.21310046889583', '/', '^', '+'], ['x1', 'x2', 'x2', '^', 'ln', '+', '36.935489913133324', '+', 'x2', '23.40082662226554', 'ln', '*', '+', 'x2', '23.391208102564352', 'ln', '*', '+'], ['x2', '84.73169218095624', 'x2', '64.81606325260263', 'x1', '-', '-', '+', '49.91503145207817', '^', 'ln', '+', '84.71736242522931', 'x2', '+', 'ln', '+', '84.71736242522931', 'x2', '64.37086710625925', 'x3', '*', '/', '^', '^']]
##best_individual=[['x1', '0', '*', '1', '+'], ['x3', 'ln', 'exp', 'ln', 'x2', '-', 'U-', 'ln', 'x2', '-', 'U-', '12.317427166200886', '*', 'x1', 'ln', '42.43221909266771', 'ln', 'U-', '81.87187249918186', '*', '-', '+'], ['x2', '69.00143826711593', 'x2', '+', '^', 'ln', '73.0823504168323', '31.402414396980696', '+', '+', '74.51619425071395', '+', '73.94576158809905', '50.41311783792886', '-', '46.88498819158077', '+', '50.41311783792886', '-', '46.91906690806133', '+', '+']]
##best_individual=[['x1', '0', '*', '1', '+'], ['x2', 'x2', '*', '84.68833026133244', 'exp', 'ln', '^', '84.87159466609211', 'ln', 'ln', '-', 'ln', 'x2', '33.99935420759307', '-', '-']]
##best_individual=[['x1', '0', '*', '1', '+'], ['x2', 'U-', '32.25363557635897', 'x2', '*', 'ln', 'x2', '65.87361908081148', '+', '*', '+'], ['x1', 'x2', 'x1', '36.899179309682', '/', '/', '30.809441174165354', '61.555916906278675', '+', '*', '+'], ['x3', 'x2', '52.82408289782073', '31.56921017331655', '67.80501675339116', '+', '/', '*', '52.75215115295083', 'ln', 'ln', '19.121709252060583', '+', '+', '*'], ['x2', 'ln', 'U-', 'x1', '/', 'x2', '24.972645696424454', '*', 'x2', '*', 'ln', '51.01839579928352', '*', '+', 'x2', 'ln', 'ln', '50.883016811151236', '*', '+'], ['x1', 'x2', 'x1', 'U-', '/', '59.09500481323676', '*', '-', 'x2', 'x1', 'ln', '/', '59.09500481323676', '*', '+']]
#best_individual=[['x1', '0', '*', '1', '+'], ['x3', '38.49994875967064', 'x2', 'x2', 'ln', 'ln', '^', '^', '+', 'ln', 'x3', '-', '49.40258022395258', '51.13404009972736', '/', '-', 'x3', '25.880974329201745', '*', '+'], ['x3', 'ln', 'U-', 'x2', '*', 'U-', 'U-', '44.225252366338474', 'x1', '63.174904657659866', '+', '+', 'U-', 'exp', 'x2', '19.52523288298337', '+', '14.901747568622362', '*', '-', '-'], ['x3', 'ln', 'x2', '-', '13.20632230331481', 'x2', '*', '+', 'x1', 'x2', '-', 'U-', '-'], ['x2', 'x3', 'ln', '+', 'ln', 'x3', 'ln', '54.31855984078441', '*', 'x2', '+', '*']]
#best_individual=[['x1', '0', '*', '1', '+'], ['x2', 'ln', '17.202591412663903', '-', 'U-', 'x2', '*', 'ln', 'exp', '60.384810553126776', 'x1', '-', '-'], ['x2', '35.8656610396232', 'x2', 'exp', 'ln', 'ln', '/', '+', 'x3', '*', '20.58613320583748', '-'], ['x3', 'U-', 'x2', '*', 'U-', '66.98802135046904', '29.774532747221386', 'ln', '44.74231629224278', '*', '+', '+']]
best_individual=[['x1', '0', '*', '1', '+'], ['x3', '63.52900190673749', '^', 'ln', 'x2', 'ln', '^', '34.35854518276749', 'U-', 'U-', '^', 'ln'], ['x3', 'x1', 'x2', '/', '+', '21.851440223478246', '+', 'x3', 'U-', 'x2', '*', 'U-', '67.45153988642079', '+', '63.971988376838915', 'ln', 'exp', '+', '+']]
best_individual=[['x1', '0', '*', '1', '+'], ['x2', '32.13290525545109', '-', 'x2', '+', '27.910143644451956', 'ln', '*', 'x3', 'exp', 'ln', 'ln', 'ln', '+', '27.923017250201312', 'x3', '*', '+']]
best_individual=[['x1', '0', '*', '1', '+'], ['x1', 'ln', 'ln', '49.65033140123362', 'x1', 'U-', '-', 'U-', 'x3', 'ln', '56.16773459677237', 'U-', '+', '+', '-'], ['x1', '51.2484192918807', '54.88496795758748', '+', '+', '44.30524540913973', 'x3', 'U-', '14.715818192904546', '-', '66.90371315178447', '50.87699949281866', '/', '*', '+', '+'], ['x1', 'ln', 'x1', '37.61171023411113', '+', '29.633386481878023', '+', '29.657198800287595', '+', '72.39467856229722', '32.30937920668134', 'x3', '/', '^', 'x2', '/', 'ln', '+', '+'], ['x1', 'x3', '/', '22.219844158032835', '*'], ['x1', 'x3', '78.12881515986635', 'x2', '+', '-', '-', 'x1', '5.253467030482258', '*', '63.55739221437503', 'exp', '66.65241581347959', 'x3', 'ln', '61.72174397991985', '^', '*', '/', '*', '+']]
best_individual=[['x1', '0', '*', '1', '+'], ['x1', 'ln', 'ln', '22.68969809730428', 'x1', '31.194612604120774', '+', '30.73526220776818', '+', 'x3', '18.192050742995246', '+', '/', '*', '*'], ['x1', '27.811584396610016', '+', '18.660874640955488', '+', 'x1', 'ln', '+', '21.756717559404795', '+', '45.87544581447075', 'x3', '-', '+'], ['x1', '31.929632024269125', 'x2', 'U-', 'exp', '39.96648982814817', '+', '+', '31.824651248705475', '+', '+'], ['x1', '43.74018881594536', '34.09817420776413', '44.85717646809658', 'x1', '31.94241899528769', 'x3', '*', '/', '*', '+', '+', '+'], ['x1', '46.31450465620024', '19.738974882837553', '+', '+', 'x3', '-', '42.31214386679903', 'U-', 'U-', 'ln', 'exp', '27.80561698606096', '38.40647229429297', '/', '18.7024494176969', '*', '+', '+']]
#
best_individual=[['x1', '0', '*', '1', '+'], ['x1', '53.869056677372555', 'x3', '57.210528629057045', '-', '-', '+'], ['x3', 'U-', '52.923868413768275', '79.70342661764268', '40.17775469058066', '/', 'U-', '62.4636495366248', 'x1', '+', '-', '-', '+'], ['x1', '63.7772846764612', '+', 'x3', 'x2', 'U-', 'x2', '+', '+', 'x3', 'U-', '29.679133383508614', '75.84239880234921', '-', '+', '+', '-'], ['x1', '64.83306111933446', 'U-', 'exp', '35.475471968174155', '37.55799142780382', '-', '+', '+', '46.2968118000892', 'x3', '79.59976372874733', 'x3', '30.109089814202306', '/', '-', '-', '-', '+'], ['x1', 'x3', '/', 'x1', 'x1', '73.62471408614078', 'U-', '+', '+', 'ln', 'x1', '73.60031507617907', '-', '*', 'ln', 'x1', '74.9440578461595', 'U-', '-', '+', '+']]
#
best_individual=[['x1', '0', '*', '1', '+'], ['x3', 'ln', 'ln', '25.015575358216296', '21.733925680813318', '*', '^', 'ln', 'U-', '18.57310171665837', '21.733925680813318', '*', '+'], ['x1', '21.218237835096627', '19.25493635455135', 'exp', 'ln', '14.930626131667005', '-', '*', '+'], ['x1', '19.561562811608592', 'x3', '/', '29.34538845037721', 'x1', '*', 'ln', 'U-', 'x1', 'U-', '-', 'ln', 'U-', '15.346169813126139', '*', '*', '-'], ['x1', '16.353430643917573', 'ln', '26.697156095014975', 'ln', '+', 'x3', '20.97097793384778', '-', '/', '26.697156095014975', 'ln', '30.23411754538501', '*', '+', '+'], ['x1', 'x3', '-', 'x1', '16.77695314407741', '^', 'ln', 'ln', '24.124679137709805', '21.4046872292083', '16.6387381669784', '-', '*', '+', '+']]
#################################
best_individual=[['x1', '0', '*', '1', '+'], ['x2', '10.36122859866811', '^', 'ln', '30.01519398453371', '9.276620193151112', '6.169105080453104', '5.98278953661611', '/', '^', '+', '9.46574314013868', 'ln', '/', '*'], ['x2', 'x2', '*', 'x3', 'x1', 'x2', '-', '+', '+', '3.006759855578316', '/', 'x2', '+', 'x3', 'x1', 'x3', '-', '+', '+'], ['x3', 'x2', '11.900904467652262', 'x3', 'x2', 'U-', 'x3', '8.78382284382298', 'U-', '+', '-', '/', '-', '+', '*'], ['x2', '26.086420178758594', 'x2', 'x2', '17.614515265541055', '/', '/', '*', 'x3', 'x2', '+', '+', '+'], ['x1', '15.554473826607158', 'x3', 'x3', 'U-', '^', '+', 'ln', 'U-', 'x2', '16.643452959516885', 'ln', 'ln', 'x3', '27.174308189870047', '-', '-', '*', '+', '+']]
best_individual=[['x1', '0', '*', '1', '+'], ['x3', 'ln', 'ln', 'x1', '^', 'ln', 'x3', '*', 'x3', 'x2', 'x3', '*', '-', '-'], ['x3', 'x2', 'ln', '52.52800399877734', '*', '55.600144949819', 'U-', 'exp', '19.869686257439803', '+', '+', '^', 'x3', '/', '80.87972005504302', '*', '36.932352341528784', 'ln', '/', 'ln'], ['x3', 'x3', 'x2', '+', '83.34434751518123', 'ln', '-', '*'], ['x2', 'U-', 'x3', '*', 'U-', 'x2', '23.227835634616298', 'x1', 'ln', '-', 'ln', '77.35883461525425', '*', 'U-', '+', '-'], ['x2', 'ln', 'exp', '74.21715491394288', 'ln', 'U-', 'U-', '/', 'x3', 'U-', '/', 'U-', 'exp', '54.90811787090503', 'x1', '+', 'x3', '+', '*']]


best_individual=[['x1', '0', '*', '1', '+'], ['x1', 'x2', 'ln', 'U-', '+', 'x1', 'x3', '/', '+', '43.53285154185595', 'x1', '/', '59.1200558403478', '/', '78.3926314920336', '^', '86.71874802447266', '+', '+'], ['x1', 'U-', 'U-', '30.24286680140764', 'ln', '87.77226355452355', 'ln', '*', '88.10230500690912', '+', '+'], ['x1', '70.25422759632197', '32.059757965196006', '64.2391856915819', '-', '-', '+'], ['x1', '57.772447524334396', '23.388137630896075', '58.9630490832074', 'x3', '-', '/', '/', '+'], ['x2', '25.17729773498366', '-', 'U-', 'U-', 'U-', '33.470699169953235', 'x1', '56.82248782965497', 'ln', 'x1', '-', 'exp', '70.26873802969286', 'x3', '-', '+', '+', '+', '+']]
best_individual=[['x1', '0', '*', '1', '+'], ['x1', '74.83135291465287', '+', '29.589425626152273', '+', 'x3', 'x3', 'x3', 'ln', 'exp', 'U-', '41.95779742617839', 'U-', 'x3', '/', '-', '/', '^', '+'], ['x1', '49.417405171094785', 'x3', '-', '38.43429649663258', 'x1', 'U-', '/', '+', '77.76382821983623', '+', '+'], ['x1', '69.03016254295039', '62.49232754847746', '31.282674876990157', '78.64308026368977', '68.85170590748181', 'U-', '47.599849623746564', 'U-', '*', 'x3', 'ln', 'exp', '/', '*', '-', '+', '-', '+', '37.663801251480486', '/'], ['x1', 'ln', 'exp', 'ln', 'exp', '82.31121259383303', '+', 'x3', 'U-', '+', 'x3', 'U-', 'x3', '-', '+', 'x2', 'x3', '-', '66.60873394843192', '+', '+', 'x3', 'U-', 'x3', '-', '+', 'x3', 'ln', 'x3', '-', '66.54248179710572', '+', '+'], ['x1', '34.50224736008452', '61.17642994874919', '52.32651207548697', 'ln', 'x3', 'U-', '-', '/', '*', '+', 'x2', 'x2', 'x1', '36.67649716456822', 'x2', '45.539274515814526', '*', '36.21563614159538', '/', '-', '*', '/', '^', '^'], ['x1', 'ln', 'ln', 'x1', '80.27299035426496', 'x3', '/', '52.54086512469135', '30.278408157157585', '-', '*', '+', '+'], ['x1', '85.93418917101768', 'x3', 'U-', 'U-', 'x2', 'ln', '*', '74.04943039093757', '/', '29.4620482137732', '+', 'U-', 'x3', '+', '-', 'x3', 'U-', '29.597834196433773', '+', 'U-', 'x3', '+', '-', 'x2', '74.04943039093757', '/', '29.597834196433773', '+', 'U-', 'x3', '+', '-', '+']]
best_individual=[['x1', '0', '*', '1', '+'], ['x1', 'x3', '-', '56.618123111120354', 'U-', '-', '54.92298162728743', '-', 'x3', '-', '56.312015762259136', 'U-', '-', 'x3', '-', 'x3', '-', '56.312015762259136', 'U-', '-', '55.28999540012778', 'x1', '58.00041165085283', 'U-', 'U-', 'x1', '/', 'x3', '/', 'x1', '+', '/', '+', '+'], ['x1', '39.63007544193186', '+', '39.64761516644167', 'x3', 'x3', 'U-', '30.435888410990895', '+', '-', 'x3', 'U-', '30.435888410990895', '+', '-', 'x3', 'ln', '*', 'x3', 'U-', '30.436192984087523', '/', '-', '-', '+'], ['x1', '28.408242348456017', 'x3', '-', '+', '28.02440745252913', 'x3', '-', '+', '28.02440745252913', 'x3', '-', 'x3', '-', '+', '28.02440745252913', 'x3', '-', 'x3', '-', '+', '28.04511995685897', 'x3', '-', '+', '28.054514650041288', 'U-', '52.04522690452899', '58.94467321719355', 'x3', '/', '+', '-', '-'], ['x1', '36.804121586955794', 'x3', '-', 'x3', '-', '+', '37.74748398162145', 'x3', '-', 'x3', '-', '+', '37.74748398162145', 'x3', '-', 'x3', '-', '+', '55.62444210578373', 'x3', '-', 'x3', '-', 'x3', '-', '+', '55.62444210578373', '33.47847636862859', '+', '+'], ['x1', '51.324552184921146', '+', '51.324552184921146', '+', 'ln', 'exp', '45.92192255356877', '+', '43.65905560092608', '+', 'ln', 'exp', '45.34193022841377', 'U-', '12.28406113129995', '+', '+', 'x2', '12.28406113129995', '+', '+', 'x3', '/', '42.104998659812594', 'ln', 'ln', 'exp', '*', 'x1', '+'], ['x1', '43.3867359481837', '54.65573902985321', 'x1', '+', 'x1', '19.585988201324206', 'x3', '-', '-', '-', '/', '/'], ['x3', '23.269406787424657', 'x1', '30.772403973035992', '+', '+', '25.438140786230083', '-', '-', 'U-', '56.1001507681452', 'x3', '3.818647477179092', '+', '-', '15.39692894165585', '40.85099193764191', '+', '+', '+']]


for i in range(len(best_individual)):
    print(((rss(eqn_to(get_eqn(best_individual[i])))/len(y))/mean)*100, '%')
    #best_individual[i]=clear_unar_minus(best_individual[i])
X=np.array(model_of_best_individual(best_individual))
num_eqn=len(best_individual)-1

y1,bb,SS2,r2=mnk(X,y)
#print(bb,'Коэф. детерминации: ',R2)
#print('PREDICT: ', y1)
string=''
for i in range(num_eqn+1):
    #print('')
    #print(best_individual[i])
    #print(get_eqn(best_individual[i]))
    if i!=0:
        string=string+'+'+str(bb[i])+'*'+'('+get_eqn(best_individual[i])+')'
    else:
        string=string+str(bb[i])+'*'+'('+get_eqn(best_individual[i])+')'
string=string.replace('+-','-')
string=string.replace('-+','-')
string=string.replace('--','+')
string=string.replace('++','+')
#eqn=simplify(string)
eqn=string
print('Eqn: ',eqn)
print('RSS: ',rss(eqn_to(eqn)))






F01='34.80119833847266*(((x1*0)+1))-0.9025224168624104*((x2+(((((44.65120370241437+x1)-34.59216915889572)+(x2+28.952408640121167))/34.594778614880774)*(x2+30.121562888695983))))+1.6899765235800948*((((x2*ln((x1+x1)))+40.16780822891951)+(27.92754861610171+(x1+(x3/x1)))))'
# 0.75
F02='267.2372608750597*(((x1*0)+1))+0.35350275236739637*(((((-ln(ln(ln(((x1+35.534221622469204)-x1)))))+x1)+((x2-18.5274502337934)*ln(x1)))+((x2+18.517247688104522)*ln(x1))))'
#0.19
F03='246.07882115091058*(((x1*0)+1))-0.3152078581760285*((ln(((x3*x2)+x3))*(69.6213081851031+x2)))+0.8382654556700235*(((((-x2)+((ln(16.532360763401147)*x2)*ln(78.56226877168623)))+((ln(16.532360763401147)/x2)*ln(78.54558914321603)))+(x1+(17.032446818332453/29.097270172742284))))-0.0006951221053562646*((x3*((x3/((-ln((((84.27272270151613-25.878844319259475)/76.41728331987666)*(15.298664440834823+17.512798632200095))))+x2))+(42.66944270658032-(31.685241515817165-x2)))))-0.14489098266406802*((x3*((x2+x3)-(76.57257232607094/16.43363228691422))))'
#0.93
F04='-262.6382598228479*(((x1*0)+1))+2.943617636203686*((((x1-((((53.14606458901729+54.87121182106304)-x1)+54.87027473031292)-x2))-x2)+(x2**(58.33095743341281/39.21310046889583))))-5.825589212892279*(((((x1+ln((x2**x2)))+36.935489913133324)+(x2*ln(23.40082662226554)))+(x2*ln(23.391208102564352))))+4.576344885908382*((((x2+ln(((84.73169218095624+(x2-(64.81606325260263-x1)))**49.91503145207817)))+ln((84.71736242522931+x2)))**(84.71736242522931**(x2/(64.37086710625925*x3)))))'
# 0.71
F05='532.6197581505184*(((x1*0)+1))-0.3101651564308274*((((((x2-32.13290525545109)+x2)*ln(27.910143644451956))+ln(ln(ln(exp(x3)))))+(27.923017250201312*x3)))'
# 0.60
F06='237.7039530052938*(((x1*0)+1))-0.6453635010544403*((ln(ln(x1))-((-(49.65033140123362-(-x1)))+(ln(x3)+(-56.16773459677237)))))-7.312398655542371*(((x1+(51.2484192918807+54.88496795758748))+(44.30524540913973+(((-x3)-14.715818192904546)*(66.90371315178447/50.87699949281866)))))+6.877550224047402*((ln(x1)+((((x1+37.61171023411113)+29.633386481878023)+29.657198800287595)+ln(((72.39467856229722**(32.30937920668134/x3))/x2)))))+1.2253343387754947*(((x1/x3)*22.219844158032835))+0.3293002604212103*(((x1-(x3-(78.12881515986635+x2)))+((x1*5.253467030482258)*(exp(63.55739221437503)/(66.65241581347959*(ln(x3)**61.72174397991985))))))'
# 0.92
F07='1490.81095600745*(((x1*0)+1))+1.6568731033649264*((ln(ln(x1))*(22.68969809730428*(((x1+31.194612604120774)+30.73526220776818)/(x3+18.192050742995246)))))+233.13284390315278*((((((x1+27.811584396610016)+18.660874640955488)+ln(x1))+21.756717559404795)+(45.87544581447075-x3)))-0.45167705487261856*((x1+((31.929632024269125+(exp((-x2))+39.96648982814817))+31.824651248705475)))+13.443203383708465*((x1+(43.74018881594536+(34.09817420776413+(44.85717646809658*(x1/(31.94241899528769*x3)))))))-249.17059120064775*((((x1+(46.31450465620024+19.738974882837553))-x3)+(exp(ln((-(-42.31214386679903))))+((27.80561698606096/38.40647229429297)*18.7024494176969))))'
# 0.92
F08='1579.5556582602076*(((x1*0)+1))-3.211448379668248*(((-ln((ln(ln(x3))**(25.015575358216296*21.733925680813318))))+(18.57310171665837*21.733925680813318)))-3.354298579879014*((x1+(21.218237835096627*(ln(exp(19.25493635455135))-14.930626131667005))))+11.860902794611013*((x1-((19.561562811608592/x3)*((-ln(((-ln((29.34538845037721*x1)))-(-x1))))*15.346169813126139))))+0.06027828475190597*((x1+(((ln(16.353430643917573)+ln(26.697156095014975))/(x3-20.97097793384778))+(ln(26.697156095014975)*30.23411754538501))))-8.439580591295083*(((x1-x3)+(ln(ln((x1**16.77695314407741)))+(24.124679137709805*(21.4046872292083-16.6387381669784)))))'
# 0.90

F00=F01+'+'+F02+'+'+F03+'+'+F04
#0.94
F11='(12766.08809473073*x3*(0.02452844676346446*x3 - 1)**2*log(x1) + x3*(x3 + 13.616379314590645)*(1.0377659484148584*x1 + 0.12561395998836633*x2 + 6.014813501589664*x3 - 7.680649249946354*log(x2) + 297.16648019197761) - (x3 + 13.616379314590645)*(11.653096181583347*x1 + 7.128509059131217*x3 + 11.77871014157171*log(exp(-x2)) + 15083.280558808076))/(x3*(x3 + 13.616379314590645))'
# 0.93
F12='(-1951.8455678195908*x1*x3**2*(x2 + 21.595745172833293)*(x3 - 18.07721252328992) - 35267.858660730633*x1*x3*(x2 + 21.595745172833293)*(x3 + 5.307917447661647e-9) + x3**2*(x2 + 21.595745172833293)*(x3 + 5.307917447661647e-9)*(1951.7020869341984*x1 - 0.88052256295512602*x2 + 0.74753638455399786*x3 + 0.5256458586448435*(x3*log(x1) - 58.755339339877594)*log(log(x1)) - 0.24917879485133262*log(x1) + 1508.0834390158052*log(log(log(75.72391298594927 - log(x2)))) - 517.10446572595072) - x3**2*(11.716092041982779*x3 + 6.2188049368050164e-8) - 787.691005410659*(x2 + 21.595745172833293)*(0.02831098670534173*x3 - 1)**2*(x3 - 117.89336718763058)*(x3 + 5.307917447661647e-9))/(x3**2*(x2 + 21.595745172833293)*(x3 + 5.307917447661647e-9))'
# 0.95 (dolgo rabotaet)
F13='21.882905564928706*(((x1*0)+1))+0.06762763946658944*((x1+(((49.31317267416327+51.3916936859223)+(27.983081068869325/(x3-23.525180756066618)))+(27.987658529337242-x3))))-3.841922019669763*((x1+((52.948151094864684+50.5995605673569)-(ln((((((exp(x3)/28.609349099135674)/28.627902033038946)/28.627902033038946)/51.16072833374824)/50.77996690367783))*ln(x3)))))-2.0255828137567473*((x1+((-(x3+15.052389218031264)+x3)/(55.92580849948869*(-x3/(x1*19.818629072490932))))))+4.983259067987561*((x1+(-x3-(-(35.16669609715807+15.536188408166195)-ln((ln((ln((ln((51.93927542996871**15.85850803623748))**16.524696064920022))**16.524696064920022))**16.524696064920022))))))+1.7443884664280964*((ln((exp(x1)*ln(16.11434556204978)))*(22.94013222190214*(36.381414058586/((x2+(36.02847474944257*x3))+(36.118615230043375/x3))))))'
# 0.91

F14='-15341.026260329794*(((x1*0)+1))+0.05674898527195521*((ln((x2**10.36122859866811))*((30.01519398453371+(9.276620193151112**(6.169105080453104/5.98278953661611)))/ln(9.46574314013868))))+2.9025758935441903*((((((x2*x2)+(x3+(x1-x2)))/3.006759855578316)+x2)+(x3+(x1-x3))))-3.289801820507482*((x3*(x2+(11.900904467652262-(x3/((-x2)-(x3+(-8.78382284382298))))))))+34.11097810683487*((x2+((26.086420178758594*(x2/(x2/17.614515265541055)))+(x3+x2))))-3.0097949433359616*((x1+((-ln((15.554473826607158+(x3**(-x3)))))+(x2*(ln(ln(16.643452959516885))-(x3-27.174308189870047))))))'
# 0.78
F15='(12766.08809473073*x3*(0.02452844676346446*x3 - 1)**2*log(x1) + x3*(x3 + 13.616379314590645)*(1.0377659484148584*x1 + 0.12561395998836633*x2 + 6.014813501589664*x3 - 7.680649249946354*log(x2) + 297.16648019197761) - (x3 + 13.616379314590645)*(11.653096181583347*x1 + 7.128509059131217*x3 + 11.77871014157171*log(exp(-x2)) + 15083.280558808076))/(x3*(x3 + 13.616379314590645))'
# 0.93
F16='21.882905564928706*(((x1*0)+1))+0.06762763946658944*((x1+(((49.31317267416327+51.3916936859223)+(27.983081068869325/(x3-23.525180756066618)))+(27.987658529337242-x3))))-3.841922019669763*((x1+((52.948151094864684+50.5995605673569)-(ln((((((exp(x3)/28.609349099135674)/28.627902033038946)/28.627902033038946)/51.16072833374824)/50.77996690367783))*ln(x3)))))-2.0255828137567473*((x1+((-(x3+15.052389218031264)+x3)/(55.92580849948869*(-x3/(x1*19.818629072490932))))))+4.983259067987561*((x1+(-x3-(-(35.16669609715807+15.536188408166195)-ln((ln((ln((ln((51.93927542996871**15.85850803623748))**16.524696064920022))**16.524696064920022))**16.524696064920022))))))+1.7443884664280964*((ln((exp(x1)*ln(16.11434556204978)))*(22.94013222190214*(36.381414058586/((x2+(36.02847474944257*x3))+(36.118615230043375/x3))))))'
# 0.91
F17='(-3.3254470986936413*x3*(x3 + (x2 + 11.900904467652262)*(x2 + x3 - 8.78382284382298)) + (x2 + x3 - 8.78382284382298)*(-0.15347066826492397*x1 + 0.72146229742871914*x2**2 + 3.0441968389156187*x2*(x3 - 28.208210182097752) + 67.306895914930717*x2 + 33.651009466997449*x3 + 13.27064938546032*log(x2**10.36122859866811) + 3.0441968389156187*log(x3**(-x3)*(15.554473826607158*x3**x3 + 1)) + 454.37523109396945))/(x2 + x3 - 8.78382284382298)'
#0.93
F18='((-14676.740051478639*((x1*0)+1))+((0.746484836823839*(ln((x2**10.36122859866811))*((30.01519398453371+(9.276620193151112**(6.169105080453104/5.98278953661611)))/ln(9.46574314013868))))+(((2.1692638732219756*(((((x2*x2)+(x3+(x1-x2)))/3.006759855578316)+x2)+(x3+(x1-x3))))-(3.3254470986936413*(x3*(x2+(11.900904467652262-(x3/(-x2-(x3-8.78382284382298))))))))+((32.92954716956873*(x2+((26.086420178758594*(x2/(x2/17.614515265541055)))+(x3+x2))))-(3.0441968389156187*(x1+(-ln((15.554473826607158+(x3**-x3)))+(x2*(ln(ln(16.643452959516885))-(x3-27.174308189870047))))))))))'
# 0.78
F19='1096.166975517887*(((x1*0)+1))+7.934065316507343*(((ln((ln(ln(x3))**x1))*x3)-(x3-(x2*x3))))-1.3844405590590343*(ln(((((x3**((ln(x2)*52.52800399877734)+(exp((-55.600144949819))+19.869686257439803)))/x3)*80.87972005504302)/ln(36.932352341528784))))-33.256084578875324*((x3*((x3+x2)-ln(83.34434751518123))))+30.007129741775742*(((-((-x2)*x3))-(x2+(-(ln((23.227835634616298-ln(x1)))*77.35883461525425)))))-3.9171904684203125*((exp((-((exp(ln(x2))/(-(-ln(74.21715491394288))))/(-x3))))*((54.90811787090503+x1)+x3)))'
#0.94
F20='0.4814301843903332*((x1*exp(((ln(x2)**2)/x3))))-0.22619260114866666*(x1)-3.22352061991986*(x2)-7.831098729966989*(ln(x2))-2.993231481961061*(ln(x2))+474.22104595501986*(1)-2494.033880360498*((x3/x1))'
#0.74
F21='1013.896242928098*(((x1*0)+1))+3.314980777573236*((x2+(((((44.65120370241437+x1)-34.59216915889572)+(x2+28.952408640121167))/34.594778614880774)*(x2+30.121562888695983))))-4.590337292476218*((((x2*ln((x1+x1)))+40.16780822891951)+(27.92754861610171+(x1+(x3/x1)))))'
#0.75
F22='48.05230120417268*(((x1*0)+1))+0.9176332844536984*(((((-ln(ln(ln(((x1+35.534221622469204)-x1)))))+x1)+((x2-18.5274502337934)*ln(x1)))+((x2+18.517247688104522)*ln(x1))))'
# 0.2
F23='252.69358071540265*(((x1*0)+1))+5.094250515621791*(ln(((((x2/((42.142524905734575-x2)/(x2+x3)))*ln(exp(x2)))*x3)**57.18270006537522)))+2.8238553132738264*(((x2+ln(((77.67920432700937**x2)**(x2/27.128674312593223))))+(x1+(59.57799442232747+(((x2+x2)+x2)+x2)))))+2.8526907907472028*((((((x1-x2)-x2)+(61.35070038654061-(ln(exp(x2))+x3)))-x2)-(60.797946684867675-(ln(exp(x2))*x3))))-4.697244954454292*(((x2*x3)-(((((x2-(42.54418067584284+51.39178492369522))-46.128811217159075)-(46.12308497620506+50.09102621952748))+(48.69675550805538+x2))-46.070918383381)))-5.835571082857385*(((((x2-13.59492120576958)-((-82.97371329243872)*ln(x2)))-(32.102022076593045-x2))-(32.10370545105927-x1)))'
# 0.76
F24='13.322922401413848*(((x1*0)+1))+0.9157763786640558*((ln(((x3*x2)+x3))*(69.6213081851031+x2)))-0.47525705368457183*(((((-x2)+((ln(16.532360763401147)*x2)*ln(78.56226877168623)))+((ln(16.532360763401147)/x2)*ln(78.54558914321603)))+(x1+(17.032446818332453/29.097270172742284))))+0.7225146534247503*((x3*((x3/((-ln((((84.27272270151613-25.878844319259475)/76.41728331987666)*(15.298664440834823+17.512798632200095))))+x2))+(42.66944270658032-(31.685241515817165-x2)))))-0.2024827081800567*((x3*((x2+x3)-(76.57257232607094/16.43363228691422))))'
# 0.54
F25='161.04945429747607*(((x1*0)+1))+0.13603266886082865*((((x1-((((53.14606458901729+54.87121182106304)-x1)+54.87027473031292)-x2))-x2)+(x2**(58.33095743341281/39.21310046889583))))-0.8145000936881159*(((((x1+ln((x2**x2)))+36.935489913133324)+(x2*ln(23.40082662226554)))+(x2*ln(23.391208102564352))))+1.3868257067425482*((((x2+ln(((84.73169218095624+(x2-(64.81606325260263-x1)))**49.91503145207817)))+ln((84.71736242522931+x2)))**(84.71736242522931**(x2/(64.37086710625925*x3)))))'
# 0.7
F26='372.75554668394034*(((x1*0)+1))-0.7011336160305143*(((-x2)+(ln((32.25363557635897*x2))*(x2+65.87361908081148))))+3.967825106249726*((x1+((x2/(x1/36.899179309682))*(30.809441174165354+61.555916906278675))))-0.37784750643686044*((x3*((x2*(52.82408289782073/(31.56921017331655+67.80501675339116)))+(ln(ln(52.75215115295083))+19.121709252060583))))-1.1727372371071363*(((((-ln(x2))/x1)+(ln(((x2*24.972645696424454)*x2))*51.01839579928352))+(ln(x2)*50.883016811151236)))-1.486148065687141*(((x1-((x2/(-x1))*59.09500481323676))+((x2/ln(x1))*59.09500481323676)))'
# 0.76
F27='116.21196044078272*(((x1*0)+1))-0.40752648041489353*((exp(ln(((-(ln(x2)-17.202591412663903))*x2)))-(60.384810553126776-x1)))+0.39994361909455795*((((x2+(35.8656610396232/ln(ln(exp(x2)))))*x3)-20.58613320583748))+0.742385272853045*(((-((-x3)*x2))+(66.98802135046904+(ln(29.774532747221386)*44.74231629224278))))'
# 0.63
F28='-321.62944364196517*(((x1*0)+1))+1.7031637163816953*((((((x2-32.13290525545109)+x2)*ln(27.910143644451956))+ln(ln(ln(exp(x3)))))+(27.923017250201312*x3)))'
# 0.6


F1=F03+'+'+F11
FF='3*(1/x2)+exp(x2-x3)+ln(ln(x3))+124*x1-4.0525917570105782*x3-9.0453039075642042*x3/(2.1256187385566602*x3-32.906563600063667)+124.86348726910705+286.94834124265074/(2.1256187385566602*x3-32.906563600063667)'

F1='761.1892879289538*(((x1*0)+1))+28.491590456355276*((((x1+(-ln(x2)))+(x1/x3))+((((43.53285154185595/x1)/59.1200558403478)**78.3926314920336)+86.71874802447266)))+4.737927090714189*(((-(-x1))+((ln(30.24286680140764)*ln(87.77226355452355))+88.10230500690912)))-26.792871330444797*((x1+(70.25422759632197-(32.059757965196006-64.2391856915819))))-0.5070978423535231*((x1+(57.772447524334396/(23.388137630896075/(58.9630490832074-x3)))))-6.676382113713274*(((-(-(-(x2-25.17729773498366))))+(33.470699169953235+(x1+(exp((ln(56.82248782965497)-x1))+(70.26873802969286-x3))))))'
#^ tipa De=const

F1=F1+'+'+FF

F1='-0.7751894809470771*x1+27.93838042674635*x1/x3+2.5816742449684256*x2+8.366945277620255*x3+6.220767390619756*x3/(2.1256187385566602*x3-32.906563600063667)+0.1721894915108182*exp(x2)*exp(-x3)-14.039326019218862*ln(x2)-69.59971448417127/(2.1256187385566602*x3-32.906563600063667)'


F1=F1+'+'+'101267.87033142726*(((x1*0)+1))-1084.07269921663*((((x1+74.83135291465287)+29.589425626152273)+(x3**(x3/((-exp(ln(x3)))-((-41.95779742617839)/x3))))))+116.10573117329943*((x1+(((49.417405171094785-x3)+(38.43429649663258/(-x1)))+77.76382821983623)))-260.36624871680306*(((x1+(69.03016254295039-(62.49232754847746+(31.282674876990157-(78.64308026368977*(((-68.85170590748181)*(-47.599849623746564))/exp(ln(x3))))))))/37.663801251480486))-3.8938621124048503*(((((((exp(ln(exp(ln(x1))))+82.31121259383303)+(-x3))+((-x3)-x3))+((x2-x3)+66.60873394843192))+((-x3)-x3))+((ln(x3)-x3)+66.54248179710572)))+5.445458405742189*(((x1+(34.50224736008452*(61.17642994874919/(ln(52.32651207548697)-(-x3)))))**(x2**(x2/(x1*(36.67649716456822-((x2*45.539274515814526)/36.21563614159538)))))))+999.6966351923065*((ln(ln(x1))+(x1+((80.27299035426496/x3)*(52.54086512469135-30.278408157157585)))))-26.18552190009892*((x1+(((85.93418917101768-((-((((-(-x3))*ln(x2))/74.04943039093757)+29.4620482137732))+x3))-((-((-x3)+29.597834196433773))+x3))-((-((x2/74.04943039093757)+29.597834196433773))+x3))))'

F1=F03+'-111637.10851851228*(((x1*0)+1))-141009.38950347103*((((((((((x1-x3)-(-56.618123111120354))-54.92298162728743)-x3)-(-56.312015762259136))-x3)-x3)-(-56.312015762259136))+(55.28999540012778+(x1/((((-(-58.00041165085283))/x1)/x3)+x1)))))-0.9213722905163131*(((x1+39.63007544193186)+(39.64761516644167-((((x3-((-x3)+30.435888410990895))-((-x3)+30.435888410990895))*ln(x3))-((-x3)/30.436192984087523)))))+339.93121747623263*(((((((x1+(28.408242348456017-x3))+(28.02440745252913-x3))+((28.02440745252913-x3)-x3))+((28.02440745252913-x3)-x3))+(28.04511995685897-x3))-((-28.054514650041288)-(52.04522690452899+(58.94467321719355/x3)))))+40041.506748388376*((((((x1+((36.804121586955794-x3)-x3))+((37.74748398162145-x3)-x3))+((37.74748398162145-x3)-x3))+(((55.62444210578373-x3)-x3)-x3))+(55.62444210578373+33.47847636862859)))-5.06731942408984*((((((exp(ln(((exp(ln(((x1+51.324552184921146)+51.324552184921146)))+45.92192255356877)+43.65905560092608)))+((-45.34193022841377)+12.28406113129995))+(x2+12.28406113129995))/x3)*exp(ln(ln(42.104998659812594))))+x1))+2.3109733116170696*((x1/(43.3867359481837/((54.65573902985321+x1)-(x1-(19.585988201324206-x3))))))+100633.98631325498*(((-(x3-((23.269406787424657+(x1+30.772403973035992))-25.438140786230083)))+((56.1001507681452-(x3+3.818647477179092))+(15.39692894165585+40.85099193764191))))'


F1=F1.replace('+-','-')
F1=F1.replace('-+','-')
F1=F1.replace('--','+')

#F1='0.0084959721207776*((x1*x3))+1142988.4536007242*((x1/(x1+(58.00041165085283/(x1*x3)))))+6.833524499549682*((x1/x3))+0.6886691064630173*((x2*x3))+18.780084162921256*((x2*ln(((x2*x3)+x3))))-180.40473456387622*(x2)+851.887103495899*((x2/x3))-0.00015724204646610628*(((x3**2)/(x2-3.221780108334508)))+34.61056223191634*((x3*ln(x3)))-156.78574918200994*(x3)+88.2078474411671*(ln(((x2*x3)+x3)))-1141954.2459589373*(1)+5.321397021927702*((1/x2))'
F1=F11


F=analys(F1,y)
print('////////////////  TEST   \\\\\\\\\\\\\\\\\\')
data=pd.read_excel('data_W.xlsx')
#F1='x2**(-3 - 0.57452212577414014*x3/x2)*(x2**3*(205.7565724686077*x2 + 13471.596529988875) + x2**(0.57452212577414014*x3/x2)*(x2 + 65.473468810063103)*(-19.322430164501*x1 + 19.322430164501*x2 - 388.984140252244*x3) + x2**(3 + 0.57452212577414014*x3/x2)*(x2 + 65.473468810063103)*(-152.37714602747786*(1/(x3*log(log(x2 + x3))))**(x3/x2) + 35.4939507421475*(log(x1*log(x2)) - 0.4675659060773293)**((-25.779664890746375*x1 + 41.38101644139716*x2)/(x1*x2)) + 9.661215082250498*log(log(log(x1))) - 35.91844506635492) - 55.32970290154393*x2**(4 + 0.57452212577414014*x3/x2))/(x2 + 65.473468810063103)'
#F='154089592105988.94*((x1/((x2**4)+(65.473468810063103*(x2**3)))))-25389257920.25107*((x1/((x2**3)+(65.473468810063103*(x2**2)))))+2165.864193314261*((x2*(((1/(x3*ln(ln((x2+x3)))))**(x3/x2))/(x2+65.473468810063103))))+9.04085448087558*((x2*(((ln((x1*ln(x2)))-0.4675659060773293)**(41.38101644139716/x1))/((x2*((ln((x1*ln(x2)))-0.4675659060773293)**(25.779664890746375/x2)))+(65.473468810063103*((ln((x1*ln(x2)))-0.4675659060773293)**(25.779664890746375/x2)))))))-1354.5297223802881*((x2/((x2*(x2**(0.57452212577414014*(x3/x2))))+(65.473468810063103*(x2**(0.57452212577414014*(x3/x2)))))))+20.950845099629454*((x2*(ln(ln(ln(x1)))/(x2+65.473468810063103))))-833.3680001348503*((x2/(x2+65.473468810063103)))-1951595342370114.0*((x3/((x2**4)+(65.473468810063103*(x2**3)))))+218974857229.0796*((x3/((x2**3)+(65.473468810063103*(x2**2)))))+30696169.268169798*((((1/(x3*ln(ln((x2+x3)))))**(x3/x2))/(x2+65.473468810063103)))-139565.00053995333*((((ln((x1*ln(x2)))-0.4675659060773293)**(41.38101644139716/x1))/((x2*((ln((x1*ln(x2)))-0.4675659060773293)**(25.779664890746375/x2)))+(65.473468810063103*((ln((x1*ln(x2)))-0.4675659060773293)**(25.779664890746375/x2))))))-30566018.21804356*((1/((x2*(x2**(0.57452212577414014*(x3/x2))))+(65.473468810063103*(x2**(0.57452212577414014*(x3/x2)))))))'
F='0.0084959721207776*((x1*x3))+1142988.4536007242*((x1/(x1+(58.00041165085283/(x1*x3)))))+6.833524499549682*((x1/x3))+0.6886691064630173*((x2*x3))+18.780084162921256*((x2*ln(((x2*x3)+x3))))-180.40473456387622*(x2)+851.887103495899*((x2/x3))-0.00015724204646610628*(((x3**2)/(x2-3.221780108334508)))+34.61056223191634*((x3*ln(x3)))-156.78574918200994*(x3)+88.2078474411671*(ln(((x2*x3)+x3)))-1141954.2459589373*(1)+5.321397021927702*((1/x2))'
test(F,data)

plt.show()
