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

unknown=['x1','x2','x3','x4','x5']
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

data_full=pd.read_excel('data_train_do12.xlsx')
#data_full=pd.read_excel('data_do_12MPa.xlsx')
#data=data_full[(data_full.index % 5) !=0]
#data_full=pd.read_excel('data_norm_alpha007.xlsx')
#data_full=pd.read_excel('data_delete.xlsx')
##data=data_full.sample(frac=0.8)
data=data_full

x1=data['T']
x2=data['P']
x3=data['Molar_mass']
x4=data['CO2']
x5=data['N2']

y=data['Z']

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
                if eqn[q]=='*': #and (eqn[q+1]!='(' and eqn[q+2]!='1' and eqn[q+3]!='0'):
                    j=q+1
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
    eqn=simplify(string)
    #eqn=string
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
        for i in range(length):
            X.append(predict(spl[i],y))
        #print(np.linalg.det(np.dot(X,X.T)))
        X=np.array(X)
        det=np.linalg.det(np.dot(X,X.T))
        if det==0:
           #length=length+1
           spl.pop(length-1)
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
    print('ITOG: ')
    print('///////////////////////////////////////////////////////')
    print(F)
    print('///////////////////////////////////////////////////////')
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
    ax.set_xlabel('Коэффициент сжимаемости Z')
    ax.set_ylabel('Абсолютная ошибка')
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
    fig7, ax5=plt.subplots()
    ax5.scatter(x4,ae)
    ax5.set_title('Зависимость остатков от молярной доли CO2')
    ax5.set_xlabel('CO2')
    ax5.set_ylabel('Абсолютная ошибка')
    fig8, ax6=plt.subplots()
    ax6.scatter(x5,ae)
    ax6.set_title('Зависимость остатков от молярной доли N2')
    ax6.set_xlabel('N2')
    ax6.set_ylabel('Абсолютная ошибка')
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
    ax_ard.set_ylabel('Накопительная частота для коэффициента сжимаемости, (%)')
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
    return r2


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
    x1=data['T']
    x2=data['P']
    x3=data['Molar_mass']
    x4=data['CO2']
    x5=data['N2']
    print(len(x1))
    y=data['Z']

    mean=sum(abs(y))/len(y)
    print(mean)
    y_mean_y=list([])
    for i in range(len(y)):
        y_mean_y.append(y[i]-mean)
    S2_0=sum( i*i for i in y_mean_y)
    n=len(y)
    m=12
    #pred_y=predict(eqn_to(F),y)
    

    pred_y=list([])
    x11=list(x1)
    x22=list(x2)
    x33=list(x3)
    x44=list(x4)
    x55=list(x5)
    x1,x2,x3,x4,x5=symbols('x1 x2 x3 x4 x5')
    eqn=parse_expr(F)
    for i in range(len(y)):
        eqn_predict=eqn
        eqn_predict=eqn_predict.subs(x1,x11[i])
        eqn_predict=eqn_predict.subs(x2,x22[i])
        eqn_predict=eqn_predict.subs(x3,x33[i])
        eqn_predict=eqn_predict.subs(x4,x44[i])
        eqn_predict=simplify(eqn_predict.subs(x5,x55[i]))
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
    ax.set_xlabel('Коэффициент сжимаемости Z')
    ax.set_ylabel('Абсолютная ошибка')


###
best_individual=[['x1', '0', '*', '1', '+'], ['x2', '59.4391132107437', '26.957909071646064', '/', '27.859513595746296', '-', 'x1', '/', '^'], ['x4', 'ln', 'x4', '*', 'x2', '+', 'ln', 'x4', '*', 'x2', 'x1', '/', 'x4', 'U-', '-', 'x4', 'U-', '-', '^'], ['x5', '64.71990513901052', 'exp', 'x3', 'exp', '+', 'ln', '*', 'x5', '56.04197426467418', '55.58752202084953', '19.55019705929017', 'x3', '-', '/', '-', '*', '^'], ['x4', 'x4', '33.56339531095271', '*', '*', 'x2', '40.72370811429743', '/', 'U-', 'exp', '+'], ['x1', 'ln', 'exp', 'exp', 'ln', '31.423828759940843', 'x1', 'x4', '+', 'x2', '+', 'x2', '+', 'x2', '+', '+', '/'], ['x2', '60.93687492617446', '/', 'exp', 'x4', 'exp', 'U-', '19.853007199966473', 'U-', '22.578956970666276', '+', '50.69004181274559', '64.26919584184166', '-', '/', '+', '19.853328906080087', 'U-', '22.368689366848944', '+', '51.48931447299082', '64.05379647790181', '-', '/', '+', '^'], ['x2', 'exp', 'x2', '-', 'x2', '*', 'x2', 'ln', 'U-', '14.485631907899126', 'x2', '*', '/', '^']]
best_individual=[['x1', '0', '*', '1', '+'], ['x1', '43.25759580494125', '-', 'x2', '24.698124465926423', 'x4', '49.14789082398188', '^', '^', 'x1', '*', '+', '/'], ['x2', 'x2', '*', '13.268409331704477', 'U-', 'x1', '/', '^'], ['x3', 'x2', '60.00498083584014', 'U-', 'x1', '26.385737671560577', '*', '26.38025230625069', '-', '/', '*', '^'], ['x2', '48.45917664781548', 'U-', '/', 'exp'], ['x4', 'x4', 'x1', '/', '-', 'x2', 'x1', '/', '^'], ['x2', 'x3', 'ln', '24.411300695801316', '-', 'x1', '58.23288516264744', '-', '/', '^', 'x4', '35.55002553974635', 'x3', '/', 'ln', '^', '+']]
best_individual=[['x1', '0', '*', '1', '+'], ['x2', '54.12353893345937', 'ln', '*', 'ln', 'x3', '/', 'U-', 'exp', 'x4', '-'], ['x5', 'x2', 'x3', 'x5', '^', '42.29930595614707', 'x2', '/', '^', '^', 'x1', '/', '^'], ['x4', 'ln', 'x1', '31.61653162417812', '/', '/', '68.13660288682678', '/', 'exp', 'x2', '49.362480103880216', '/', 'exp', '/'], ['x5', 'exp', 'exp', 'x3', 'ln', 'x5', '+', 'x3', 'U-', 'x1', '23.738223976503896', '+', '/', '*', '^'], ['x3', 'exp', 'x4', '^', 'x2', '+', '51.102955748190126', 'U-', '/', 'x4', '*', 'exp', '30.813042312218343', '*', 'x5', '*', 'exp', 'x2', '+', '51.11626193244735', 'U-', '/', 'exp']]
best_individual=[['x1', '0', '*', '1', '+'], ['x2', 'x4', '^', '67.14589585564076', 'x1', '38.24831854549455', '-', '/', 'U-', 'exp', '*'], ['x1', 'x1', '+', 'x4', '*', 'x4', '+', 'x1', '*', 'x1', '+', 'x4', '*', 'x5', '*', 'exp', 'x1', '*', 'x1', '75.06354809996076', '+', '/'], ['x1', 'x1', '63.1788580091158', '+', '/'], ['x1', 'U-', 'U-', 'x1', '60.35400831529689', '+', '/'], ['x4', 'exp', 'U-', 'exp', '78.15418976588987', '48.42334512887181', '67.44718331384267', 'x1', '*', '/', '*', '^'], ['x2', 'x1', 'ln', '73.79258526182392', '/', '^', '80.65566673823368', 'x5', '+', 'exp', '-', 'exp', 'x2', 'ln', '57.00264811071095', 'U-', '/', '+', 'x4', 'x4', '^', '+'], ['x3', 'ln', 'ln', '74.55423786715318', 'x3', '42.888811075887496', '69.1294494290081', '-', '/', '+', '81.31111932502762', '83.10493211500571', 'x1', '-', '-', '/', '-']]

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


##################################################################

F1='2.8738510259628773*(((x1*0)+1))+8.983642816860002*((x2**(((59.4391132107437/26.957909071646064)-27.859513595746296)/x1)))+0.08291696054611009*(((ln(((ln(x4)*x4)+x2))*x4)**(((x2/x1)-(-x4))-(-x4))))-0.23194124506202524*(((x5*ln((exp(64.71990513901052)+exp(x3))))**(x5*(56.04197426467418-(55.58752202084953/(19.55019705929017-x3))))))+0.1834166153281256*(((x4*(x4*33.56339531095271))+exp((-(x2/40.72370811429743)))))-5.278060356766509*((ln(exp(exp(ln(x1))))/(31.423828759940843+((((x1+x4)+x2)+x2)+x2))))+0.9191206744844256*((exp((x2/60.93687492617446))**(((-exp(x4))+(((-19.853007199966473)+22.578956970666276)/(50.69004181274559-64.26919584184166)))+(((-19.853328906080087)+22.368689366848944)/(51.48931447299082-64.05379647790181)))))-7.036306887871249*((((exp(x2)-x2)*x2)**((-ln(x2))/(14.485631907899126*x2))))'
#
F1='0.4829128986469927*(((x1*0)+1))-0.36315364368077374*(((x1-43.25759580494125)/(x2+((24.698124465926423**(x4**49.14789082398188))*x1))))+0.5600084669301489*(((x2*x2)**((-13.268409331704477)/x1)))+4.4591811784934245*((x3**(x2*((-60.00498083584014)/((x1*26.385737671560577)-26.38025230625069)))))-4.054444184045858*(exp((x2/(-48.45917664781548))))+0.09164303445690077*(((x4-(x4/x1))**(x2/x1)))-0.21462319522809062*(((x2**((ln(x3)-24.411300695801316)/(x1-58.23288516264744)))+(x4**ln((35.55002553974635/x3)))))'
#F1='0.4829128986469927*(((x1*0)+1))-0.36315364368077374*(((x1-43.25759580494125)/(x2+((24.698124465926423**(x4**49.14789082398188))*x1))))+2.0578382469502515*(x1**2)*ln(ln(x2))+0.5600084669301489*(((x2*x2)**((-13.268409331704477)/x1)))+4.4591811784934245*((x3**(x2*((-60.00498083584014)/((x1*26.385737671560577)-26.38025230625069)))))-4.054444184045858*(exp((x2/(-48.45917664781548))))+0.09164303445690077*(((x4-(x4/x1))**(x2/x1)))-0.21462319522809062*(((x2**((ln(x3)-24.411300695801316)/(x1-58.23288516264744)))+(x4**ln((35.55002553974635/x3)))))'
#F1='0.4829128986469927*(((x1*0)+1))-0.36315364368077374*(((x1-43.25759580494125)/(x2+((24.698124465926423**(x4**49.14789082398188))*x1))))+2.0578382469502515*((x1**2)*ln(ln(x2)))+0.5600084669301489*(((x2*x2)**((-13.268409331704477)/x1)))+4.4591811784934245*((x3**(x2*((-60.00498083584014)/((x1*26.385737671560577)-26.38025230625069)))))-4.054444184045858*(exp((x2/(-48.45917664781548))))+0.09164303445690077*(((x4-(x4/x1))**(x2/x1)))-0.21462319522809062*(((x2**((ln(x3)-24.411300695801316)/(x1-58.23288516264744)))+(x4**ln((35.55002553974635/x3)))))+0.35297303898788196*((exp((-(ln((x2*ln(54.12353893345937)))/x3)))-x4))+2.4727533673591253*((x5**((x2**((x3**x5)**(42.29930595614707/x2)))/x1)))+0.18615051137189825*((exp(((ln(x4)/(x1/31.61653162417812))/68.13660288682678))/exp((x2/49.362480103880216))))+2.250404832345879*((exp(exp(x5))**((ln(x3)+x5)*((-x3)/(x1+23.738223976503896)))))-1.7945387919386058*(exp(((exp(((exp(((((exp(x3)**x4)+x2)/(-51.102955748190126))*x4))*30.813042312218343)*x5))+x2)/(-51.11626193244735))))'
#

F1='-2.114202360194956*(((x1*0)+1))+0.35297303898788196*((exp((-(ln((x2*ln(54.12353893345937)))/x3)))-x4))+2.4727533673591253*((x5**((x2**((x3**x5)**(42.29930595614707/x2)))/x1)))+0.18615051137189825*((exp(((ln(x4)/(x1/31.61653162417812))/68.13660288682678))/exp((x2/49.362480103880216))))+2.250404832345879*((exp(exp(x5))**((ln(x3)+x5)*((-x3)/(x1+23.738223976503896)))))-1.7945387919386058*(exp(((exp(((exp(((((exp(x3)**x4)+x2)/(-51.102955748190126))*x4))*30.813042312218343)*x5))+x2)/(-51.11626193244735))))'
#F1='-2.836294083863734*(((x1*0)+1))+0.6181205625401072*((ln(((((ln(((((exp(x1)*x5)*(ln(58.20705647268794)+(-x4)))-(ln(58.20705647268794)**41.97059512247089))*x3))+x3)*(ln(57.673119590206504)/41.97059512247089))**x3)*x3))/67.54729800710263))-1.5996111602967478*(exp((-(exp(ln(ln((exp((-exp((x2**x2))))*((x3/(x1+x3))*x2)))))/34.12221215209524))))+1.0183038940325275*((x4+(x5*(x1/(x5*(43.10424708126337+x1))))))+1.0069889349857781*(((ln(exp(x1))/(x1-((x4-(x2-x5))*x5)))-x5))-4.2130154731166165*((exp(((-(-(-(exp(ln(exp((-(ln(x1)/(-((x5*x5)-52.28024321764542)))))))/x1))))*47.42302972276099))+x5))+8.180833887914343*((((x3**x4)**x4)-((((((x5+(x3+x3))/(ln((x3-x5))+x1))*x3)+(x3+x3))+ln((x3-x5)))/(ln(ln(x3))+x1))))-0.9216412212636214*((exp(x5)**(((-x2)/x1)-(37.85919980867566**(x5-(x1/49.23391152378004))))))'
#F1='-6.543620841184885*(((x1*0)+1))+11.020825003975805*((exp((-exp((-x3))))-(38.663699656026786/(x1-44.0405721365802))))+0.47421380700612525*((x4**(x2/(43.04976904745683*(68.68058237498849*(52.65212485445105+71.11449120143945))))))+2.2746097094671516*((((((((-(-x1))-x3)**exp((ln((x3/x2))/ln(55.89845699710906))))+x4)*x5)**x4)+(ln((x3/x2))/55.91693578179563)))+3.2272371787953027*((ln((x2/x1))**(ln(x5)/(x3*(51.362338968665455/x3)))))-6.943204167316427*((((ln(x5)+x1)+(x1/38.342885583606346))/(x1+47.24127461620375)))-1.4049471631322197*((x5**((x4+(x2/(x1-(61.53270764842784+x4))))/(31.718563697975515*44.945581097107066))))'
# хорошая формула

#F1='-10.237922639821512*(((x1*0)+1))+1.7043370290438766*(((exp((ln(x5)/41.386005361430165))-x5)*(x1**(44.13154295785829/x2))))-1.5596116068637471*((exp(x5)/(36.13431209211308/31.694931007170236)))+30.57664930641198*((ln(x2)**((32.53911029594708-49.791392401323776)/x1)))-1.8065367694283705*(((x2+(x2-((x2/(45.69539893339916/x4))+(40.38615576678313*(47.8664019357436-x1)))))**(x2/((x1/(47.2999817174788/x5))-(40.38615576678313*(47.8664019357436*x1))))))-34.41211571498984*((x1/(49.14563650895685+(ln(exp(x1))+(x3/x1)))))+16.78652689849787*((x2**(ln(38.778648893790184)/((-x1)+(x5/46.77593062264902)))))+0.07251557267116682*((exp(x5)**(x2/(x4-x1))))'
#F1='-31.886057325537898*(((x1*0)+1))-5.553393873060122*((x1/(x1+(43.15586433033245*(x2/(x1*32.12606386785891))))))+0.35443988764328305*(((exp(x5)*(37.3758578306126/43.517570809497464))+(39.471885369975666/(x2-((x5+32.9896541358815)*48.66961570158783)))))+64.12178011344389*((x2**ln(((x1-ln(46.20666453958798))/x1))))-33.708010521230825*(((x1**((35.22902250826692+x1)/(39.141444942624105*x2)))+(49.310478380546705/(x4-(33.97779255434196+(x1-x4))))))+10.576449232625375*(((x3/x2)**((ln(x2)/ln(x2))/42.599240709341245)))+0.07936696590524939*((((x4+(35.710026476011535/40.8341755243675))-(31.80844837540896/(44.68435933348356-x2)))-((ln(ln(31.93980413659144))*39.82238679434423)/(44.68435933348356-x2))))'
#^R2=0.996

#F1='-6.811590476939822*(((x1*0)+1))-5.472895309762319*(((x2/(x5-(-x1)))**(11.584694517357423/(-(x1-x3)))))+0.7368158142502752*(((x1/18.349262047744332)/(((((-(x2*x4))*exp(19.0543874977051))/(x2*((((15.794942662373673+x2)-x4)-x4)-x3)))/(18.88779044875449*15.794942662373673))+18.88779044875449)))-0.05371960312346297*((ln(ln((exp((-x5))*17.12576854658968)))**(x2/((-11.60733866614392)/x5))))+0.38188683104138654*(((((x1-(((-x1)+15.149299584035875)/x2))-(12.288815510184055/(13.199471188713318+14.970369811385137)))-(12.422650882609414-(17.22889277618782-((x1*17.659827082525837)**x5))))/(12.422650882609414*(13.199471188713318+14.970369811385137))))+8.058073820261441*(((((x1*x3)/(x2+x4))+15.405273360945493)/18.49593483882183))+0.3117285691096283*((ln((exp(((x5*x2)+(x4-(14.08011783977349/x1))))-11.978698732558101))**(x5-(14.067498438909281/x1))))+4.939086385435292*((exp((ln((exp(x3)-(12.08185351432107-exp(x4))))/x1))**((ln(x3)/x3)*(-x3))))'
#F1='1.4581188277515824*(((x1*0)+1))+0.5546535697083361*(((x1/(exp(exp((x1-(x2-14.234629608276421))))-((-x5)-18.077744520060683)))/(exp(ln((15.586649081051196/x3)))-((-x5)-17.832418862685838))))+3.7660788969310204*((ln(x2)**(x3/(x3-x1))))-0.11678506873987338*((ln((-(((-(-((ln(exp(exp((x3/(-ln(11.401151748067132))))))*x2)*ln(11.385557626980688))))*x4)-ln(11.153646541624962))))**(x2**(x4**x3))))-4.362584293278378*((exp((-x3))+((12.615773550279267/16.967841951435396)+(ln((19.65228191786348/18.60975846943178))+(x1/x2)))))-1.149985889576258*((ln(((x2-x4)-(x4*(x5**(x4+12.019921855804267)))))**(14.718877472971995/(x3-(x1-11.451264462865238)))))+0.30438071466895433*(((ln(exp(x5))*((-x5)*(-x1)))+(15.595994257136852**((-15.25498927458025)/x1))))'
##F1='-5.285193117418444*(((x1*0)+1))+0.28449891041535214*((x5**(ln(x2)/(((x1-(-19.88472918911195))-ln((19.878943018338788/x1)))+11.208928859456671))))+0.7509876059102116*(((((((((x3+x4)**x4)-17.388136091813553)+17.251900687135418)+x5)+(17.244688473784073*(19.309857287446782/x2)))+((17.244688473784073*(19.231741025253207/x2))/(13.303716243927019-15.55398339146635)))-x4))-0.5868858124872336*((((ln(x2)/(((10.315748623549231*x5)*(x1*ln(x5)))+18.944964342697652))+17.644247377625565)/(ln(16.68262346928645)+17.000630516423477)))+7.027145254672143*(((ln((x2-18.606112542529154))**((x4-18.305574428949136)/x1))/((15.202357926522117/15.296296352781392)+exp((ln(19.75562806408235)-x2)))))-0.6359212164692649*(((exp((x2/((-(-(13.83026437652837*(x1+(exp((-18.216902392799316))**13.372913351442651)))))*(-14.830354550816475))))-(-x5))*(17.06019159356653/(exp(x4)*18.166776510239483))))'
##F1='4.748903037840183*(((x1*0)+1))+1.7003490291956695*(((-(-((x1/18.728874360220154)/18.845535704484224)))**(x3/18.21165821537276)))+0.17085925810629277*((x5**(((((((((18.200804336092716/x2)+x1)/x3)/11.99760850264991)+12.000945385720973)+x3)*12.000945385720973)/exp(12.000945385720973))*12.1230938632386)))-1.1049850603564466*((ln(exp(x4))+((15.071486604074817/(16.679660963985594+(18.94897705511373-x3)))+(16.661878845092758*(19.038806968822623/x2)))))-5.198900158006589*((ln((x2/(ln((x2**14.006669045617176))*(ln(19.40329415173085)/19.651268363099774))))**(x3/(12.47056422391823-x1))))'

Fo1='-2.8746754350260186*(((x1*0)+1))+5.108681099788843*((x1/((32.38459837987492+(34.23210891897487+x1))-(22.22171034083161/x2))))-0.6066399581732421*(((ln(ln((x1-x3)))+x3)/22.815431853789356))'
Fo1='-1.1221419515510496*(((x1*0)+1))+3.3333424442763446*((exp((exp(exp(x4))/x1))**(-x3)))-3.55108111051376*(exp((-(ln((x2+(26.274044566876206*33.71401097674961)))/33.723526825936155))))+2.1154983251181547*(((-x1)/((x4-(x1+33.965041966248066))-(x5+33.95356008935353))))'
Fo1='9.654357622134233*(((x1*0)+1))-70.21084885476316*((x3**(20.8258498909017/((x5-x1)+x4))))-26.41705250131978*(exp((ln((ln((x2/x3))+30.304537373404056))/(-x3))))+76.46471068156514*(((exp((exp(x5)/21.367816795042362))/21.58174607067223)**(22.059275106944266/(x1+20.180928704761676))))+9.784840496324517*((x1/((ln(37.59530075710155)-(-36.001227976635136))/(32.716241031060456/x1))))'
Fo1='125.61130504223001*(((x1*0)+1))-12.409886785891558*((ln((exp(ln((exp(x3)-24.42338019850099)))-24.42338019850099))/20.64694972006998))+16.269288308807518*(((x3/22.314712806166934)**(20.55100056216242/28.167937095598198)))-205.9448954145014*(exp(((x4**(x5/(27.94655200292724+35.781539264726185)))/(ln((x5+ln(ln(x1))))-ln(x1)))))+69.97409429870851*(((exp(ln(exp(exp((-ln(exp((x3/32.3996996959575))))))))/32.81732861978454)**(22.974767861629843/(35.06570132437192+x1))))-21.267178146154308*((x3**(18.9429209883298/(24.07865115456301-x1))))+2.197267714240466*((((x4*x5)**(x4/26.49581386600017))*(27.108173401786612/32.7740265929917)))'
#F1='311.5800049615898*(((x1*0)+1))+54.51000072404655*(exp((ln(ln((x3*x3)))/((-x1)/34.30470256837824))))-38.02776730041026*(exp((x3/(((-x1)/ln(26.580727143370808))+ln(x2)))))+116.01567356248489*((x3**(21.06781474992139/((((((x1/34.268326581378524)-x3)-x4)/x2)-x5)-x1))))-20.552986183556357*((exp(ln((ln(ln(x2))+x3)))**(ln(26.192024787254827)-ln(28.09361323584863))))+9.341035602291214*((((x3**((x1/(x3+29.90951570446966))/x1))**x4)-((27.368821718562508+32.392173034044525)/x1)))-105.39668899064262*((ln((exp(ln(x2))*ln((34.28398535826342**x3))))**(30.723165875757804/(-x1))))-330.6417094618715*((x2**((27.10737022130983/(27.21919094018781-34.44769294926292))/(x2+33.434873732897856))))'
# ^ R2=0.992
Fo1='21.76326735143103*(((x1*0)+1))-0.14295810592208424*((x3/(ln(9.444483384170635)**ln(42.32020858202117))))-0.9453320615439651*(((exp(ln((exp(exp(ln((exp((x1/77.99058553120888))/59.9873082217172))))/59.9873082217172)))/34.08351292540209)**(ln((12.384744111524533/61.52695439693147))/(x3-73.4791652212227))))+4.834878835483044*((x3**((-(21.903344637381903-x5))/((x1-36.84257252102647)+35.571480036826514))))-23.918712001370196*(((exp((-(-(-x4))))+x3)**(x3/(x2-((x1-x3)-x3)))))'
Fo1='28.753361647433167*(((x1*0)+1))+11.610598917623463*(exp((-exp((-exp((exp((-exp((x4+x4))))+x4)))))))-151.13922372271463*((((ln((((x4+x1)-59.07404347036938)+(-x3)))+x1)-58.96159515394046)/x1))+184.8573626989987*(((exp((-x4))+41.2800199591034)**(14.225204797330747/((66.16939261996825-41.40316736368232)-x1))))+11.02827497125136*((ln(ln(ln(((x1*(ln(25.108391735436122)*x1))/x3))))+(x2**(x4-77.32750962191413))))-91.39566784177964*((ln(ln(ln(x1)))**(ln(62.87016819875507)/11.838711701952525)))'

Fn1='0.0099070718333366545*x2 - 0.064862749991825841*x3 - 9.0200417513619407*(0.087117335417587058*x1 + 1)**(-0.020914216842078989) + 21.468157230561715 - 26.858956820827107*(10.83684915050508*x3 - 10.589883641159599)**(-11.684194184171247/(x2 + 11.591870486726181)) + 9.135671516945788*x1**(-log(x2)/(x2 - 42.970522664487596))'
#^singular


#F0='x2**(-3 - 0.57452212577414014*x3/x2)*(x2**3*(205.7565724686077*x2 + 13471.596529988875) + x2**(0.57452212577414014*x3/x2)*(x2 + 65.473468810063103)*(-19.322430164501*x1 + 19.322430164501*x2 - 388.984140252244*x3) + x2**(3 + 0.57452212577414014*x3/x2)*(x2 + 65.473468810063103)*(-152.37714602747786*(1/(x3*log(log(x2 + x3))))**(x3/x2) + 35.4939507421475*(log(x1*log(x2)) - 0.4675659060773293)**((-25.779664890746375*x1 + 41.38101644139716*x2)/(x1*x2)) + 9.661215082250498*log(log(log(x1))) - 35.91844506635492) - 55.32970290154393*x2**(4 + 0.57452212577414014*x3/x2))/(x2 + 65.473468810063103)'
#^ супер формула с нормальным распределением

#F1='-1.2333313304792055*(x1**1.8673873825529029)**(-0.011198325008742172) + 0.4421628362690665 + 1.717736547498913*x1**(-6.587599830386579/x2)'
#F1='1231*x1+4251*x2+1.7003490291956695*(((-(-((x1/18.728874360220154)/18.845535704484224)))**(x3/18.21165821537276)))+0.17085925810629277*((x5**(((((((((18.200804336092716/x2)+x1)/x3)/11.99760850264991)+12.000945385720973)+x3)*12.000945385720973)/exp(12.000945385720973))*12.1230938632386)))-1.1049850603564466*((ln(exp(x4))+((15.071486604074817/(16.679660963985594+(18.94897705511373-x3)))+(16.661878845092758*(19.038806968822623/x2)))))-5.198900158006589*((ln((x2/(ln((x2**14.006669045617176))*(ln(19.40329415173085)/19.651268363099774))))**(x3/(12.47056422391823-x1))))+1.9683189827588907*x2/(x2 + 67.26211779128894)+3.3333424442763446*((exp((exp(exp(x4))/x1))**(-x3)))-3.55108111051376*(exp((-(ln((x2+(26.274044566876206*33.71401097674961)))/33.723526825936155))))+2.1154983251181547*(((-x1)/((x4-(x1+33.965041966248066))-(x5+33.95356008935353)))) +x3 +9.661215082250498*log(x2+x1) + 0.6685241221355813*(x1 + x3)**(-0.01958839843541684) + 1.717736547498913*x1**(-6.587599830386579/x2) - 1.3372249644701104'
F1='-5.553393873060122*((x1/(x1+(43.15586433033245*(x2/(x1*32.12606386785891))))))+0.35443988764328305*(((exp(x5)*(37.3758578306126/43.517570809497464))+(39.471885369975666/(x2-((x5+32.9896541358815)*48.66961570158783)))))+64.12178011344389*((x2**ln(((x1-ln(46.20666453958798))/x1))))-33.708010521230825*(((x1**((35.22902250826692+x1)/(39.141444942624105*x2)))+(49.310478380546705/(x4-(33.97779255434196+(x1-x4))))))+10.576449232625375*(((x3/x2)**((ln(x2)/ln(x2))/42.599240709341245)))+0.07936696590524939*((((x4+(35.710026476011535/40.8341755243675))-(31.80844837540896/(44.68435933348356-x2)))-((ln(ln(31.93980413659144))*39.82238679434423)/(44.68435933348356-x2))))-52.84371649127911*(((x1*0)+1))+3.784542340080558*(((x2**x4)*exp((-(67.14589585564076/(x1-38.24831854549455))))))+0.04153364634297585*(((exp((((((((x1+x1)*x4)+x4)*x1)+x1)*x4)*x5))*x1)/(x1+75.06354809996076)))-1351.5679878470228*((x1/(x1+63.1788580091158)))+1436.667068304435*(((-(-x1))/(x1+60.35400831529689)))-40.17344614036813*((exp((-exp(x4)))**(78.15418976588987*(48.42334512887181/(67.44718331384267*x1)))))+6.228954021992455*(((exp(((x2**(ln(x1)/73.79258526182392))-exp((80.65566673823368+x5))))+(ln(x2)/(-57.00264811071095)))+(x4**x4)))-0.9804882493785527*((ln(ln(x3))-((74.55423786715318+(x3/(42.888811075887496-69.1294494290081)))/(81.31111932502762-(83.10493211500571-x1)))))'



F1=F1.replace('+-','-')
F1=F1.replace('-+','-')
F1=F1.replace('--','+')





#analys(F1,y)


print('////////////////  TEST   \\\\\\\\\\\\\\\\\\')
data=pd.read_excel('data_do_12MPa.xlsx')
#F1='x2**(-3 - 0.57452212577414014*x3/x2)*(x2**3*(205.7565724686077*x2 + 13471.596529988875) + x2**(0.57452212577414014*x3/x2)*(x2 + 65.473468810063103)*(-19.322430164501*x1 + 19.322430164501*x2 - 388.984140252244*x3) + x2**(3 + 0.57452212577414014*x3/x2)*(x2 + 65.473468810063103)*(-152.37714602747786*(1/(x3*log(log(x2 + x3))))**(x3/x2) + 35.4939507421475*(log(x1*log(x2)) - 0.4675659060773293)**((-25.779664890746375*x1 + 41.38101644139716*x2)/(x1*x2)) + 9.661215082250498*log(log(log(x1))) - 35.91844506635492) - 55.32970290154393*x2**(4 + 0.57452212577414014*x3/x2))/(x2 + 65.473468810063103)'
F='1.654714155006305*x1/(x1 + 1.343328722368273*x2/x1) - 413.4195153599808*x1/(x1 + 63.1788580091158) + 410.1236695794414*x1/(x1 + 60.35400831529689) + 24.573016064170886*x1**(0.90004399581844128/x2)*x1**(0.025548366992221683*x1/x2) + 0.38250482960151544*x2**x4*exp(-67.14589585564076/(x1 - 38.24831854549455)) + 49.4492042952632*x2**log(1 - 3.833124041787706/x1) - 27.303588872281978*x3/(x1 - 1.793812789978089) - 6.710466390375586*x4 + 448.7183811145337*(x3/x2)**0.023474596808499407 + 0.4683682600741399*exp(x5) - 28.523058461976305*exp(-exp(x4))**(56.11038324152387/x1) + 9.599465103332037*log(x2) - 22.52960191268849*log(log(x3)) - 495.77504042071814 - 20.259521626142316/(-x2 + 48.66961570158783*x5 + 1605.59378892165) - 2220.6459698414305/(x2 - 44.68435933348356)'

test(F,data)

plt.show()
