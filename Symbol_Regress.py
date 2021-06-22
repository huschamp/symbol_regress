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
warnings.filterwarnings("ignore")

unknown=['x1','x2','x3','x4','x5']
binary_op=['+','-','*','/','^']
unar_op=['U-','ln','exp']
a=random.randint(0,45)
b=random.randint(a+1,100)
start_length=11  #10
max_length=35   #25

num_eqn=8#7 - norm
#num_generation=5
num_generation=10 #20-30 norm
#num_generation=3

data_full=pd.read_excel('data_do_25MPa.xlsx')
data=data_full[data_full.index %5 != 0]

x1=data['T']
x2=data['P']
x3=data['Molar_mass']
x4=data['CO2']
x5=data['N2']

y=data['Z']

##print(list(x1))
##print(list(x2))
##print(list(x3))
##print(list(x4))
##print(list(x5))
##print(list(y))

data_test=data_full[data_full.index.map(lambda x: x not in data.index)]
data_test.to_excel('data_test.xlsx')

x1.index=np.arange(len(x1))
x2.index=np.arange(len(x2))
x3.index=np.arange(len(x3))
x4.index=np.arange(len(x4))
x5.index=np.arange(len(x5))

y.index=np.arange(len(y))

x1=np.array(x1)
x2=np.array(x2)
x3=np.array(x3)
x4=np.array(x4)
x5=np.array(x5)

y=np.array(y)

#X=np.array([(x1*0+1),x1,np.power(x1,2),np.sin(x1),x2,x3,x4]).T

data.to_excel('data_train.xlsx')

mean=sum(abs(y))/len(y)
print('Среднее значение: ', mean)


def is_digit(string):
    if string.isdigit():
       return True
    else:
        try:
            float(string)
            return True
        except ValueError:
            return False


def is_int(str):
    try:
        int(str)
        return True
    except ValueError:
        return False

##def is_int(str):
##    num=float(str)
##    if int(num)==num:
##        return True
##    else:
##        return False

def mnk(X,y):
    b=np.dot(np.dot(np.linalg.inv(np.dot(X,X.T)),X),y)
    print
    y1=np.dot(X.T,b)
    e=y-y1
    S2=np.dot(e.T,e)
    sigma2=S2/(np.shape(X)[0]-np.shape(X)[1])
    S0=np.dot((y-np.mean(y)).T,(y-np.mean(y)))
    R2=1-S2/S0
    return y1,b,R2

def predict(chromosoma,y):
    pred=list([])
    for i in range(len(y)):
        stack=[]
        for gen in chromosoma:
            if gen in binary_op and gen!='^':
                #print(chromosoma)
                #print(stack[-2])
                #print(gen)
                #print(stack[-1])
                elem=eval(str(stack[-2])+gen+str(stack[-1]))
                stack.pop()
                stack.pop()
                stack.append(elem)
            elif gen=='^':
                elem=pow(stack[-2],stack[-1])
                stack.pop()
                stack.pop()
                stack.append(elem)
            elif gen in unar_op:
                if gen=='U-':
                    stack[-1]=-stack[-1]
                if gen=='ln':
                    elem=math.log(stack[-1])
                    stack.pop()
                    stack.append(elem)
                if gen=='exp':
                    elem=math.exp(stack[-1])
                    stack.pop()
                    stack.append(elem)
                if gen=='tanh':
                    elem=math.tanh(stack[-1])
                    stack.pop()
                    stack.append(elem)
            else:
                if is_digit(gen):
                    stack.append(float(gen))
                else:
                    stack.append(eval(gen+'[i]'))
        pred.append(stack[0])
    return np.array(pred)

def rss(chromosoma):
    RSS=0
    for i in range(len(y)):
        stack=[]
        for gen in chromosoma:
            if gen in binary_op and gen!='^':
                try:
                    if stack[-2]=='inf':
                        RSS=float('Inf')
                        return RSS
                    elem=eval(str(stack[-2])+gen+str(stack[-1]))
                    stack.pop()
                    stack.pop()
                    stack.append(elem)
                except:
                    RSS=float('Inf')
                    return RSS
                    #print(gen, 'problem')
                    #break
            elif gen=='^':
                #print(chromosoma)
                if stack[-2]<0 and is_int(str(stack[-1]))==False:
                    RSS=float('Inf')
                    return RSS
                try:
                    elem=pow(stack[-2],stack[-1])
                    if elem=='inf':
                        RSS=float('Inf')
                        return RSS
                    stack.pop()
                    stack.pop()
                    stack.append(elem)
                except:
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
                        return RSS
                    try:
                        elem=math.log(stack[-1])
                        stack.pop()
                        stack.append(elem)
                    except:
                        RSS=float('Inf')
                        return RSS
                        #print('ln problem')
                        #break
                if gen=='exp':
                    try:
                        elem=math.exp(stack[-1])
                        if elem=='inf':
                            RSS=float('Inf')
                            return RSS 
                        stack.pop()
                        stack.append(elem)
                    except:
                        RSS=float('Inf')
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
            return RSS
        #print(type(stack[0]))
        if stack[0]<0:
            RSS=float('Inf')
            return RSS
        try:
            RSS=RSS+abs(stack[0]-y[i])
        except:
            RSS=float('Inf')
            return RSS
    if RSS<0:
        RSS=float('Inf')
    try:
        if math.isnan(RSS):
            RSS=float('Inf')
    except:
        RSS=float('Inf')
    return RSS

def model_of_best_individual(best_individ):
    X=list([])
    #print(np.shape(best_individ))
    for i in range(np.shape(best_individ)[0]): #мб 0 вместо 1
        X.append(predict(best_individ[i],y))
    return np.array(X)


################################################################


def num_operation(chromosoma):
    binar=0
    unar=0
    num=0
    for gen in chromosoma:
        if gen in binary_op:
            binar=binar+1
        elif gen in unar_op:
            unar=unar+1
        else:
            num=num+1
    return [num,binar,unar]



def generate_chromosoma(length):  #уже не будет неправильных хромосом
    I=0
    while I!=1:
        chromosoma=[random.choice(unknown)]
        for i in range(1,length):
            no=num_operation(chromosoma)
            if no[2]//2!=no[2]/2:
                unar=1
            else:
                unar=0
            if no[0]-no[1]-unar<length-i:
                r2=random.random()
                if r2<0.25:
                    chromosoma.append(str(random.uniform(a,b)))
                elif r2>=0.25 and r2<0.5:
                    chromosoma.append(random.choice(unknown))
                elif r2>=0.5 and r2<0.75 and no[0]!=no[1]+1:
                    chromosoma.append(random.choice(binary_op))
                else:
                    chromosoma.append(random.choice(unar_op))
            else:
                chromosoma.append(random.choice(binary_op))
        if rss(chromosoma)!=float('Inf'):
            I=1
    no=num_operation(chromosoma)
    if no[0]!=no[1]+1:
        return generate_chromosoma(length)
    elif no[0]==no[1]+1:
        return chromosoma
        

def cross(chromosoma1,chromosoma2):
    #print('1', len(chromosoma1))
    #print('2', len(chromosoma2))
    r1=random.randint(1,min([len(chromosoma1),len(chromosoma2)]))
    left1=chromosoma1[0:r1]
    right1=chromosoma1[r1:len(chromosoma1)]
    i=0
    no1=num_operation(left1)
    no2=[10,1]
    index=[]
    for i in range(1,len(chromosoma2)):
        no2=num_operation(chromosoma2[0:i])
        if no1[0]-no1[1]==no2[0]-no2[1]:
            index.append(i)
    if len(index)==0:
        return [chromosoma1,chromosoma2]
    index=random.choice(index)
    left2=chromosoma2[0:index]
    right2=chromosoma2[index:len(chromosoma2)]
    chromosoma1=left2+right1
    chromosoma2=left1+right2
    return chromosoma1,chromosoma2


def mutate(chromosoma):
    R=random.random()
    if R>=0.2:
        r1=random.randint(0,len(chromosoma)//4)
        index=random.sample(range(len(chromosoma)),r1)
        for i in index:
            if is_digit(chromosoma[i]):
                r3=random.random()
                if r3<0.33:
                    chromosoma[i]=str(eval(chromosoma[i]+'+'+str(np.random.normal(0,1)*0.01)))
                elif r3>=0.33 and r3<0.67:
                    chromosoma[i]=str(eval(chromosoma[i]+'+'+str(np.random.normal(0,1)*0.1)))
                else:
                    chromosoma[i]=str(eval(chromosoma[i]+'+'+str(np.random.normal(0,1))))
            elif chromosoma[i] in unknown:
                chromosoma[i]=random.choice(unknown)
            elif chromosoma[i] in binary_op:
                chromosoma[i]=random.choice(binary_op)
            else:
                chromosoma[i]=random.choice(unar_op)            
    return chromosoma


##def del_clon(population,RSS):
##    RSS_copy=list(RSS)
##    RSS_new=[]
##    #print(RSS)
##    #print(len(population))
##    population_new=[]
##    for i in range(len(RSS)):
##        RSS_copy.pop(0)
##        if RSS[i] not in RSS_copy:
##            RSS_new.append(RSS[i])
##            #print(len(population))
##            population_new.append(population[i])
##    return [population_new, RSS_new]


def del_clon(population,RSS):
    #RSS_copy=list(RSS)
    RSS_new=list([])
    #print(RSS)
    #print(len(population))
    population_new=list([])
    for i in range(len(RSS)):
        #RSS_copy.append(RSS[i])
        if RSS[i] not in RSS_new:
            RSS_new.append(RSS[i])
            #print(len(population))
            population_new.append(population[i])
    return [population_new, RSS_new]


def sort(population,RSS):
    RSS_sort=sorted(RSS)
    index=[]
    for i in range(len(RSS_sort)):
        if RSS_sort[i]!=float('Inf'):
            index.append(RSS.index(RSS_sort[i]))
        else:
            for i in range(len(RSS_sort)):
                if i not in index:
                    index.append(i)
    #print(RSS)
    #print(index)
    population_sort=list(population)
    #print(index)
    #print(len(population_sort))
    for i in range(len(index)):
        try:
            population_sort[i]=population[index[i]]
        except:
            print(index)
            print(RSS_sort)
            print(len(population))
            print(len(population_sort))
    return population_sort, RSS_sort



def maximum(a):
    max = a[0]
    pos = 0
    for i in range(len(a)):
        if a[i]>max:
            max=a[i]
            pos=i
    return [max,pos]

def minimum(a):
    min = a[0]
    pos = 0
    for i in range(len(a)):
        if a[i]<min:
            min=a[i]
            pos=i
    return [min,pos]

def del2(population,RSS):
    for i in range(2):
        MAX=maximum(RSS)
        population.pop(MAX[1])
        RSS.pop(MAX[1])
    return [population,RSS]


def otbor(population,RSS,m,opts):
    if opts=='rang':
        n=m
        m=n//2
        DELCLON=del_clon(population,RSS)
        population=DELCLON[0]
        RSS=DELCLON[1]
        if len(population)>n:
            n=len(population)
        population=sort(population,RSS)
        RSS=sorted(RSS)
        index1=range(m)
        index_end=range(len(population)-m,len(population),1)
        index=list(index1)
        for i in range(m):
            index.append(index_end[i])
        for k in range(len(index)):
            for j in range(len(index)):
                if k!=j:
                    potomstvo=cross(population[index[k]],population[index[j]])
                    #DEL=del2(population,RSS)
                    #population=DEL[0]
                    #RSS=DEL[1]
                    m1=mutate(potomstvo[0])
                    m2=mutate(potomstvo[1])
                    population.append(m1)
                    RSS.append(rss(m1))
                    population.append(m2)
                    RSS.append(rss(m2))
        population_new=[]
        RSS_new=[]
        for i in range(m):
            population_new.append(generate_chromosoma(start_length))
            RSS_new.append(rss(population_new[i]))
        #unique=unique_generate(population_new,RSS_new,m)
        #population_new=unique[0]
        #RSS_new=unique[1]
        for k in range(len(index)):
            for j in range(m):
                    potomstvo=cross(population[index[k]],population_new[j])
                    DEL=del2(population,RSS)
                    population=DEL[0]
                    RSS=DEL[1]
                    m1=mutate(potomstvo[0])
                    m2=mutate(potomstvo[1])
                    population.append(m1)
                    RSS.append(rss(m1))
                    population.append(m2)
                    RSS.append(rss(m2))
        return [population,RSS]
    elif opts=='inbriding':
        for i in range(m):
            r=random.choice(range(len(population)))
            raznRSS=[]
            for j in range(len(RSS)):
                if j!=r:
                    raznRSS.append(abs(RSS[j]-RSS[r]))
                else:
                    raznRSS.append(float('Inf'))
            j=minimum(raznRSS)[1]
            potomstvo=cross(population[r],population[j])
            DEL=del2(population,RSS)
            population=DEL[0]
            RSS=DEL[1]
            m1=mutate(potomstvo[0])
            m2=mutate(potomstvo[1])
            population.append(m1)
            RSS.append(rss(m1))
            population.append(m2)
            RSS.append(rss(m2))
        return population
    elif opts=='outbriding':
        for i in range(m):
            raznRSS=[]
            r=random.choice(range(len(population)))
            #print(len(RSS))
            #print(r)
            for j in range(len(RSS)):
                if j!=r:
                    raznRSS.append(abs(RSS[j]-RSS[r]))
                else:
                    raznRSS.append(0)
            j=maximum(raznRSS)[1]
            #print('r:',r)
            #print('j:',j)
            #print(len(population))
            potomstvo=cross(population[r],population[j])
            DEL=del2(population,RSS)
            population=DEL[0]
            RSS=DEL[1]
            m1=mutate(potomstvo[0])
            m2=mutate(potomstvo[1])
            population.append(m1)
            RSS.append(rss(m1))
            population.append(m2)
            RSS.append(rss(m2))
        return population



def unique_generate(population,RSS,kolvo):
    while len(population)<kolvo:
        new=generate_chromosoma(start_length)
        population.append(new)
        RSS.append(rss(new))
        DEL=del_clon(population,RSS)
        population=DEL[0]
        RSS=DEL[1]
    return [population,RSS]


def get_eqn(formula):
    stack=deque()
    for gen in formula:
        if is_digit(gen) or gen in unknown:
            stack.append(gen)
        elif gen in unar_op:
            elem=stack[-1]
            stack.pop()
            if gen!='U-':
                stack.append(gen+'('+elem+')')
            else:
                stack.append('-'+elem)
        elif gen in binary_op:
            #print(stack)
            elem2=stack[-1]
            elem1=stack[-2]
            stack.pop()
            stack.pop()
            stack.append('('+elem1+gen+elem2+')')
    eqn=stack[0].replace('--','+')
    eqn=stack[0].replace('+-','-')
    eqn=stack[0].replace('-+','-')
    eqn=stack[0].replace('^','**')
    return eqn



def del_maxL(population,RSS):
    count=0
    i=0
    while i<(len(population)):
        if len(population[i])>max_length:
            count=count+1
            population.pop(i)
            RSS.pop(i)
        else:
            i=i+1
    return [population, RSS]

    
############################################################
def selection(population,RSS,m,opts):
    population_new=list([])
    RSS_new=list([])
    if opts=='rang':
        dooo=(len(population))
        d=del_clon(population, RSS)
        population=d[0]
        RSS=d[1]
        population_s, RSS_s=sort(population,RSS)
        pooosle=(len(population_s))
        if dooo!=pooosle:
            print('не равно')
        index=list(range(m))
        #print(index)
        #index.append(len(population)-2)
        index.append(len(population_s)-1)
        for i in range(10):
            population_s.append(generate_chromosoma(start_length))
            RSS_s.append(population_s[len(population_s)-1])
            index.append(len(population_s)-1)
        #index.append(len(population_s)-4)
        #index.append(len(population_s)-3)
        #index.append(len(population_s)-2)
        #index.append(len(population_s)-1)
        for i in index:
            for j in index:
                if i!=j:
                    
                    chromosoma1, chromosoma2=cross(population_s[i],population_s[j])
                    chromosoma1=mutate(chromosoma1)
                    chromosoma2=mutate(chromosoma2)
##                    population.append(chromosoma1)
##                    RSS.append(rss(chromosoma1))
##                    population.append(chromosoma2)
##                    RSS.append(rss(chromosoma2))
                    population_new.append(chromosoma1)
                    RSS_new.append(rss(chromosoma1))
                    population_new.append(chromosoma2)
                    RSS_new.append(rss(chromosoma2))
                    #print(i)
        population_new.append(population_s[0])
        RSS_new.append(rss(population_s[0]))
        population_new.append(population_s[len(population_s)-1])
        RSS_new.append(rss(population_s[len(population_s)-1]))
    return population_new, RSS_new




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
    eqn=eqn.replace('e-','*10^-')
    eqn=eqn.replace('e+','*10^')
    print('eqn_to: ',eqn)
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




for qqq in range(10):

    best_individual=[]
    best_individual.append(['x1','0','*','1','+'])
    for i in range(num_eqn):
        elite_population=list([])
        elite_RSS=list([])
        population=list([])
        RSS=list([])
        for j in range(30):
            population.append(generate_chromosoma(start_length))
            RSS.append(rss(population[j]))
        count=0
        RSS_prev=0
        for j in range(num_generation):
       #     RSS=list([])
       #     for k in range(len(population)):
       #         RSS.append(rss(population[k]))
            #print(len(population))
            #print(len(RSS))
            DEL=del_clon(population,RSS)
            population=DEL[0]
            RSS=DEL[1]
            
       #     RSS=list([])
       #     for k in range(len(population)):
       #         RSS.append(rss(population[k]))
                
            delMaxL=del_maxL(population,RSS)
            population=delMaxL[0]
            RSS=delMaxL[1]

       #     RSS=list([])
       #     for k in range(len(population)):
       #         RSS.append(rss(population[k]))
            
            unique=unique_generate(population,RSS,30)
            population=unique[0]
            RSS=unique[1]

       #     RSS=list([])
       #     for k in range(len(population)):
       #         RSS.append(rss(population[k]))

            if min(RSS)==RSS_prev:
                count=count+1
            if count>=2:
                for i in range(len(population)):
                    population[i]=mutate(population[i])
                    RSS[i]=rss(population[i])
                count=0

            RSS=list([])
            for k in range(len(population)):
                RSS.append(rss(population[k]))
            print(((min(RSS)/len(y))/mean)*100, '%')
            #print(population[minimum(RSS)[1]])
            #print(len(population))
            elite_RSS.append(min(RSS))
            best_chromosoma=tuple(population[minimum(RSS)[1]])
            elite_population.append(list(best_chromosoma))
            #print(elite_RSS)
            #print(elite_population)
            #for k in range(len(elite_population)):
            #    print('elite: ', (((rss(elite_population[k]))/len(y))/mean)*100, '%')

            RSS_prev=min(RSS)

            #m=len(population)
            m=8
            m=len(population)
            m=len(population)//5
            m=12
    ##        OTBOR=otbor(population,RSS,m,'rang')
    ##        population=OTBOR[0]
    ##        RSS=OTBOR[1]
            population, RSS=selection(population,RSS,m,'rang')
            #LENPOP=30
            #while len(population)>LENPOP:
            #    maxim=maximum(RSS)
            #    population.pop(maxim[1])
            #    RSS.pop(maxim[1])
         #   RSS=list([])
         #   for k in range(len(population)):
         #       RSS.append(rss(population[k]))
                
        best_individual.append(elite_population[minimum(elite_RSS)[1]])
        print(((min(elite_RSS)/len(y))/mean)*100, '%')
        print(elite_population[minimum(elite_RSS)[1]])
    print(best_individual)
    ##
    ##best_RSS=[]
    ##for i in range(1,num_eqn+1):
    ##    rss1=rss(best_individual[i])
    ##    best_RSS.append(rss1)
    ##    print(((rss1/len(y))/mean)*100, '%')

    #best_individual=[['x1', '0', '*', '1', '+'], ['x4', 'x1', 'x4', '+', 'x2', '+', 'x2', '+', 'x2', '/', 'x2', '+', 'x2', '/', '+'], ['x1', 'x1', '+', 'x4', '+', 'x4', '+', 'x2', '/', 'x2', '/'], ['x1', 'x2', 'tanh', 'x4', 'x3', '-', '-', 'x2', 'x2', '/', 'x2', '/', 'x2', '+', 'x2', '/', 'x2', '/', 'x4', '+', '-', '-']]
    #best_individual=[['x1', '0', '*', '1', '+'], ['x4', '8.340002281494696', 'ln', '*', 'x1', '+', 'x1', '+', 'x2', '/', 'x2', '/'], ['x1', 'x1', '+', 'x4', '+', 'x4', '+', 'x2', '/', 'x2', '/'], ['x2', 'ln', 'U-', 'U-', 'U-', 'exp', 'exp', 'x1', '*'], ['x1', '15.065565656604521', 'ln', '*', 'x3', '/', 'x2', '33.333639397896896', '/', '/', 'x4', '+'], ['x4', 'x1', 'x4', '+', 'x1', '+', 'x2', '/', 'x2', '/', '+'], ['x1', 'x4', '+', 'x4', '+', 'x1', '+', 'x2', '/', 'x2', '/'], ['x1', 'x1', '+', 'x4', '+', 'x4', '+', 'x2', '/', 'x2', '/']]
    best_RSS=[]
    for i in range(1,num_eqn+1):
        rss1=rss(best_individual[i])
        best_RSS.append(rss1)
        print(((rss1/len(y))/mean)*100, '%')
    print(best_individual)
    print('DELETE_CLON')

    DEL=del_clon(best_individual,best_RSS)
    best_individual=DEL[0]
    best_RSS=DEL[1]
    print(best_individual)

    for i in range(len(best_individual)):
        print(((best_RSS[i]/len(y))/mean)*100, '%')

    X=np.array(model_of_best_individual(best_individual))
    num_eqn=len(best_RSS)-1


    print(X)
    y1,bb,R2=mnk(X,y)
    print(bb,'Коэф. детерминации: ',R2)
    #print('PREDICT: ', y1)
    string=''
    for i in range(num_eqn+1):
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
    print(eqn)
    print(rss(eqn_to(eqn)))

##from sympy import *
##from sympy.parsing.sympy_parser import parse_expr
##
##x11=list(x1)
##x22=list(x2)
##x33=list(x3)
##x44=list(x4)
##x1,x2,x3,x4=symbols('x1 x2 x3 x4')
###eqn=parse_expr(eqn)
##
##RSS_all1=0
##for i in range(len(y)):
##    eqn_predict=eqn
##    eqn_predict=eqn_predict.subs(x1,x11[i])
##    eqn_predict=eqn_predict.subs(x2,x22[i])
##    eqn_predict=eqn_predict.subs(x3,x33[i])
##    eqn_predict=simplify(eqn_predict.subs(x4,x44[i]))
##    eqn_predict=str(eqn_predict)
##    eqn_predict=eval(eqn_predict)
###    #print('/////////////////')
###    print(eqn_predict)
##    RSS_all1=simplify(RSS_all1+abs(eqn_predict-float(y[i])))
##RSS_all=0
##for i in range(len(y)):
##    RSS_all=RSS_all+abs(y1[i]-y[i])
##print('Итоговая ошибка модели: ',((RSS_all/len(y))/mean)*100, '%')
##print('Итоговая ошибка модели по eqn: ',((RSS_all1/len(y))/mean)*100, '%')
