# -*- coding: utf-8 -*-
"""
Created on Mon Sep 11 13:54:29 2017

@author: serge_000


---test part for github
TEST
---

global wTexp
global wEexp
global wSBexp
wTexp=[1,0.024,0.021,0.019,0.018,0.019,0.019,0.015,0.011];
wEexp=[0,0.898,0.905,0.911,0.910,0.909,0.913,0.918,0.925];
wSBexp=[0,0.086,0.082,0.078,0.075,0.071,0.072,0.069,0.064];
"""

# -*- coding: utf-8 -*-
"""
Created on Mon Jun  5 17:09:12 2017

@author: serge_000
"""

import numpy as np
import matplotlib.pyplot as plt
import scipy.optimize as sciopt
from scipy.optimize import fmin
from scipy.integrate import odeint


#Нумерация начинается с 0
#Начальные значения масла (с0[5]), спирта (с0[1]), катализатора (с0[8]) в % соотношении.
c0=[0,0,0,0,0,0,0,0,0];
#Константы скорости (k[0-3]) и коэфициенты массопереноса (B[4-6]). начальные значения
k=[0.1,0.1,0.1,0.1,0.1,0.1,0.1];
#Начальное время реакции
t0=0;
#Время протекания реакции и количество итераций
t = np.linspace(0, 130, 9)   #sec
#Молярные массы веществ        
Mt=872.93;
Me=306.47;
Msb=482.69;

Oil=16.72;
Alc=1.057;
ROH=0.442;
summa_elem=Oil+Alc+ROH;          #18.219
c0[1]=Alc/summa_elem;
c0[5]=Oil/summa_elem;
c0[8]=ROH/summa_elem;
#Функция системым дифференциальных уравнений
def f(c,t,k):
    MYDATA=c0[5]*c0[8];
    cdot=np.zeros(9);
    cdot0=(k[4]*(c[1]-c[0])-k[0]*c[0]*c[5]*c[8]+k[1]*c[2]*c[6]*c[8]);                                      
    cdot1=(-k[4]*(c[1]-c[0]));                                        
    cdot2=(k[5]*(c[3]-c[2])+k[0]*c[0]*c[5]*c[8]-k[1]*c[2]*c[6]*c[8]-k[2]*c[2]*c[5]*c[8]+k[3]*c[4]*c[6]*c[8]);                      
    cdot3=(-k[5]*(c[3]-c[2]));                                         
    cdot4=(k[2]*c[2]*c[5]*c[8]-k[3]*c[4]*c[6]*c[8]);
    cdot5=(-k[0]*c[0]*c[5]*c[8]+k[1]*c[2]*c[6]*c[8]-2*k[2]*c[2]*c[5]*c[8]+2*k[3]*c[4]*c[6]*c[8]);                     
    cdot6=(k[6]*(c[7]-c[6])+k[0]*c[0]*c[5]*c[8]-k[1]*c[2]*c[6]*c[8]+2*k[2]*c[2]*c[5]*c[8]-2*k[3]*c[4]*c[6]*c[8]);                
    cdot7=(-k[6]*(c[7]-c[6]));                                        
    cdot8=(MYDATA/c[5]);
    cdot=[cdot0,cdot1,cdot2,cdot3,cdot4,cdot5,cdot6,cdot7,cdot8]
    return cdot  
#Решение системы ДУ
c1=odeint(f,c0,t,(k,));

         #Пересчет значений из % в массовые единицы (моль/л)
mV=(c1[:,1]+c1[:,0])*Mt+(c1[:,3]+c1[:,2])*Msb+(c1[:,7]+c1[:,6])*Me;
wTa0=c1[:,1]*Mt;
wTp0=c1[:,0]*Mt;
wSBa0=c1[:,3]*Msb;
wSBp0=c1[:,2]*Msb;
wEa0=c1[:,7]*Me;
wEp0=c1[:,6]*Me
       
wTa=wTa0/mV[:];
wTp=wTp0/mV[:];
wSBa=wSBa0/mV[:];
wSBp=wSBp0/mV[:];
wEa=wEa0/mV[:];
wEp=wEp0/mV[:];


wTr=wTa+wTp;
wEr=wEa+wEp;
wSBr=wSBa+wSBp;
#Построение кинетических кривых для исходной системы ду
print('wTr=\n',wTr)
print('wEr=\n',wEr)
print('wSBr=\n',wSBr)
plt.plot(t, wSBr, "b", label="SB")
plt.plot(t, wEr, "g", label="E")
plt.plot(t, wTr, "r", label="T")
plt.legend()
plt.grid(True)
plt.show()

#Оптимизация         
wTexp=np.array([1,0.024,0.021,0.019,0.018,0.019,0.019,0.015,0.011]);
wEexp=np.array([0,0.898,0.905,0.911,0.910,0.909,0.913,0.918,0.925]);
wSBexp=np.array([0,0.086,0.082,0.078,0.075,0.071,0.072,0.069,0.064]);
               
def Objective_function(ode):
   
   ode=odeint(f,c0,t,(k,))
   mV=(ode[:,1]+ode[:,0])*Mt+(ode[:,3]+ode[:,2])*Msb+(ode[:,7]+ode[:,6])*Me;
   wTa0=c1[:,1]*Mt;
   wTp0=c1[:,0]*Mt;
   wSBa0=c1[:,3]*Msb;
   wSBp0=c1[:,2]*Msb;
   wEa0=c1[:,7]*Me;
   wEp0=c1[:,6]*Me
       
   wTa=wTa0/mV[:];
   wTp=wTp0/mV[:];
   wSBa=wSBa0/mV[:];
   wSBp=wSBp0/mV[:];
   wEa=wEa0/mV[:];
   wEp=wEp0/mV[:];

   wTr=wTa+wTp;
   wEr=wEa+wEp;
   wSBr=wSBa+wSBp;

   diffmat0=((wTr-wTexp))**2
   diffmat1=((wEr-wEexp))**2
   diffmat2=((wSBr-wSBexp))**2                
   diffmat=diffmat0,diffmat1,diffmat2
#   print('1=',diffmat)
   Objectivity=np.sum(diffmat)
   print('sum=',Objectivity)
   return Objectivity

optimised=fmin(Objective_function,k,xtol=0.000001)
 
"""
часть из сайлаба с оптимизацией. 
tau=0:2:16
function fun=cost(k)

dc=ode("stiff",c0,t0,tau,f)  
//массовая доля для каждого компонента в двух фазах (а,р)
wTa=dc(2,:)*Mt./mV(1,:);
wTp=dc(1,:)*Mt./mV(1,:);
wSBa=dc(4,:)*Msb./mV(1,:);
wSBp=dc(3,:)*Msb./mV(1,:);
wEa=dc(8,:)*Me./mV(1,:);
wEp=dc(7,:)*Me./mV(1,:);

wTr=wTa+wTp;
wEr=wEa+wEp;
wSBr=wSBa+wSBp;

dc=[wTr;wSBr;wEr];
mess=[wTexp;wSBexp;wEexp];

  diffmat=((dc-mess)).^2
  fun=sum(diffmat);
  disp(fun)
endfunction

function [f1,g,ind]=f3(k, ind)
    f1=cost(k)
    g=numderivative(cost,k)
endfunction

[fopt,k,gopt]=optim(f3,'b',[0;0;0;0;0;0;0],[100;100;100;100;100;100;100],k,'gc')
"""
