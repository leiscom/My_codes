# -*- coding: utf-8 -*-
"""
Created on Fri Jun  2 14:10:54 2023

@author: Nicolas
"""


import numpy as np
import matplotlib.pyplot as plt
import scipy.stats as stats
import scipy


with open("homework_2.txt") as f:
    a=[i.strip().split() for i in f.readlines()]
a2=[]
for i in range(0,len(a)):
    for j in range(0,len(a[i])):
        a2.append(a[i][j])
data=np.array(a2).astype(float)

nombre= "seccion_0"
#data=data[210:]
#data=data[140:210]
#data=data[:140]
#%
fig=plt.figure()
fig.subplots_adjust(wspace=0.5, hspace=0.5)

#SCATTERPLOT
plt.subplot(2,2,1)
x=np.linspace(0,len(data),num=len(data)+1)
y=data
#y=2*np.cos(((2*(np.pi)*(x))/50)+0.6*np.pi)+ np.random.randn(101)
plt.plot(y,".")

#%
#LAG PLOT
plt.subplot(2,2,2)
lag=1
y=data
#y=2*np.cos(((2*(np.pi)*(x))/50)+0.6*np.pi)+ np.random.randn(101)
plt.scatter(y[:-lag],y[lag:],s=5)
#%

#HISTOGRAM
plt.subplot(2,2,3)
y=data
#y=2*np.cos(((2*(np.pi)*(x))/50)+0.6*np.pi)+ np.random.randn(101)
plt.hist(y)
#%
# NORMAL PROBABILITY
plt.subplot(2,2,4)
N=len(data)+1
x=np.linspace(0,len(data),num=N)
y=data
#y=2*np.cos(((2*(np.pi)*(x))/50)+0.6*np.pi)+ np.random.randn(N)

#stats.probplot(y, dist="norm", plot=plt,marker='o', markersize=3)
plt.scatter(*stats.probplot(data, plot=None)[0], marker='o', s=5)

#%%
#plt.acorr(data, maxlags = 10,marker="o")
#plt.grid(True)
from statsmodels.graphics.tsaplots import plot_pacf
plot_pacf(data, lags=30)
#%%
# ------------------Gaussiana-------------------------------

_, bins, _ = plt.hist(data, density = True)
mu, std = scipy.stats.norm.fit(data) 
best_fit_line = scipy.stats.norm.pdf(bins, mu, std)
    

params = stats.maxwell.fit(bins)
plt.plot(bins,best_fit_line)
plt.plot(bins,stats.maxwell.pdf(bins, *params))
#plt.hist(data)

plt.show()
#%%

params = stats.maxwell.fit(data, floc=0)
_, bins, _ = plt.hist(data, density = True)
plt.plot(bins, stats.maxwell.pdf(bins, *params))
plt.show()


#%%

lag=6
y=data
plt.scatter(y[:-lag],y[lag:],s=5)
#%%

diff = []
days_in_year = 1
for i in range(days_in_year, len(data)):
 value = data[i] - data[i - days_in_year]
 diff.append(value)

plt.plot(diff)
plt.show()
#%%


import seaborn as sns

sns.lineplot(x=np.linspace(0,len(data),num=len(data)),y=data,marker="o")
plt.axhline(np.mean(data), color='red', linestyle='--', label='Mean')
plt.axhline(np.mean(data) + np.std(data), color='green', linestyle='--', label='Mean + Std Dev')
plt.axhline(np.mean(data) - np.std(data), color='blue', linestyle='--', label='Mean - Std Dev')

plt.show()
#%%
plt.savefig('homework_2'+nombre+'.png')


'''
plt.acorr(data, maxlags = 10,marker="o")
'''






