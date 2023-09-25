# -*- coding: utf-8 -*-
"""
Created on Mon Jun 19 10:10:06 2023

@author: Nicolas
"""


import numpy as np
import matplotlib.pyplot as plt
import scipy.stats as stats
#import scipy
import pandas as pd
from datetime import datetime

#niveles_camar_2.txt
#fechas_camar_2.txt

# with open("fechas_camar_2.txt") as f:
#     a=[i.strip().split() for i in f.readlines()]
# a2=[]
# for i in range(0,len(a)):
#     for j in range(0,len(a[i])):
#         a2.append(a[i][j])
# data=np.array(a2).astype(float)



def days_passed(year,month,day):
    fixed_date = datetime(year=int(date[0][0]),month=int(date[0][1]),day=int(date[0][2]))
    given_date =datetime(year, month, day)
    days = (given_date - fixed_date).days
    return days
    



data = pd.read_excel("AGUAS DE QUELANA_CAMAR-2.xlsx")

date_list=[]
date=[]
for i in range(0,len(data.iloc[:, 0])):
    date_str, time_str = data.iloc[i, 0].split(' ')
    date_list.append(date_str)
    date_list_split = date_list[i].split('-')
    date.append(date_list_split)

converted_date = []
for i in range(0,len(date)):
    year = int(date[i][0])
    month = int(date[i][1])
    day = int(date[i][2])
    converted_date.append(days_passed(year,month,day))    

#data = converted_date
#%
#plt.plot(np.diff(data))



#%

#import pandas as pd

#niveles_camar_2.txt
#fechas_camar_2.txt

with open("niveles_camar_2.txt") as f:
    a=[i.strip().split() for i in f.readlines()]
a2=[]
for i in range(0,len(a)):
    for j in range(0,len(a[i])):
        a2.append(a[i][j])
data=np.array(a2).astype(float)



#data =data[~np.isnan(data)]

#nombre= "seccion_0"

#data=data[:170]

#data=data[:22]
#data=data[22:120]
#data=data[121:131]
#data=data[132:141]

#data=data[152:]
#converted_date=converted_date[152:]
#data=data[142:]
#data=data[150:]
data=data[165:]

converted_date=converted_date[22:120]

print(np.nanstd(np.diff(converted_date))


#%%
fig=plt.figure()
fig.subplots_adjust(wspace=0.5, hspace=0.5)




#SCATTERPLOT
plt.subplot(2,1,1)
x=np.linspace(0,len(data),num=len(data))
y=data

slope, intercept, r_value, p_value, std_err = stats.linregress(x,data)
line = slope * x + intercept

#plt.plot(x, line, color='red')
#plt.plot(y,".", markersize=1.4)

data = data-line
y=data
plt.plot(x, line, color='red')
plt.plot(y,".", markersize=1.4)


#%
#LAG PLOT
plt.subplot(2,1,2)
#plt.plot(np.diff(converted_date))
plt.plot(np.diff(data))
#%%

data=data[22:120]

converted_date=converted_date[22:120]

print(np.nanstd(np.diff(converted_date))

y = data
#plt.plot(np.diff(y),linewidth=0.5)
#plt.plot(y,".", markersize=1.4)
plt.plot(y- np.nanmean(y)+24,".", markersize=1.4)
plt.plot(np.diff(converted_date)-np.mean(np.diff(converted_date))+24,".",markersize=1.4)


#plt.plot(np.diff(converted_date)-np.mean(np.diff(converted_date))+24,color='red',linewidth=0.4)
plt.xlabel("data sequence")
plt.ylabel("days between data")
plt.xlim(5, 257)
plt.ylim(0, 60)
#plt.ylabel("water level [masl]")
#plt.xlabel("data sequence")
plt.show()


#color='red',linewidth=0.5
