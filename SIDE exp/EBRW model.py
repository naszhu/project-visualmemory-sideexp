#!/usr/bin/env python
# coding: utf-8

# # Set up 

# In[1]:


import numpy as np
import pandas as pd
import time
import matplotlib.pyplot as plt
import sys
# from rpy2.robjects.packages import importr, data
# from rpy2.robjects.lib.dplyr import DataFrame


# In[2]:


path="./" + "sidedata/"
filename = path + "Alldata.csv"

df=pd.read_csv(filename, index_col=None)
df['Probtype'] = np.where(df['Stimkind']==1, "CM",
                   np.where(df['Stimkind']==0, "AN",
                   np.where(df['Stimkind']==2, 'VM',"wrong")))

df['Oldnew'] = np.where(df['Old']==1, "old",
                   np.where(df['Old']==2, "new","wrong"))
df['Error'] = 1-df['Correctness']
# df[["Correctness","Error"]]
# df.groupby("Error")
df.columns


# # Calc A (activation)

# In[87]:


def calcA(df, name):
    
#     reveal_param(param_dic)
#     for key,val in param_dic.items(): exec(key + '=val')
    
    dnow = df[df['FileCondi']==name]
    
    m = [(j**(-beta) + alpha) for j in np.arange(1,9)]

    a = np.zeros((dnow.shape[0], 8)) 
    for j in np.arange(1,9):

        indexj = j-1 
        a[dnow['Lag'] == j,indexj] = m[indexj] # if tested
        a[dnow['Lag'] != j,indexj] = m[indexj] *s #if not tested
        
    for i in range(a.shape[0]): a[i,dnow['Setsize'].iloc[i]:] = 0 
    
    A = a.sum(axis = 1)
    
#     print(dnow['Lag']==dnow['Serpos'])
    
    return(A)

param_dic['s']=0.2
for key,val in param_dic.items(): exec(key + '=val')
print(calcA(df, "CMat"))

param_dic['s']=0.8
for key,val in param_dic.items(): exec(key + '=val')
print(calcA(df, "CMat"))


# In[85]:


now = np.arange(1,10)
print(now)
print(now[6:])

# now=[0 for i in [1,2,3]]
a.shape


# # Calc p (drift rate)

# In[89]:



def calcp(df, name):
    
#     reveal_param(param_dic)
#     for key,val in param_dic.items(): exec(key + '=val')
    
    A = calcA(df, name)
    
    dnow = df[df['FileCondi']==name]

    p = np.repeat(3.0, A.shape[0])
    
    
    for i in np.arange(0,A.size):
        
        IR_old_current = eval("IR_old_"+dnow['Probtype'].iloc[i]+"_"+dnow['Oldnew'].iloc[i])
        IR_new_current = eval("IR_new_"+dnow['Probtype'].iloc[i]+"_"+dnow['Oldnew'].iloc[i])
#         print(IR_new_current,IR_old_current)
        
#         print((A[i] + IR_old_current)/(A[i] + IR_old_current + c + IR_new_current))
        p[i] = (A[i] + IR_old_current)/(A[i] + IR_old_current + c + IR_new_current)
#     
#     p = np.array([((A[i] + IR_old[pt.iloc[i]].iloc[dnow['Old'].iloc[i]-1])/(\
#                                             A[i] + IR_old[pt.iloc[dnow['Old'].iloc[i]-1]] +\
#                                                c + IR_new[pt.iloc[dnow['Old'].iloc[i]-1]])) for i in np.arange(0,A.size)])
    return(p)


param_dic['s']=0.2
for key,val in param_dic.items(): exec(key + '=val')
print(calcp(df, "CMat"))

param_dic['s']=0.8
for key,val in param_dic.items(): exec(key + '=val')
print(calcp(df, "MIX"))


# In[108]:


p[2]


# # Random Walk

# In[90]:




def Random_Walk(df, name):
    
#     reveal_param(param_dic)
#     for key,val in param_dic.items(): eval(key=val)

    p = calcp(df, name)
    
#     print("p is:",p)
#     print("new",New_crit)
    
    maxreps = 100
    Predi_resp = np.zeros(p.size)
    process_time = np.zeros(p.size)
    nstep = np.zeros(p.size)

    for i in range(p.size):

        start = time.time()
        Acc = 0
        for irep in range(maxreps):
            
            if p[i]>1: print("wrong p")
                
#             print(i,"trial;",irep,"step;","ACC",Acc,";p:",p[i])
            
            Acc += np.random.choice([1,-1], p= [p[i],1-p[i]])

            Predi_respi = np.random.choice([1,2]) # defult randomly assign if not found

            if Acc == Old_crit: 

                Predi_respi = 1 #'Old'
                break
            elif Acc == New_crit:

                Predi_respi = 2 #"new"
                break
        end = time.time()
        
        nstep[i] = irep
        process_time[i] = (end-start)*1000
        Predi_resp[i] = int(Predi_respi)


    return(process_time, Predi_resp, nstep)

# print(Random_Walk(df, "CMpure", param_dic))


# ## Theoretical Random Walk

# In[91]:


def calc_theoretical_RW(df, name):

    A = calcA(df, name)
    p = calcp(df, name)

    dnow = df[df['FileCondi']==name]
#     print(p.shape)

    IR_old_current = [eval("IR_old_"+dnow['Probtype'].iloc[i]+"_"+dnow['Oldnew'].iloc[i]) for i in range(p.size)]
    IR_new_current = [eval("IR_new_"+dnow['Probtype'].iloc[i]+"_"+dnow['Oldnew'].iloc[i]) for i in range(p.size)]
    
#     print(IR_old_current, IR_new_current)
    
    p_resp_old = np.zeros(p.size)
    pred_correct = np.zeros(p.size)

    q = 1-p

    AA = Old_crit
    BB = -New_crit

    exp_nstep = np.zeros(p.size)
    exp_nstep[p!=q] = [BB/(q[i]-p[i])-            ((BB+AA)/(q[i]-p[i]))*                        ((1-(q[i]/p[i])**BB)/(1-(q[i]/p[i])**(BB+AA)))                  for i in range(p.size)]
    exp_nstep[p==q] = AA* BB

    exp_duration = [t0 + 1/(A[i] + IR_old_current[i] + c + IR_new_current[i]) for i in range(p.size)]
    exp_raw_t = [(exp_nstep[i] * exp_duration[i]) for i in range(p.size)]

    p_resp_old[p!=q] = [ (1-(q[i]/p[i])**BB)/(1-(q[i]/p[i])**(AA+BB)) for i in range(p.size)]
    p_resp_old[p==q] = [BB/(AA+BB)]


    pred_correct[dnow["Old"]==1] = p_resp_old[dnow["Old"]==1]
    pred_correct[dnow["Old"]==2] = 1-p_resp_old[dnow["Old"]==2]
    
#     print(np.array(exp_raw_t).size)
    pred_rt = kao * np.array(exp_raw_t)

    return(np.array(pred_rt), np.array(pred_correct))

param_dic['s'] = 0.2
for key,val in param_dic.items(): exec(key + '=val')
# print(calcA(df, "MIX", param_dic))
calc_theoretical_RW(df, "MIX")

param_dic['s'] = 0.8
for key,val in param_dic.items(): exec(key + '=val')
calc_theoretical_RW(df, "MIX")
# print(calcA(df, "MIX", param_dic))
# pred_rt.shape


# In[43]:


something = np.array([ (1-(q[i]/p[i])**(New_crit+Old_crit)) for i in range(p.size)])
something[something==0]
something


# # Parameter set

# In[101]:


Numitem = 16
param_dic = {
       "alpha" : 0.2, 
       "beta" : 1.2,
       "s" : 0.2,
       "c" : 0.58,
       "LMT_old" : np.zeros(shape = Numitem),
       "LMT_new" : np.zeros(shape = Numitem),
       "Old_crit" : 3,
       "New_crit" : -3,
       "t0" : 500,
       "kao" : 38,
       "IR_old_CM_old" : 0.23,
       "IR_old_AN_old" : 0.14,
       "IR_old_VM_old": 0.4,

       "IR_old_CM_new" : 0,
       "IR_old_AN_new" : 0,
       "IR_old_VM_new": 0,

       "IR_new_CM_old" : 0,
       "IR_new_AN_old" : 0,
       "IR_new_VM_old" : 0,

       "IR_new_CM_new" : 0.62,
       "IR_new_AN_new" : 0,
       "IR_new_VM_new" : 0}

for key,val in param_dic.items(): exec(key + '=val')

# walk unit to time unit}
# pred_raw_rt, predi_resp = Random_Walk(df, "CMpure", param_dic)

# actual_resp = df[df['FileCondi']=='CMpure']['Response']
# predi_resp = predi_resp

# for key,val in param_dic.items(): eval(key + '=val')


# # Calculate WSSE

# In[93]:


wssd_weight = 2

def calc_wssd(df, name):
#     name = "MIX"
    #pred_raw_rt, pred_resp, nstep = Random_Walk(df, name, param_dic)
    pred_rt, pred_correct = calc_theoretical_RW(df, name)
#     print(pred_rt)

    # add columns to df
    df_sub = df[df["FileCondi"] == name]
    df_sub["pred_rt"] = pred_rt
    df_sub["pred_error"] = 1- pred_correct
    #1-((pred_resp == df_sub["Old"]).astype(int))

    # Calculate wwsd 
    df_sub_agg = df_sub.groupby(["Oldnew","Probtype"])[["RT","Error","pred_rt","pred_error"]].agg(["mean"])
#     print(df_sub.groupby(["Oldnew","Probtype"]).mean())
    
    WSSD = ((df_sub_agg["RT"]-df_sub_agg["pred_rt"])**2 +    wssd_weight*(df_sub_agg["Error"]-df_sub_agg["pred_error"])**2)/3

    return(WSSD)

calc_wssd(df, "MIX")


# # Parameter search

# ## Search s

# In[103]:


name = "MIX"
param_name = 's'
parm_range = [0.01, 0.99]
niter=5


# navg = 1
all_s = np.linspace(parm_range[0], parm_range[1],num = niter)

wsse_s = pd.DataFrame()

i = 0
# now = np.zeros((niter, calc_wssd(df, name).size))
for i_s in all_s:
    
    param_dic[param_name] = i_s
    for key,val in param_dic.items(): exec(key + '=val')

    wsse_s["iter"+str(i)] = calc_wssd(df, name)["mean"].to_numpy()
        
    i+=1
    
wsse_s

# wsse_s * 10**-7
i=0
for i in range(4):
    plt.figure()
    plt.plot(range(wsse_s.iloc[i].size), wsse_s.iloc[i])
    plt.show()
    


# In[44]:


minindex = wsse_s.iloc[0].to_numpy().argmin()
print(all_s[minindex])

param_dic[param_name] = all_s[minindex]
for key,val in param_dic.items(): exec(key + '=val')
param_name


# In[ ]:





# ## c

# In[104]:


name = "MIX"
param_name = 'c'
parm_range = [0.01, 0.3]
niter=5


# navg = 1
all_s = np.linspace(parm_range[0], parm_range[1],num = niter)

wsse_s = pd.DataFrame()

i = 0
# now = np.zeros((niter, calc_wssd(df, name).size))
for i_s in all_s:
    
    param_dic[param_name] = i_s
    for key,val in param_dic.items(): exec(key + '=val')

    wsse_s["iter"+str(i)] = calc_wssd(df, name)["mean"].to_numpy()
        
    i+=1
    
wsse_s

# wsse_s * 10**-7
i=0
for i in range(4):
    plt.figure()
    plt.plot(range(wsse_s.iloc[i].size), wsse_s.iloc[i])
    plt.show()
    


# In[54]:


minindex = wsse_s.iloc[0].to_numpy().argmin()
print(all_s[minindex])

param_dic[param_name] = all_s[minindex]
for key,val in param_dic.items(): exec(key + '=val')


# ## Old_crit

# In[105]:


name = "MIX"
param_name = 'Old_crit'
parm_range = [2, 10]
niter= parm_range[1]-parm_range[0]+1


# navg = 1
all_s = np.linspace(parm_range[0], parm_range[1],num = niter)

wsse_s = pd.DataFrame()

i = 0
# now = np.zeros((niter, calc_wssd(df, name).size))
for i_s in all_s:
    
    param_dic[param_name] = i_s
    for key,val in param_dic.items(): exec(key + '=val')

    wsse_s["iter"+str(i)] = calc_wssd(df, name)["mean"].to_numpy()
        
    i+=1
    
wsse_s

# wsse_s * 10**-7
i=0
for i in range(4):
    plt.figure()
    plt.plot(range(wsse_s.iloc[i].size), wsse_s.iloc[i])
    plt.show()
    


# In[65]:


minindex = wsse_s.iloc[0].to_numpy().argmin()
print(all_s[minindex])

param_dic[param_name] = all_s[minindex]
for key,val in param_dic.items(): exec(key + '=val')


# ## New_crit

# In[69]:


name = "MIX"
param_name = 'New_crit'
parm_range = [-10, -2]
niter= parm_range[1]-parm_range[0]+1


# navg = 1
all_s = np.linspace(parm_range[0], parm_range[1],num = niter)

wsse_s = pd.DataFrame()

i = 0
# now = np.zeros((niter, calc_wssd(df, name).size))
for i_s in all_s:
    
    param_dic[param_name] = i_s
    for key,val in param_dic.items(): exec(key + '=val')

    wsse_s["iter"+str(i)] = calc_wssd(df, name)["mean"].to_numpy()
        
    i+=1
    
wsse_s

# wsse_s * 10**-7
i=0
for i in range(4):
    plt.figure()
    plt.plot(range(wsse_s.iloc[i].size), wsse_s.iloc[i])
    plt.show()
    


# In[68]:


minindex = wsse_s.iloc[0].to_numpy().argmin()
print(all_s[minindex])

param_dic[param_name] = all_s[minindex]
for key,val in param_dic.items(): exec(key + '=val')


# ## t0

# In[71]:


name = "MIX"
param_name = 't0'
parm_range = [100, 1000]
niter= 100


# navg = 1
all_s = np.linspace(parm_range[0], parm_range[1],num = niter)

wsse_s = pd.DataFrame()

i = 0
# now = np.zeros((niter, calc_wssd(df, name).size))
for i_s in all_s:
    
    param_dic[param_name] = i_s
    for key,val in param_dic.items(): exec(key + '=val')

    wsse_s["iter"+str(i)] = calc_wssd(df, name)["mean"].to_numpy()
        
    i+=1
    
wsse_s

# wsse_s * 10**-7
i=0
for i in range(4):
    plt.figure()
    plt.plot(range(wsse_s.iloc[i].size), wsse_s.iloc[i])
    plt.show()
    


# In[ ]:





# In[ ]:





# In[ ]:




