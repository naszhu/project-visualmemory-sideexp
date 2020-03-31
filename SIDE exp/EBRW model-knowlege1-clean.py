#!/usr/bin/env python
# coding: utf-8

# # Set up 
# 

# In[1]:


import numpy as np
import pandas as pd
import time
import matplotlib.pyplot as plt
import sys
import scipy.optimize as optimize
import threading
from scipy.optimize import LinearConstraint
from datetime import datetime


# In[2]:


path="./" + "sidedata/"
filename = path + "Alldata2.csv"

df=pd.read_csv(filename, index_col=None)
df['Probtype'] = np.where(df['Stimkind']==1, "CM",
                   np.where(df['Stimkind']==0, "AN",
                   np.where(df['Stimkind']==3, 'VM',"wrong")))

df['Oldnew'] = np.where(df['Old']==1, "old",
                   np.where(df['Old']==2, "new","wrong"))
df['Error'] = 1-df['Correctness']
df_org = df[df['RT']<2000]
df = df_org.groupby(["Oldnew","Setsize","Probtype","Lag","Error","FileCondi"])[["RT"]].agg(["mean"])
# df[["Correctness","Error"]]
# df.groupby("Error")
df.index.name = 'Setsize'
df.reset_index(inplace=True)
# df["Error"]


df_err = df_org.groupby(["Oldnew","Setsize","Probtype","Lag","FileCondi"])[["Error"]].agg(["mean"]).reset_index()
df_crt = df_org[df_org["Error"]==0].groupby(["Oldnew","Setsize","Probtype","Lag","FileCondi"])[["RT"]].agg(["mean"]).reset_index()

df_err.columns = df_err.columns.droplevel(1)
df_crt.columns = df_crt.columns.droplevel(1)

# df_err.reindex(np.arange(1,df_err.shape[0]))
# df_crt.reindex(np.arange(1,df_crt.shape[0]))


# In[3]:


df


# ## Some global setting

# In[28]:


global vary_ss
vary_ss=1


# # Parameter set

# In[31]:


global boost, alpha, beta, s, c ,Old_crit ,New_crit ,    t0 ,kao , F, L
  
alpha={}; beta ={}; s = {}; F = {}; L = {}


boost = 1.05
alpha["all"] = 0.19
# alpha['CM'] = 0.1969 
# alpha['VM'] = 0.1969 
# alpha["AN"]=0.84

beta["all"] =1.2
# beta['CM'] = 1.228
# beta['VM'] = 1.228
# beta["AN"]= 1.708

s["all"]=0.01
s["ss24"]=0.01
s["ss8"]=0.01
# s["AN"] = 0.01
# s["CM"] = 0.02
# s["VM"] = 0.02
c = 0.3938

Old_crit = 1.9197
New_crit = -2.33
t0 = 699.98
kao = 37

#---CM
F["CM_oldiold_oldinew"] = 0.2

L["CM_oldiold_newinew"] = 0.2
L["CM_oldinew_newiold"] = 0.2

#---AN
F["AN_oldiold_oldinew"] = 0.2


# In[ ]:





# In[ ]:





# # Assign LTM global

# In[5]:


def assign_LTM_global(item_condi, walk, item):
    
    global F,L,Fnow,Lnow
    
    if item_condi=="CM":
        
        if walk+"i"+item == "oldiold" or walk+"i"+item == "oldinew":
            
            Fnow = F["CM_oldiold_oldinew"]
        else: Fnow = 0
        
        if walk+"i"+item == "oldiold" or walk+"i"+item == "newinew":
            
            Lnow = L["CM_oldiold_newinew"]
        elif walk+"i"+item == "oldinew" or walk+"i"+item == "newiold":
            
            Lnow = L["CM_oldinew_newiold"]
        else: Lnow=0
            
    if item_condi == "AN":
        
        if walk+"i"+item == "oldiold" or walk+"i"+item == "oldinew":
            
            Fnow = F["AN_oldiold_oldinew"]
        else: Fnow=0
        
        Lnow = 0
    
    return(Fnow + Lnow)

assign_LTM_global("CM","new","old")
        


# # Calc A (activation)

# In[41]:


# seperate albes
def calcA(df, name):

    
    dnow = df[df['FileCondi']==name]
    a = np.zeros((dnow.shape[0], 8))  #activation
    

    betanow=beta["all"] 
    alphanow=alpha["all"]
    
    m = np.array([(j**(-betanow) + alphanow) for j in np.arange(1,9)])
    
    for j in range(1,9):

        indexj = j-1
        
        
        if vary_ss==0:
            
            a[dnow["Lag"]==j,indexj] = m[indexj]
            a[dnow["Lag"]!=j,indexj] = m[indexj] * s["all"]
        else:
            a[np.logical_and(dnow["Lag"]!=j,dnow["Setsize"]!=8),indexj] = m[indexj] * s["ss24"]
            a[np.logical_and(dnow["Lag"]!=j,dnow["Setsize"]==8),indexj] = m[indexj] * s["ss8"]
            a[dnow["Lag"]==j,indexj] = m[indexj]
        
        
#         print(m[indexj] * s["all"])
    for i in range(a.shape[0]): a[i,dnow['Setsize'].iloc[i]:] = 0 #a_ij suit for the correct amount of setsize
    debug = 0
    if debug==1:
        for i in range(a.shape[0]):
            print(i, "begin\n","a is",a[i].round(5),                  "\n m is", np.array(m).round(3),"\n Probtype is",                  "\n Lag is ",dnow['Lag'].iloc[i],                  "\n Setsize is", dnow["Setsize"].iloc[i],                  "\n Probtype:",                  dnow["Probtype"].iloc[i],                  "\n Oldnew: ",dnow["Oldnew"].iloc[i],                  "\n Ai is", a[i,:].sum().round(2),                  "\n snow", snow,                 "\n------------------------------------------------" )

    
    A = a.sum(axis = 1)
    
    return(A)

print(calcA(df[df["Error"]==0], "CMat"))


# # Calc p (drift rate)

# In[7]:



def calcp(df, name):
    
    A = calcA(df, name)
    dnow = df[df['FileCondi']==name]
    p = np.repeat(3.0, A.shape[0])
    
    probs = dnow["Probtype"].astype("category").cat.categories.to_numpy() #get categories name
    

    for iprob in probs:
        
        for ion in ["old","new"]:
            
            tf_ion = dnow['Oldnew'] == ion
            tf_iprob = dnow['Probtype'] == iprob
            tf_all = np.logical_and(tf_ion, tf_iprob)
            
            IR_old_current = assign_LTM_global(iprob,"old",ion)
            IR_new_current = assign_LTM_global(iprob,"new",ion)
        
            p[tf_all] = (A[tf_all] + IR_old_current)/(A[tf_all] + IR_old_current + c + IR_new_current)
        

    return(np.array(p))



ok1=calcp(df, "MIX")
on = df[df['FileCondi']=="MIX"]["Oldnew"]
onn = df[df['FileCondi']=="MIX"]["Probtype"]
ss = df[df['FileCondi']=="MIX"]["Setsize"]
lagg = df[df['FileCondi']=="MIX"]["Lag"]
# print(np.array([str(ok1[i].round(2))+"-"+str(on.iloc[i]) for i in range(ok1.shape[0])]))
print(np.array([str(ok1[i].round(8))+"-"+str(on.iloc[i]) +"-"+str(onn.iloc[i])+"-ss:"+str(ss.iloc[i])                + "-lag:" + str(lagg.iloc[i]) for i in range(ok1.shape[0])]))


# # RW 

# In[8]:



def calc_theoretical_RW(df, name):

    A = calcA(df, name)
    p = calcp(df, name)

#     dnow
    dnow = df[df['FileCondi']==name]

    p_resp_old = np.zeros(p.size)
    pred_correct = np.zeros(p.size)

    q = 1-p

    AA = Old_crit
    BB = -New_crit # bb is a postive number

    exp_nstep = np.zeros(p.shape[0])

    
    theta1 = ((p/q)**(AA+BB)+1)/((p/q)**(AA+BB)-1)
    theta2 = ((p/q)**BB+1)/((p/q)**BB-1)

    
    tf = np.logical_and(p!=q , dnow["Oldnew"]=='old')
    exp_nstep[tf] = (1/(p[tf]-q[tf]))*(theta1[tf]*(AA+BB) - theta2[tf]*BB)
    
    tf = np.logical_and(p==q , dnow["Oldnew"]=='old')
    exp_nstep[tf] = (AA/3)*(2*BB+AA)
    
    #-------------------
    theta1 = ((p/q)**(-(AA+BB))+1)/((p/q)**(-(AA+BB))-1)
    theta2 = ((p/q)**-AA+1)/((p/q)**-AA-1)

    tf = np.logical_and(p!=q , dnow["Oldnew"]=='new')
    exp_nstep[tf] = (1/(q[tf]-p[tf]))*(theta1[tf]*(AA+BB) - theta2[tf]*AA)
    tf = np.logical_and(p==q , dnow["Oldnew"]=='new')
    exp_nstep[tf] = (BB/3)*(2*AA+BB)
    
    
    p_resp_old[p!=q] = ((1-(q/p)**BB)/(1-(q/p)**(AA+BB)))[p!=q]
    p_resp_old[p==q] = BB/(AA+BB)
    
    pred_correct[dnow["Oldnew"]=='old'] = p_resp_old[dnow["Oldnew"]=="old"]
    pred_correct[dnow["Oldnew"]=="new"] = 1-p_resp_old[dnow["Oldnew"]=="new"]
    
    pred_rt = t0 + kao * exp_nstep
    


    return(np.array(pred_rt), np.array(pred_correct))
          

df_sub_err = df_org.groupby(["Oldnew","Setsize","Probtype","Lag","FileCondi"])[["Error"]].agg(["mean"]).reset_index()


ok1=calcp(df_sub_err, "MIX")
on = df_sub_err[df_sub_err['FileCondi']=="MIX"]["Oldnew"]
onn = df_sub_err[df_sub_err['FileCondi']=="MIX"]["Probtype"]
# print(np.array([str(ok1[i].round(2))+"-"+str(on.iloc[i]) +"-"+str(onn.iloc[i]) for i in range(ok1.shape[0])]))

jjj, ok1 =calc_theoretical_RW(df_sub_err, "MIX")
on = df_sub_err[df_sub_err['FileCondi']=="MIX"]["Oldnew"]
onn = df_sub_err[df_sub_err['FileCondi']=="MIX"]["Probtype"]
one = (df_sub_err[df_sub_err['FileCondi']=="MIX"][[("Error","mean")]]).to_numpy().round(3)

# print(np.array([str((ok1[i]).round(2))+"-"+str(on.iloc[i])+"-"+str(onn.iloc[i]) for i in range(ok1.shape[0])]))

np.array([str(one[i])+"-"+str(on.iloc[i])+"-"+str(onn.iloc[i]) for i in range(ok1.shape[0])])


# # Calculate WSSE

# In[9]:


w={}
w["new_rt"] = 4*2
w["old_rt"] = 1*2
w["new_err"] = 4
w["old_err"] = 1

def calc_wssd(df, name):

    df_sub_org = df[df["FileCondi"] == name]
    df_sub_err = df_err[df_err["FileCondi"] == name]
    df_sub_crt = df_crt[df_crt["FileCondi"] == name]

    df_sub_err = df_sub_err.reset_index()
    df_sub_crt = df_sub_crt.reset_index()

    pred_rt, pred_correct = calc_theoretical_RW(df_sub_err, name)
    pred_crt, pred_ccorrect = calc_theoretical_RW(df_sub_crt, name)

    pred_rt, pred_correct = calc_theoretical_RW(df_sub_err, name)
    pred_crt, pred_ccorrect = calc_theoretical_RW(df_sub_crt, name)

    df_sub_crt["pred_crt"] = pred_crt/1000 #translate crt to seconds. 
    df_sub_err["pred_error"] = 1- pred_correct

    df_sub_crt["RT"]= df_sub_crt["RT"]/1000

    df_sub_crt["SSD_RT"] = ((df_sub_crt["RT"] - df_sub_crt["pred_crt"])**2).to_numpy()
    df_sub_err["SSD_err"] = ((df_sub_err["Error"] - df_sub_err["pred_error"])**2).to_numpy()

    df_sub_crt["wSSD_RT"] = [df_sub_crt.loc[i,"SSD_RT"] *                             w[df_sub_crt.loc[i,"Oldnew"]+"_rt"] for i in range(df_sub_crt.shape[0])]

    df_sub_err["wSSD_err"] = [df_sub_err.loc[i,"SSD_err"] *                             w[df_sub_err.loc[i,"Oldnew"] + "_err"] for i in range(df_sub_err.shape[0])]

    WSSD = df_sub_err["wSSD_err"].sum() + df_sub_crt["wSSD_RT"].sum()

    return(WSSD)

# time1=time.time()
# calc_wssd(df_org, "MIX2")
# print(time.time()-time1)
# df_sub_agg, WSSD = calc_wssd(df, "MIX")


# # Parameter search

# ## Asign parm_dic

# In[42]:


if vary_ss==0:
     param_dic=np.array(( alpha["all"], beta["all"], s["all"],
                         c, Old_crit, New_crit, t0, kao,
                         F["CM_oldiold_oldinew"], L["CM_oldiold_newinew"], L["CM_oldinew_newiold"],
                         F["AN_oldiold_oldinew"]))
else:
     param_dic=np.array((alpha["all"], beta["all"], s["ss24"], s["ss8"],
                         c, Old_crit, New_crit, t0, kao,
                         F["CM_oldiold_oldinew"], L["CM_oldiold_newinew"], L["CM_oldinew_newiold"],
                         F["AN_oldiold_oldinew"]))
         
param_dic


# ## optim_wsse

# In[43]:


def optim_wsse(params_dic):
    global alpha, beta, s, c ,Old_crit ,New_crit ,        t0 ,kao ,F, L

    alpha={}; beta ={}; s = {}; F = {}; L={}
    
    # print(params)  # <-- you'll see that params is a NumPy array
    if vary_ss==0:
        [alpha["all"], beta["all"], s["all"],
                             c, Old_crit, New_crit, t0, kao,
                             F["CM_oldiold_oldinew"], L["CM_oldiold_newinew"], L["CM_oldinew_newiold"],
                             F["AN_oldiold_oldinew"]] = params_dic
    else:
        [alpha["all"], beta["all"], s["ss24"], s["ss8"],
                             c, Old_crit, New_crit, t0, kao,
                             F["CM_oldiold_oldinew"], L["CM_oldiold_newinew"], L["CM_oldinew_newiold"],
                         F["AN_oldiold_oldinew"]] = params_dic
    


    return calc_wssd(df_org, "MIX2")
optim_wsse(param_dic)


# ## random start

# In[55]:


def random_start():
    global alpha, beta, s, c ,Old_crit ,New_crit ,    t0 ,kao , F, L
  
    alpha={}; beta ={}; s = {}; F = {}; L = {}

    alpha["all"] = np.random.uniform(0.1,1)
    beta["all"] =np.random.uniform(1,2)
    s["all"] = np.random.uniform(0.01,1)
    s["ss24"] = np.random.uniform(0.01,1)
    s["ss8"] = np.random.uniform(0.01,1)
    c = np.random.uniform(0.05,0.99)

    Old_crit = np.random.uniform(1,5)
    New_crit = np.random.uniform(-5,-1)
    t0 = 445.5
    kao = 50.4

    #---CM
    F["CM_oldiold_oldinew"] = np.random.uniform(0.01,1)

    L["CM_oldiold_newinew"] = np.random.uniform(0.01,1)
    L["CM_oldinew_newiold"] = np.random.uniform(0.001,1)

    #---AN
    F["AN_oldiold_oldinew"] = np.random.uniform(0.001,1)
    if vary_ss==0:
        param_dic=[alpha["all"], beta["all"], s["all"],
                             c, Old_crit, New_crit, t0, kao,
                             F["CM_oldiold_oldinew"], L["CM_oldiold_newinew"], L["CM_oldinew_newiold"],
                             F["AN_oldiold_oldinew"]] 
    else:
        param_dic=[alpha["all"], beta["all"], s["ss24"], s["ss8"],
                             c, Old_crit, New_crit, t0, kao,
                             F["CM_oldiold_oldinew"], L["CM_oldiold_newinew"], L["CM_oldinew_newiold"],
                         F["AN_oldiold_oldinew"]]

    return(param_dic)

if vary_ss==0:
    bdd = (

        (0.00001, 1.0), #alphaall
        (1, 3.0), #beta_all
        (0.001, 1.0), #s all

        (0.001, 0.99), #c

        (1, 10.0), #old_cr
        (-10.0, -1), #new_cr
        (500,500),#(400,700), #t0
        (22,22),#(10, 70), #kao

        (0,0.99),#F["CM_oldiold_oldinew"],

        (0,0.99),#L["CM_oldiold_newinew"],
        (0,0.99),#L["CM_oldinew_newiold"],

        #---AN
        (0,0.50)#F["AN_oldiold_oldinew"]
          )
else:
    bdd = (

        (0.00001, 1.0), #alphaall
        (1, 3.0), #beta_all
        (0.001, 1.0), #s 24
        (0.001, 1.0), #s 8

        (0.001, 0.99), #c

        (1, 10.0), #old_cr
        (-10.0, -1), #new_cr
        (500,500),#(400,700), #t0
        (22,22),#(10, 70), #kao

        (0,0.99),#F["CM_oldiold_oldinew"],

        (0,0.99),#L["CM_oldiold_newinew"],
        (0,0.99),#L["CM_oldinew_newiold"],

        #---AN
        (0,0.50)#F["AN_oldiold_oldinew"]
          )
    

def Jcstrain():
    if vary_ss==0:
        linear_constraint = LinearConstraint(            [[ 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, -1,  0],[ 0, 0, 0, 0, 0, 0, 0, 0, 1, 0,  0, -1]],            [0.01, 0.01],[np.inf,np.inf])
    else: 
        linear_constraint = LinearConstraint(            [[ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, -1,  0],[0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0,  0, -1]],            [0.01, 0.01],[np.inf,np.inf])
    
    return(linear_constraint)


# ## Actual search

# In[56]:


def actual_search_trust3():# 

    time1 = time.time()
    numit=5
    global fitted_params3,result
    fitted_params3 = np.zeros((numit,12))
    for i in range(0,numit):
        param_dic = random_start()
        result = optimize.minimize(optim_wsse, param_dic,                               bounds=bdd,options={'verbose': 1,"maxiter" : 1500},
                                   constraints=[Jcstrain(), ],
                                   method='trust-constr')
        print(i)
        if result.success:
            fitted_params3[i,:] = result.x
        else: 
            print(result.message)


    time2 = time.time()
    print(time2-time1)


# In[ ]:


actual_search_trust3()

def write_file(which_work):
    
    if vary_ss==0:
        Parnames = ["alpha","beta","s","c","old_crt","new_crt","t0","kappa","F_cmooon","L_cmoonn","L_cmonno","F_anooon"]
    else:
        Parnames = ["alpha","beta","s_ss24","s_ss8","c","old_crt","new_crt","t0","kappa","F_cmooon","L_cmoonn","L_cmonno","F_anooon"]
    fitted_params3
    fitpdf=pd.DataFrame(fitted_params3)
    fitpdf.columns=Parnames
    now = datetime.now()
    current_time = now.strftime("%m_%d_%H_%M")
#     print("Current Time =", current_time)
    fitsel=fitpdf[np.logical_and(fitpdf["alpha"]>0 , fitpdf["alpha"]<1)]
    #------
#     print(bdd)
    bddstr=np.array(bdd).astype(str)
    bddbdd=pd.DataFrame([[bddstr[i,0]+" ~ "+bddstr[i,1] for i in range(bddstr.shape[0])]])
    bddbdd.columns=Parnames
    fitsel=fitsel.append(bddbdd, sort=False)
    # print(fitsel)
    #---------
    
    filen = which_work+"_" + current_time +".csv"
    
    fitsel.to_csv(filen)

    return(fitsel)
    
global Which_file
Which_file = "MIX2"
current = "s_vary_ss"

fitsel=write_file(Which_file + "_" + current)
fitsel


# In[ ]:


fitsel


# ## Assignresults

# In[21]:


[
       alpha["all"],
       beta["all"],
       s["all"],
       c,
       Old_crit,
       New_crit,
       t0,
       kao,
       F["CM_oldiold_oldinew"],
       L["CM_oldiold_newinew"],
       L["CM_oldinew_newiold"],
       F["AN_oldiold_oldinew"]] = fitted_params3[0]


# In[24]:


[
       alpha["all"],
       beta["all"],
       s["all"],
       c,
       Old_crit,
       New_crit,
       t0,
       kao,
       F["CM_oldiold_oldinew"],
       L["CM_oldiold_newinew"],
       L["CM_oldinew_newiold"],
       F["AN_oldiold_oldinew"]]


# # Final Random Walk

# In[ ]:


[alpha["all"],
beta["all"],
s["all"],
c,
Old_crit,
New_crit,
t0,
kao,
F["CM_oldiold_oldinew"],
L["CM_oldiold_newinew"],
L["CM_oldinew_newiold"],
F["AN_oldiold_oldinew"]]=
   [0.24941660611987385,
2.2862066094582514,
0.029145705532782337,
0.23158409905186297,
3.231419242494317,
-3.386426204534701,
638.0157651219173,
37.139415940933524,
0.040672079639333676,
0.027151542160548207,
0.017151333293794205,
0.019855019723203894]


# In[18]:


s_AN=0.01
s_others = 0.07
def finalRW(df_org):
    names = df_org["FileCondi"].astype("category").cat.categories.to_numpy() #get categories name
    df_all_crt= pd.DataFrame()
    df_all_err= pd.DataFrame()

    df_sub_err = df_org.groupby(["Oldnew","Setsize","Probtype","Lag","FileCondi"])[["Error"]].agg(["mean"]).reset_index()
    df_sub_crt = df_org[df_org["Error"]==0].    groupby(["Oldnew","Setsize","Probtype","Lag","FileCondi"])[["RT"]].agg(["mean"]).reset_index()

    for iname in names:

        pred_rt, pred_correct = calc_theoretical_RW(df_sub_err, iname)
        pred_crt, pred_ccorrect = calc_theoretical_RW(df_sub_crt, iname)

        df_temp_crt = df_sub_crt[df_sub_crt["FileCondi"] == iname]
        df_temp_err = df_sub_err[df_sub_err["FileCondi"] == iname]


        df_temp_crt["pred_rt"] = pred_rt   
        df_temp_err["pred_error"] = 1- pred_correct

        df_all_crt = df_all_crt.append(df_temp_crt)
        df_all_err = df_all_err.append(df_temp_err)


    return(df_all_crt, df_all_err)


# ## Error plot

# In[23]:



# print(L["CM_oldinew_newiold"])
df_all_crt, df_all_err = finalRW(df_org)

def all_plot(df_all_aggnew,plotwhaty,ylim):
    fig, axes = plt.subplots(nrows=1, ncols=6, figsize=(16,8))
    fig.canvas.set_window_title(plotwhaty)
    line_width = 2.5
    i=0
    names = df_all_aggnew["FileCondi"].astype("category").cat.categories.to_numpy() #get categories name
    for iname in names:
        df_all_aggnew[df_all_aggnew["FileCondi"]==iname].plot.line(x="Setsize", y=plotwhaty,ax=axes[i],                                                                   title=iname,ylim=ylim )
        i+=1

# all_plot(df_all_aggnew, "RT",(600,1100))
# all_plot(df_all_aggnew, "pred_rt",(600,1100))
# all_plot(df_all_aggnew, "Error",(0,0.3))
# all_plot(df_all_aggnew, "pred_error",(0,1))

df_all_agg = df_all_err.groupby(["Oldnew","Setsize","Probtype","FileCondi"]).agg(["mean"])#.apply(lambda x: x)
df_all_aggnew = df_all_agg.unstack(["Oldnew","Probtype"])
# df_all_aggnew = df_all_agg.unstack(["FileCondi","Oldnew","Probtype"])
df_all_aggnew
# df_all_aggnew= df_all_aggnew.swaplevel(i=0,j=2,axis = 1)

df_all_aggnew.index.name = 'Setsize'
df_all_aggnew.reset_index(inplace=True)


# df_all_aggnew.columns = df_all_aggnew[df_all_aggnew["FileCondi"]==iname].columns.droplevel([1,2])
df_all_aggnew.columns = df_all_aggnew.columns.droplevel([1,2])
# df_all_aggnew


all_plot(df_all_aggnew, "pred_error",(0,0.3))


# In[56]:


all_plot(df_all_aggnew, "Error",(0,0.3))


# ## correct RT plot

# In[24]:


def all_plot(df_all_aggnew,plotwhaty,ylim):
    fig, axes = plt.subplots(nrows=1, ncols=6, figsize=(16,8))
    fig.canvas.set_window_title(plotwhaty)
    line_width = 2.5
    i=0
    names = df_all_aggnew["FileCondi"].astype("category").cat.categories.to_numpy() #get categories name
    for iname in names:
        df_all_aggnew[df_all_aggnew["FileCondi"]==iname].plot.line(x="Setsize", y=plotwhaty,ax=axes[i],                                                                   title=iname,ylim=ylim )
        i+=1

# all_plot(df_all_aggnew, "RT",(600,1100))
# all_plot(df_all_aggnew, "pred_rt",(600,1100))
# all_plot(df_all_aggnew, "Error",(0,0.3))
# all_plot(df_all_aggnew, "pred_error",(0,1))

df_all_agg = df_all_crt.groupby(["Oldnew","Setsize","Probtype","FileCondi"]).agg(["mean"])#.apply(lambda x: x)
df_all_aggnew = df_all_agg.unstack(["Oldnew","Probtype"])
# df_all_aggnew = df_all_agg.unstack(["FileCondi","Oldnew","Probtype"])
df_all_aggnew
# df_all_aggnew= df_all_aggnew.swaplevel(i=0,j=2,axis = 1)

df_all_aggnew.index.name = 'Setsize'
df_all_aggnew.reset_index(inplace=True)


df_all_aggnew.columns = df_all_aggnew.columns.droplevel([1,2])
# df_all_aggnew

all_plot(df_all_aggnew, "RT",(600,1100))
all_plot(df_all_aggnew, "pred_rt",(600,1100))


# In[ ]:




