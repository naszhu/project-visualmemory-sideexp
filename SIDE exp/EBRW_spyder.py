# -*- coding: utf-8 -*-
"""
Created on Wed Oct 30 20:11:31 2019

@author: naszh
"""

import numpy as np
import pandas as pd
import time
import matplotlib.pyplot as plt
import sys
import scipy.optimize as optimize
# =============================================================================
# 
# =============================================================================
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
df_crt = df_org[df_org["Error"]==0].\
groupby(["Oldnew","Setsize","Probtype","Lag","FileCondi"])[["RT"]].agg(["mean"]).reset_index()

df_err.columns = df_err.columns.droplevel(1)
df_crt.columns = df_crt.columns.droplevel(1)

# df_err.reindex(np.arange(1,df_err.shape[0]))
# df_crt.reindex(np.arange(1,df_crt.shape[0]))
# =============================================================================
# function 
# =============================================================================
global boost, alpha, beta, s, c ,Old_crit ,New_crit ,\
    t0 ,kao ,LTM
  
alpha={}; beta ={}; s = {}; LTM = {}


boost = 1.05
alpha['CM'] = 0.1969 
alpha['VM'] = 0.1969 
alpha["AN"]=0.84

beta['CM'] = 1.228
beta['VM'] = 1.228
beta["AN"]= 1.708

s["AN"] = 0.01
s["CM"] = 0.02
s["VM"] = 0.02
c = 0.3938

Old_crit = 1.9197
New_crit = -2.33
t0 = 699.98
kao = 37

LTM["CMold_iold"] = 0.463
LTM["VMold_iold"] = 0.143
LTM["ANold_iold"] = 0

LTM["CMnew_inew"] = 0.46
LTM["VMnew_inew"] = 0.039
LTM["ANnew_inew"] = 0

#the following is 0 all the times
LTM["CMold_inew"] = 0
LTM["VMold_inew"] = 0
LTM["ANold_inew"] = 0


LTM["CMnew_iold"] = 0
LTM["VMnew_iold"] = 0
LTM["ANnew_iold"] = 0





# =============================================================================
# 
# =============================================================================




def calcA(df, name):

    
    dnow = df[df['FileCondi']==name]
    a = np.zeros((dnow.shape[0], 8))  #activation
    
    if "CM" in df["Probtype"].to_numpy():
        snow=s["CM"]
        betanow=beta["CM"]
        alphanow=alpha["CM"]
        
    elif "VM" in df["Probtype"].to_numpy():
        snow=s["VM"]
        betanow=beta["VM"] 
        alphanow=alpha["VM"]
    
    m = np.array([(j**(-betanow) + alphanow) for j in np.arange(1,9)])
    m_an = np.array([(j**(-beta["AN"]) + alpha["AN"]) for j in np.arange(1,9)])
    

    
    for j in np.arange(1,9):
        
        indexj = j-1 
        
        a[np.logical_and(dnow['Lag'] != j, dnow['Probtype'] == "AN")\
          ,indexj]  = m_an[indexj] * s["AN"]   #AN, when i is not tested
        a[np.logical_and(dnow['Lag'] == j, dnow['Probtype'] == "AN")\
          ,indexj]  = m_an[indexj]    #AN, when i is tested
        
        a[np.logical_and(dnow['Lag'] != j, dnow['Probtype'] != "AN") ,indexj]  = m[indexj] *snow
        a[np.logical_and(dnow['Lag'] == j, dnow['Probtype'] != "AN") ,indexj]  = m[indexj] 
#     print("a=",a,"\n")
        
    for i in range(a.shape[0]): a[i,dnow['Setsize'].iloc[i]:] = 0 #a_ij suit for the correct amount of setsize
    
    debug = 0
    if debug==1:
        for i in range(a.shape[0]):
            print(i, "begin\n","a is",a[i].round(2),\
                  "\n m is", np.array(m).round(3),"\n Probtype is",\
                  "\n Lag is ",dnow['Lag'].iloc[i],\
                  "\n Setsize is", dnow["Setsize"].iloc[i],\
                  "\n Probtype:",\
                  dnow["Probtype"].iloc[i],\
                  "\n Oldnew: ",dnow["Oldnew"].iloc[i],\
                  "\n Ai is", a[i,:].sum().round(2),\
                  "\n m_an", s,\
                 "\n------------------------------------------------" )

    
    A = a.sum(axis = 1)
    
    return(A)

print(calcA(df[df["Error"]==0], "CMat"))



# =============================================================================
# 
# =============================================================================

def calcp(df, name):
    
    A = calcA(df, name)
    dnow = df[df['FileCondi']==name]
    p = np.repeat(3.0, A.shape[0])
  
    for i in np.arange(0,A.size):
        
        #LTM["CMold_iold"]
        
        IR_old_current = LTM[dnow['Probtype'].iloc[i] + "old_i"+ dnow['Oldnew'].iloc[i]]
        IR_new_current = LTM[dnow['Probtype'].iloc[i] + "new_i"+ dnow['Oldnew'].iloc[i]]
        
#         print(IR_old_current,IR_new_current)
        
        p[i] = (A[i] + IR_old_current)/(A[i] + IR_old_current + c + IR_new_current)

    return(np.array(p))



ok1=calcp(df, "MIX")
on = df[df['FileCondi']=="MIX"]["Oldnew"]
onn = df[df['FileCondi']=="MIX"]["Probtype"]
# print(np.array([str(ok1[i].round(2))+"-"+str(on.iloc[i]) for i in range(ok1.shape[0])]))
print(np.array([str(ok1[i].round(2))+"-"+str(on.iloc[i]) +"-"+str(onn.iloc[i]) for i in range(ok1.shape[0])]))




# =============================================================================
# 
# =============================================================================


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

# =============================================================================
# 
# =============================================================================

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

    df_sub_crt["wSSD_RT"] = [df_sub_crt.loc[i,"SSD_RT"] *\
                             w[df_sub_crt.loc[i,"Oldnew"]+"_rt"] for i in range(df_sub_crt.shape[0])]

    df_sub_err["wSSD_err"] = [df_sub_err.loc[i,"SSD_err"] *\
                             w[df_sub_err.loc[i,"Oldnew"] + "_err"] for i in range(df_sub_err.shape[0])]

    WSSD = df_sub_err["wSSD_err"].sum() + df_sub_crt["wSSD_RT"].sum()

    return(WSSD)

time1=time.time()
calc_wssd(df_org, "MIX")
print(time.time()-time1)
# df_sub_agg, WSSD = calc_wssd(df, "MIX")


# =============================================================================
# 
# =============================================================================

param_dic=np.array((boost ,
        alpha['CM'] ,
        alpha['VM'] ,
        alpha["AN"] ,

        beta['CM'] ,
        beta['VM'] ,
        beta["AN"] ,

        s["AN"] ,
        s["CM"] ,
        s["VM"] ,
        c ,

        Old_crit ,
        New_crit ,
        t0,
        kao,

        LTM["CMold_iold"] ,
        LTM["VMold_iold"] ,
        LTM["ANold_iold"] ,
        
        LTM["CMnew_inew"] ,
        LTM["VMnew_inew"] ,
        LTM["ANnew_inew"] ,

        #the following is 0 all the times
        LTM["CMold_inew"] ,
        LTM["VMold_inew"] ,
        LTM["ANold_inew"] ,


        LTM["CMnew_iold"] ,
        LTM["VMnew_iold"] ,
        LTM["ANnew_iold"]))
param_dic

# =============================================================================
# optim_wsse
# =============================================================================
def optim_wsse(params_dic):
    global boost, alpha, beta, s, c ,Old_crit ,New_crit ,\
        t0 ,kao ,LTM

    alpha={}; beta ={}; s = {}; LTM = {}
    
    # print(params)  # <-- you'll see that params is a NumPy array
    [
        boost ,
        alpha['CM'] ,
        alpha['VM'] ,
        alpha["AN"] ,

        beta['CM'] ,
        beta['VM'] ,
        beta["AN"] ,

        s["AN"] ,
        s["CM"] ,
        s["VM"] ,
        c ,

        Old_crit ,
        New_crit ,
        t0,
        kao,

        LTM["CMold_iold"] ,
        LTM["VMold_iold"] ,
        LTM["ANold_iold"] ,
        
        LTM["CMnew_inew"] ,
        LTM["VMnew_inew"] ,
        LTM["ANnew_inew"] ,

        #the following is 0 all the times
        LTM["CMold_inew"] ,
        LTM["VMold_inew"] ,
        LTM["ANold_inew"] ,


        LTM["CMnew_iold"] ,
        LTM["VMnew_iold"] ,
        LTM["ANnew_iold"]] = params_dic
    


    return calc_wssd(df_org, "MIX")

# =============================================================================
# search
# =============================================================================
bdd = ((0,0), #boost
       (0.00001, 1.0), #alphacm
       (0, 0), #alphavm
       (0.00001, 1.0), #alpha _an
       
       (1, 3.0), #beta_cm
       (0, 0), #beta_vm
       (1, 3.0), #beta_an
       
       (0.0001, 1.0), #s an
       (0.0001, 0.99), #s cm
       (0,0), #s vm
       
       (0.00001, 0.99), #c
       
       (1, 5.0), #old_cr
       (-5.0, -1), #new_cr
       (600, 700), #t0
       (0.00001, 50), #kao
       
       (0.00001, 1.0), #IR_old_CM_old
       (0, 0), #old_vm_old
       (0, 0), #old_an_old
       (0, 1.0), #new_cm_new
       (0, 0), #new_vm_new
       (0, 0.5), #new_an_new
       (0,0), #old_cm_new
       (0,0), #old_vm_new
       (0,0), #old_an_new
       (0,0), #new_cm_old
       (0,0), #new_vm_old
       (0,0) #new_an_old
      )





initial_guess = param_dic
time1= time.time()
result = optimize.minimize(optim_wsse, initial_guess, bounds=bdd,options={"maxiter" : 1500})
if result.success:
    fitted_params = result.x
    print(fitted_params)
else:
    raise ValueError(result.message)
time.time()-time1
# params_dic
# alpha_AN


