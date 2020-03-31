# -*- coding: utf-8 -*-
"""
Created on Tue Aug 28 19:28:06 2018

@author: Shuchun Lai
"""


#3-7 3 each . Oswald et al. The development of a shor domain-general measure of WMC

from genOperationList import genOperationList
from genLetterList import genLetterList
from getExpInfo import getExpInfo
from genStimulusList import genStimulusList
import numpy as np
import random
from psychopy import locale_setup, sound, gui, visual, core, data, event, logging, clock
from psychopy.constants import (NOT_STARTED, STARTED, PLAYING, PAUSED,
                                STOPPED, FINISHED, PRESSED, RELEASED, FOREVER)


import csv
import os

# Get experimental parameters:
expInfo = getExpInfo

expName = u'updating'  # from the Builder filename that created this script
expInfo_now = {'participant':'', 'session':'001'}
dlg = gui.DlgFromDict(dictionary=expInfo_now, title=expName)
if dlg.OK == False:
    core.quit()  # user pressed cancel

endExpNow = False

win = visual.Window([2000,1000], fullscr=True, monitor="testMonitor", units="deg", color=[0,0,0])
win.setMouseVisible(False)
path=os.getcwd()
#----------------------------------------------------------------------------------------

def Message(mes):
    txt=visual.TextStim(win=win, name='text_3',
    text=mes,
    font=u'Arial',
    pos=(0, 0), #height=0.1, 
    color=u'white')
    
    txt.draw()
    win.flip(clearBuffer=True)
    


def doTrial(windowPtr, expInfo, letterStimulus, operationStimulus):
    from getExpInfo import getExpInfo
    expInfo = getExpInfo
  
    
    #itemLength = len(letterStimulus) #from input of arguement
    print("  set size: "+str(len(letterStimulus)))
    for i in range(len(letterStimulus)):
        Message(operationStimulus[i]+'  '+letterStimulus[i])
        print("    "+operationStimulus[i]+" "+letterStimulus[i])
        
        event.waitKeys(keyList=["space"])
        #core.wait(expInfo.letterDuration)
        
    
    Message('Answer:')
    event.waitKeys(keyList=["space"])
    
#*-----------------------------------------------------------------------------------
#do trial
Message(expInfo.welcomeMsg)
event.waitKeys(keyList=["space"])

# Generate stimulus list for practice:
#pracStimulusList = genStimulusList(expInfo)
letter,operation,op_ans=genStimulusList(expInfo)

Message(expInfo.pracMsg)
event.waitKeys(keyList=["space"])

# Do practice:
for i in range(expInfo.nPracTrials):
    print("prac trial "+str(i)+":   ")
    doTrial(win,expInfo,letter[i],operation[i])
    
#------------------------------------------------letter and operation list is here. 

#===========================================================================================
#do EXP
# Generate stimulus list for experiment:
letter2,operation2,op_ans2 = genStimulusList(expInfo)

Message(expInfo.expMsg)
event.waitKeys(keyList=["space"])

for i in range(expInfo.nTrials):
    print("trial "+str(i+1)+":   ")
    doTrial(win,expInfo,letter2[i],operation2[i])

#------------------------------------------
#write file

filePath = path + os.sep + u'data/%s_%s' % ('participant', expInfo_now['participant'])
myFile=open(filePath+'.csv','wb')


block_num_write=np.concatenate( (np.repeat(0,3),np.repeat(1,12)) ) #num of prac+num_trial
trial_num_write=np.concatenate( (np.arange(1,4),np.arange(1,13)) )
op_ans_write=np.concatenate( (op_ans[0:expInfo.nPracTrials],op_ans2) )
letter_write=np.concatenate( (letter[0:expInfo.nPracTrials],letter2) ) #EX: ([A,F,G],[G,G,G,S],...[Z,B,D,D])
op_write=np.concatenate( (operation[0:expInfo.nPracTrials] ,operation2) )
 
set_size_write=[]
for i in range(len(letter)):
    set_size_write.append(len(letter[i]))
for i in range(len(letter2)):
    set_size_write.append(len(letter2[i]))


 
all_write=[]
for i in range(15): #15 trials in total
    
    #---create new row_write variable
    row_write=[]
    for j in range(25 ): #change the numer if the numer of line writen is changed
        row_write.append(' ')
    #---Finish row_write vairable
     
    row_write[0]=expInfo_now['participant']
    row_write[1]=block_num_write[i]
    row_write[2]=trial_num_write[i]
    row_write[3]=set_size_write[i]

    #-----write op_ans etc.
    pr=4 #prior=4  change this line if change priors
    for j in range(7): #7 answer at most each row (each trial)
        if j<len( (op_ans_write[i]) ) :
            row_write[j+pr]=op_ans_write[i][j] 
            row_write[j+pr+7]=letter_write[i][j]
            row_write[j+pr+7+7]=op_write[i][j]
        else:
            row_write[j+pr]=' '
            row_write[j+pr+7]=' ' 
            row_write[j+pr+7+7]=' '
    
    all_write.append(row_write)
             
with myFile:
    writer=csv.writer(myFile)
    writer.writerows([['Sub_num','Block_num','Trial_num','set_size',
        'op_ans1','op_ans2','op_ans3','op_ans4','op_ans5','op_ans6','op_ans7',
        'Letter1','Letter2','Letter3','Letter4','Letter5','Letter6','Letter7',
        'op1','op2','op3','op4','op5','op6','op7']])

    writer.writerows(all_write)

Message(expInfo.thankYouMsg)
event.waitKeys(keyList=["space"])
core.quit()

