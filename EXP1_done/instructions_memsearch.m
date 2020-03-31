function instructions_memsearch(wind1,rect)
%yoffset=input('yoffset ');
yoffset=30;
%instructsize=input('instructsize ');
instructsize=20;
centerx=(rect(3)-rect(1))/2;
centery=(rect(4)-rect(2))/2;
xadjust=70;
%[wind1 rect] = Screen('OpenWindow',0,[100 100 175],[50 50 1700 900]);
% [wind1 rect] = Screen('OpenWindow',0,[100 100 175]);
sent1='This experiment tests your ability to remember pictures of objects. ';
sent2='On each trial, a list of pictures of objects will be presented.';
sent3='Next, a test picture will be presented.  (The computer will ';
sent3b=' present a plus sign right before it presents the test picture.)';
sent4='If the test picture was on the list, then press the OLD key (J) ';
sent5=' as quickly as possible.';
sent6='If the test picture was not on the list, then press the NEW key (F)';
sent7=' as quickly as possible. ';

blank=' ';

sent8='Please note that a picture is defined as OLD only if it was on ';
sent9=' the CURRENT list. Pictures that appeared on PREVIOUS lists but ';
sent10=' not on the current one are defined as NEW.';

sent11='We will be measuring both your accuracy ';
sent12=' and how long it takes you to make your responses.';
sent13='Therefore, please always rest your right index finger on the OLD key';
sent14=' and your left index finger on the NEW key.'; 
sent15='Try to make your responses as quickly as you can without making errors.';

sent16='The experiment has 9 blocks of 25 test lists each.';
sent17='Please ask the experimenter now if you have any questions.';
sentlast=' PRESS SPACE TO CONTINUE';
blank=' ';
sentence={sent1 sent2 sent3 sent3b sent4 sent5 sent6 sent7 blank sent8 sent9 sent10 blank sent11 sent12 sent13 sent14 sent15 blank sent16 sent17 sentlast};
Screen('TextSize',wind1,instructsize);
textbounds_sentlast=Screen('Textbounds',wind1,sentlast);

for i=1:21
    Screen('DrawText',wind1,sentence{i},50,1000-(21-i)*yoffset-400);
end
Screen('DrawText',wind1,sentlast,rect(3)/2-textbounds_sentlast(3)/2,rect(4)-50);
Screen('Flip',wind1);
%%
%      user presses space when ready to start
%
legal=0;
while legal == 0
    [keydown secs keycode]=KbCheck;
    key=KbName(keycode);
    if strcmp(key,'space')
        legal=1;
    end
end
Screen('Flip',wind1);
WaitSecs(.5);




