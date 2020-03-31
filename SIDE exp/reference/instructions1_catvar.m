function instructions1_catvar(wind1)
%yoffset=input('yoffset ');
yoffset=30;
%instructsize=input('instructsize ');
instructsize=20;
%[wind1 rect] = Screen('OpenWindow',0,[100 100 175],[50 50 1700 900]);
% [wind1 rect] = Screen('OpenWindow',0,[100 100 175]);
sent1='On each trial of this experiment ';
sent2=' you will be presented with a study list of 2;4 or 8 pictures.';
sent3='Following the study list, you will be presented with a test picture.';
sent4='If the test picture was a member of the study list,';
sent5=' press either F or J key depending on the side of the presentation. Left is F; right is J';
sent6='If the test picture was not a member of the study list,';
sent7=' press the other key';
sent8='After you make your response, the computer will tell you if';
sent9=' you were correct or incorrect.';
sent10='Please always rest your fingers on the F and J keys';
sent11=' and try to respond as quickly as you can without making errors.';
sent12='The experiment has 9 blocks of trials with 25 trials in each block.';
sentlast='Please press SPACE to begin';
blank=' ';

sentence={sent1 sent2 sent3 sent4 sent5 sent6 sent7 blank sent8 sent9 blank sent10 sent11 blank sent12 sentlast};
Screen('TextSize',wind1,instructsize);
for i=1:15
    Screen('DrawText',wind1,sentence{i},50,1000-(26-i)*yoffset-400)
end
Screen('DrawText',wind1,sentlast,50,1000-10*yoffset);
Screen('Flip',wind1)
%%
%      user presses space for next screen
legal=0;
while legal == 0
    [keydown secs keycode]=KbCheck;
    key=KbName(keycode);
    if strcmp(key,'space')
        legal=1;
    end
end
Screen('Flip',wind1)
%
%
WaitSecs(.5)

