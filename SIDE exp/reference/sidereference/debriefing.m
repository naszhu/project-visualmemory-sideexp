function debriefing(wind1)
%yoffset=input('yoffset ');
yoffset=30;
%instructsize=input('instructsize ');
instructsize=20;
%[wind1 rect] = Screen('OpenWindow',0,[100 100 175],[50 50 1700 900]);
% [wind1 rect] = Screen('OpenWindow',0,[100 100 175]);
sent1='DEBRIEFING:';
blank=' ';
sent2='Depending on the condition in which you were in,';
sent3=' the lists you were trying to remember either:';
sent4=' 1) stayed the same throughout the experiment,';
sent5=' 2) kept switching back and forth, or';
sent6=' 3) varied randomly from trial to trial.';
sent7='The purpose of the experiment is to test how';
sent8=' long-term memories for items affect your ability';
sent9='  to remember them in the short term.';
sent10='Your data will help decide between competing theories';
sent11=' of how long-term memory and short-term memory interact.';
sent12='Thank you for your participation!';
sentence={sent1 blank sent2 sent3 sent4 sent5 sent6 blank sent7 sent8 sent9 sent10 sent11 blank blank sent12};
Screen('TextSize',wind1,instructsize);
for i=1:16
    Screen('DrawText',wind1,sentence{i},50,1000-(21-i)*yoffset-400)
end
%%Screen('DrawText',wind1,sentence{21},50,1000-10*yoffset)
Screen('Flip',wind1)
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
Screen('Flip',wind1)
WaitSecs(.5);

