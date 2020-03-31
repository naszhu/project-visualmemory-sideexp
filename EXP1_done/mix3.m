%%mixed cmvman experiment
% mixedcmvman.m
% replicate rui's experiment
% 3 CM+, 3CM-, 6VM, unlimited AN
% 01/02/19
%
clear all;
image_location=[pwd '\images\'];
data_location=[pwd '\data2\'];
subid=input(' subject # ');
filename=[data_location 'cmvman' num2str(subid) '.txt'];
variedvars=[data_location 'cmvman' num2str(subid)];

%% debug or test on my laptop
debug=0;
mylaptop=1;


%% File starts
if ~exist(filename)
    fid=fopen(filename,'wt');
    s=RandStream('mt19937ar','Seed','shuffle');
    RandStream.setGlobalStream(s);
    
    if debug==0
        HideCursor;
    end
    %
    % allow user to input experiment parameters for pilot testing
    %corr
    %n=input(' n ');
    n=6; % set size of memory pool
    %nblocks=input('# of blocks for pilot test ');
    %ntrials=input('# of trials for pilot test ');
    nblocks=7;
    ntrials=25;
    %recsiz=input('size of color squares (40) ');
    recsiz=200;
    %textsize=input('size of text (30) ');
    textsize=30;
    %fixation_size=input('size of fixation ');
    fixation_size=20;
    %nsets=input('# of mem set sizes');
    nsets=3;
    %present_rate=input('study time duration (secs) ');
    if debug==1
        present_rate=0.1;
        retention_interval=0;
        isi=0;
    else
        present_rate=1;
        retention_interval=.1;
        isi=.1;
    end
    %isi=input('isi (secs) ');

    %retention_interval=input('retention interval');

    %
    %% set all constants
    %
    set_sizes=[2 4 6];
    images=dir([image_location '*.jpg']);
    nstim=length(images);
    for i=1:nstim
        fullset{i}=images(i).name;  %contain full set of images name
    end
    order=randperm(nstim);  %permutation stimuli
    for i=1:nstim
        pictures{i}=fullset{order(i)}; %contain cell of permuted pictures
    end
    %
    % Establish the stimulus sets for the conditions 
    %
    %% participant wide
    nvm=6; %number of vms
    cmpos_s=1;
    cmpos_r=2;
    cmneg_h=3;
    cmneg_l=4;
    for i=1:nvm
        vm(i)=4+i; %5,6,7,8,9,10
    end

    %% Define some variables
    wss={[0.12346,0.4,0.47654],[0.12346,0.8,0.07654],[0.12346,0.87654,0]};
    %prob(test|test old)
    probchart_old={[0.9,0.1],[0.9,repelem(1/30,3)],[0.9,repelem(1/50,5)];...
    [1/2,1/2], [repelem(1/4,4)],[0.228169,repelem(0.154366,5)];...
    [1/2,1/2],[repelem(1/4,4)],[0]};
    %prob(test i |test new)
    prob_c_d=[1/3,1/6];
    probchart_new={[prob_c_d,repelem(1/10,5)],[prob_c_d,repelem(1/6,3)],[prob_c_d,1/2];...
    [prob_c_d,repelem(1/10,5)],[prob_c_d,repelem(1/6,3)],[prob_c_d,1/2];...
    [prob_c_d,repelem(1/8,4)],[1/3,1/6,repelem(1/4,2)],[0]};
    %
    %   all-new stimuli will have code numbers greater than 2n
    %
    if mylaptop==1
        Screen('Preference', 'SkipSyncTests', 1);
    end
    if debug==1
        [wind1 rect] = Screen('OpenWindow',0,[175 175 175],[100 100 900 400]);
    else %if i'm testing
        [wind1 rect] = Screen('OpenWindow',0,[175 175 175]);
    end
    

    centerx=(rect(3)-rect(1))/2;
    centery=(rect(4)-rect(2))/2;
    fixation='*';
    prompt='+';
    rest_index='Rest Index Fingers on F and J';
    press_space='When Ready, Press Space to Begin Block ';
    block_end='End of Block';
    warning='Do Not Press Key Until Probed!';
    thanks='Thank You, the Experiment is Over!';
    text_correct='CORRECT!';
    text_incorrect='INCORRECT';
    percentage='Percentage Correct = ';
    Screen('TextSize',wind1,textsize);
    textbounds_rest_index=Screen('TextBounds',wind1,rest_index);
    textbounds_thanks=Screen('TextBounds',wind1,thanks);
    textbounds_press_space=Screen('TextBounds',wind1,press_space);
    textbounds_block_end=Screen('TextBounds',wind1,block_end);
    textbounds_percentage=Screen('TextBounds',wind1,percentage);
    textbounds_warning=Screen('TextBounds',wind1,warning);
    textbounds_correct=Screen('TextBounds',wind1,text_correct);
    textbounds_incorrect=Screen('TextBounds',wind1,text_incorrect);
    Screen('TextSize',wind1,fixation_size);
    textbounds_fixation=Screen('TextBounds',wind1,fixation);
    %iok=input('textbounds okay?')
    block_store=[];
    trial_store=[];
    rt_store=[];
    resp_store=[];
    correct_store=[];
    cond_store=[];
    old_new_store=[];
    setsize_store=[];
    serpos_store=[];
    lag_store=[];
    probe_store=[];
    legalkeys={'f','j'};
    %%
    %  compute location coordinates
    %   these coordinates are determined by the size of the images
    %
    x1=centerx-recsiz/2;
    y1=centery-recsiz/2;
    x2=x1+recsiz;
    y2=y1+recsiz;
    loc=[x1 y1 x2 y2];
    if debug~=1
        WaitSecs(1);
    end
    %%
    instructions_memsearch(wind1,rect);
    %  start of actual experiment
    %
    tot_trials=0;
    istim=0;
    antot=0; %all new total
    for block=1:nblocks
        percent_correct=0;
        mean_rt=0;
        block_start=[press_space num2str(block)];
        Screen('TextSize',wind1,textsize);
        Screen('DrawText',wind1,rest_index,rect(3)/2-textbounds_rest_index(3)/2,rect(4)/2-textbounds_rest_index(4)/2-200);
        Screen('DrawText',wind1,block_start,rect(3)/2-textbounds_press_space(3)/2,rect(4)/2-textbounds_press_space(4)/2);
        Screen('Flip',wind1);
        %%
        %      user presses space when ready to start
        %
        legal=0;
        if debug==0
            while legal == 0
                [keydown secs keycode]=KbCheck;
                key=KbName(keycode);
                if strcmp(key,'space')
                    legal=1;
                end
            end
        else
            legal=1;
        end
        %%
        %      present ntrials memory-search trials
        %
        for trial=1:ntrials

            index=randi(nsets); % set sizes was defined previously
            setsize=set_sizes(index); %random select a set size
            memlist=[]; %start from fresh, though not necessary because for loop clear it
            foilist=[];
            ss=setsize;  % abbrev. for convenience; ss is the current set size

            %% Memory Set (MS):
            % give what sequence we want and what stimuli we want
            % memset is a number's set 

            %% foil list:
            %the first two are always c,d
            foilist(1)=cmneg_h; %=3
            foilist(2)=cmneg_l; %=4
 

            %% Memory list           
            testkind=["A","B","nAB"];

            whatcm=randsample(testkind,1,true,wss{ss/2}); %ss/2 =1,2,3
            cmassign=0; %a variable to indicate what CM is included

            if whatcm=="A"
                memlist(1)=cmpos_s;
                ordvm=randperm(nvm);
                for i=1:(ss-1) %1:1,1:2,1:5
                    memlist(i+1)=vm(ordvm(i)); %i+1 because 1 with A
                end
                % foil list
                for i=1:(nvm-ss+1) %1:5,1:3,1:1
                    foilist(i+2)=vm(ordvm(i+(ss-1))); 
                    %ss-1 is number of vm in memlist, the followings are what left in foilset
                end
                cmassign=1;
            elseif whatcm=="B"
                memlist(1)=cmpos_r;
                ordvm=randperm(nvm);
                for i=1:(ss-1)
                    memlist(i+1)=vm(ordvm(i));
                end  
                %foil list
                for i=1:(nvm-ss+1) %1:5,1:3,1:1
                    foilist(i+2)=vm(ordvm(i+(ss-1))); 
                    %ss-1 is number of vm in memlist, the followings are what left in foilset
                end
                cmassign=2;
            elseif whatcm=="nAB" %ss=6 has 0 probabilty to reach here 
                ordvm=randperm(nvm);
                for i=1:ss %ex. 1:2,1:4,1:6
                    memlist(i)=vm(ordvm(i));
                end
                cmassign=3;
                % foil list
                for i=1:(nvm-ss) %1:5,1:3,1:1
                    foilist(i+2)=vm(ordvm(i+(ss))); 
                    %ss-1 is number of vm in memlist, the followings are what left in foilset
                end
            else 
                disp("Error ss manipulation");
            end


            %%  Randomize memset
            ordss=randperm(ss);
            memset=[];
            for i=1:ss
                memset(i)=memlist(ordss(i)); % permutation memset
            end

            %% Test item:
            % randomly choose old/new 
            %

            %% Decide old or new
            old=1;
            whattest=rand;
            if memlist(1)==cmpos_s % with special cm+
                if whattest<0.1 %test new
                    old=2;
                end
            else % not with special CM+,then all 50%
                if whattest<0.5 %test new
                    old=2;
                end
            end

            %% Draw test item
            if old==1
                disp(probchart_old{cmassign,ss/2});
                disp(memlist);
                disp(foilist);
                testitem = randsample(memlist,1,true,probchart_old{cmassign,ss/2}); %there is no 3,6 
                serpos = find(memset==testitem);
                lag=setsize-serpos+1;
            else % New
                testitem = randsample(foilist,1,true,probchart_new{cmassign,ss/2});
                serpos=0;
                lag=0;
            end


            %% Determin what condition is this
            cond=0; %clear condition
            if testitem==cmpos_s
                cond=1;
            elseif testitem==cmpos_r
                cond=2;
            elseif testitem==cmneg_h
                cond=3;
            elseif testitem==cmneg_l
                cond=4;
            else %vm
                cond=5;
            end

            %%
                        %   record trial information
            %
                tot_trials=tot_trials+1;
                block_store(tot_trials)=block;
                trial_store(tot_trials)=trial;
                serpos_store(tot_trials)=serpos;
                lag_store(tot_trials)=lag;
                old_new_store(tot_trials)=old;
                cond_store(tot_trials)=cond;
                setsize_store(tot_trials)=setsize;
                probe_store(tot_trials)=testitem;
                for i=1:setsize
                    set_store(tot_trials,i)=memset(i);
                end
            
            %
            %
            %   present the fixation point
            %
            Screen('TextSize',wind1,fixation_size);
            Screen('DrawText',wind1,fixation,rect(3)/2-textbounds_fixation(3)/2,rect(4)/2-textbounds_fixation(4)/2);
            Screen('Flip',wind1);
            if debug~=1
                WaitSecs(.5);
            end
            Screen('Flip',wind1);
            if debug~=1
                WaitSecs(.1);
            end
            %%
            %          present the memory set
            %
            for i=1:setsize
                img_data=imread([image_location pictures{memset(i)}]);
                picture=Screen('MakeTexture',wind1,img_data);
                Screen('DrawTexture',wind1,picture,[],[loc(1) loc(2) loc(3) loc(4)]);
                Screen('Close',picture);
                Screen('Flip',wind1);
                if debug~=1
                    WaitSecs(present_rate);
                end
                Screen('Flip',wind1);
                WaitSecs(isi);
            end
            %
            %   present the prompt
            %
            Screen('TextSize',wind1,fixation_size);
            Screen('DrawText',wind1,prompt,rect(3)/2-textbounds_fixation(3)/2,rect(4)/2-textbounds_fixation(4)/2);
            Screen('Flip',wind1);
            if debug~=1
                WaitSecs(.5);
            end
            Screen('Flip',wind1);
            %%
            %          retention interval; make sure
            %             that subject does not press key
            %             prematurely
            %
            flag=0;
            start=GetSecs;
            while GetSecs-start<retention_interval
                [keydown secs keycode]=KbCheck;
                if keydown==1
                    flag=1;
                end
            end
            %%
            %          present probe and collect legal response
            %
            if flag == 0
                img_data=imread([image_location pictures{testitem}]);
                picture=Screen('MakeTexture',wind1,img_data);
                Screen('DrawTexture',wind1,picture,[],[loc(1) loc(2) loc(3) loc(4)]);
                Screen('Flip',wind1);
                Screen('Close',picture);
                start=GetSecs;
                legal=0;
                if debug==0
                    while legal == 0
                        [keydown secs keycode] = KbCheck;
                        key=KbName(keycode);
                        if ischar(key)
                            if any(strcmp(key,legalkeys))
                                rt=secs-start;
                                legal=1;
                            end
                        end
                    end
                else
                    key='f';
                    rt=0;
                    legal=1;
                end
                
                
                mean_rt=mean_rt+rt;
                Screen('Flip',wind1);
                if debug~=1
                    WaitSecs(.5);
                end
                %%
                %             determine the response and whether it is correct
                %
                if strcmp(key,'j')
                    resp=1;
                else
                    resp=2;
                end
                corr=0;
                if old == 1 && resp == 1
                    corr=1;
                end
                if old == 2 && resp == 2
                    corr=1;
                end
                
                %
                % record response information
                %
                resp_store(tot_trials)=resp;
                rt_store(tot_trials)=rt;
                corr_store(tot_trials)=corr;
                
                %%
                %             present feedback
                %
                Screen('TextSize',wind1,textsize);
                if corr == 1
                    percent_correct=percent_correct+1;
                    Screen('DrawText',wind1,text_correct,rect(3)/2-textbounds_correct(3)/2,rect(4)/2-textbounds_correct(4)/2);
                    Screen('Flip',wind1);
                else
                    Screen('DrawText',wind1,text_incorrect,rect(3)/2-textbounds_incorrect(3)/2,rect(4)/2-textbounds_incorrect(4)/2);
                    Screen('Flip',wind1);
                end
                if debug==0
                    WaitSecs(1);
                else
                    WaitSecs(0.1);
                end
                Screen('Flip',wind1);
                %%
                %               record results
                %
                
                %%
                %               write to output text file
                %
                fprintf(fid,'%5d',block_store(tot_trials),trial_store(tot_trials),cond_store(tot_trials),...
                    old_new_store(tot_trials),setsize_store(tot_trials),serpos_store(tot_trials),lag_store(tot_trials),...
                    resp_store(tot_trials),corr_store(tot_trials));
                fprintf(fid,'%10d',round(1000*rt_store(tot_trials)));
                fprintf(fid,'%5d',probe_store(tot_trials));
                for i=1:setsize
                    fprintf(fid,'%5d',set_store(tot_trials,i));
                end
                fprintf(fid,'\n');
            else
                %%
                %               warn subject for premature key press
                %
                resp_store(tot_trials)=9;
                rt_store(tot_trials)=0;
                corr_store(tot_trials)=0;
                 fprintf(fid,'%5d',block_store(tot_trials),trial_store(tot_trials),cond_store(tot_trials),...
                    setsize_store(tot_trials),old_new_store(tot_trials),serpos_store(tot_trials),lag_store(tot_trials),...
                    resp_store(tot_trials),corr_store(tot_trials));
                fprintf(fid,'%10d',round(1000*rt_store(tot_trials)));
                fprintf(fid,'%5d',probe_store(tot_trials));
                for i=1:setsize
                    fprintf(fid,'%5d',set_store(tot_trials,i));
                end
                fprintf(fid,'\n');
                Screen('DrawText',wind1,warning,rect(3)/2-textbounds_warning(3)/2,rect(4)/2-textbounds_warning(4)/2);
                Screen('Flip',wind1);
                if debug~=1
                    WaitSecs(5);
                end
            end
        end

        
        percent_correct=round(100*percent_correct/ntrials);
        percentage_correct=[percentage num2str(percent_correct)];
        mean_rt=mean_rt/ntrials;
        Screen('TextSize',wind1,textsize);
        Screen('DrawText',wind1,percentage_correct,rect(3)/2-textbounds_percentage(3)/2,rect(4)/2-textbounds_percentage(4)/2-200);
        Screen('DrawText',wind1,block_end,rect(3)/2-textbounds_block_end(3)/2,rect(4)/2-textbounds_block_end(4)/2);
        Screen('Flip',wind1);
        if debug~=1 
            WaitSecs(4);
        end
    end
    fclose(fid);
    save(variedvars);
    Screen('DrawText',wind1,thanks,rect(3)/2-textbounds_thanks(3)/2,rect(4)/2-textbounds_thanks(4)/2);
    Screen('Flip',wind1);
    if debug~=1 
        WaitSecs(.5);
    end
    legal=0;
    while legal==0
        [keydown,secs,keycode]=KbCheck;
        if ischar(key)
            if strcmp(KbName(keycode),'q')
                legal=1;
            end
        end
    end
    clear screen;
    clear all;
%else
 %   disp('Error: the file already exists!')
end





