  %%  CMVMmixBlocked_Rapid.m
% VM (8 items/block) & CM (8 items/block) in mixed between blocks 
% 2 blocks, each block contains 100 trials
% filler items: 8 x ntrials x nblocks
% total number of test stimuli: 16 x #of blocks
% mem-set size is either 2 or 4
% Rapid presentation: 500ms after each item
% color condition: 1 (red as target color, green as foil color); 0 (green
% as target color)
% Side condition: 1 (target stimuli at right); -1(target stimuli at left)
% Things to do: change the filler to be a fixed set
% To compensate for lack of CM, a special version with just CM, 50 trials
% per block for 4 blocks
% 1/26/2018
%%
image_location=[pwd '\images\'];
data_location=[pwd '\data\'];
subid=input(' subject # '); 
color_con=input(' color condition(1 or 2): '); 
filename=[data_location 'CMVMmixBlocked_Rapid' num2str(subid) '.txt'];
variedvars=[data_location 'CMVMmixBlocked_Rapid' num2str(subid)];
if ~exist(filename)
    fid=fopen(filename,'wt');
    s=RandStream('mt19937ar','Seed','shuffle');
    RandStream.setGlobalStream(s);
    HideCursor;
    % 
    % allow user to input experiment parameters for pilot testing
    % 
    %nblocks=input('# of blocks for pilot test ');
    %ntrials=input('# of trials for pilot test '); 
    nblocks=4;   % have to be an even number CHANGED FOR SPECIAL CONDITION
    ntrials=50; 
    %recsiz=input('size of color squares (40) ');
    recsiz=250;  %% IMPORTANT defult 200
    %dist_cetral=input('distance to the central point');
    %textsize=input('size of text (30) ');
    textsize=30;
    %fixation_size=input('size of fixation ');
    fixation_size=20;
    %nsets=input('# of mem set sizes');
    nsets=2; %change the set size #
    %present_rate=input('study time duration (secs) ');
    present_rate=.1;
    %isi=input('isi (secs) ');
    isi=.5;  %%%%%%%%%%%%%% CHANGE ##### DON'T GO MORE THAN 0.9%%%We tried 0.5, we can use this isi but it could be difficult to attend for a long time period.
    %pre test probe time for CDA collection
    pre_test_presentation=0.9-isi;
    %retention_interval=input('retention interval ');
    retention_interval=.1;
    %dist=input('distance to the center point ');
    dist=200; %distance from the fixation point
    p = 0.5; %chance that all presentations are at the same side
    minibreak=50; % Allow subject to have a break within a block
    %
    %%
    % set all constants
    %
    set_sizes=[2 4]; 
    images=dir([image_location '*.jpg']);
    for i=1:length(images)
        fullset{i}=images(i).name;
    end
    session_order=randperm(length(fullset)); % SET THE STIMULI ORDER, DONT CHNAGE ACROSS BLOCKS
    for i=1:8*nblocks  %for the filler images
        filler_set{i}=fullset{session_order(i+8*nblocks)};
    end
    %%
    % Condition & color order detemermined before the test %
     order1=randperm(nblocks);
%     order2=randperm(nblocks);
    condition_set=1:nblocks; %we have 2 blocks, 1 CM, >1 is VM
    for i=1:nblocks
        condition_list{i}=condition_set(order1(i));  
    end
%     color_set=[zeros(1, 0.5*nblocks),ones(1,0.5*nblocks)];
%     for i=1:nblocks
%         condition_list{i}=condition_set(order1(i));  
%         color_list{i}=color_set(order2(i));
%     end
    %%    
    %[wind1 rect] = Screen('OpenWindow',0,[175 175 175],[100 100 1500 900]);
    [wind1 rect] = Screen('OpenWindow',0,[175 175 175]);  %THE SCREEN 
    centerx=(rect(3)-rect(1))/2;
    centery=(rect(4)-rect(2))/2;
    fixation='*';  %will be in every frame during study phase
    prompt='+';
    rest_index='Rest Index Fingers on F and J';
    press_space='When Ready, Press Space to Begin Block ';
    red_target_color='For this block, please only pay attention to images with red frame';
    green_target_color='For this block, please only pay attention to images with green frame';
    block_end='End of Block';
    warning='Do Not Press Key Until Probed!';
    thanks='Thank You, the Experiment is Over!';
    text_correct='CORRECT!';
    text_incorrect='INCORRECT';
    percentage='Percentage Correct = ';
    mini_break='Trials completed: ';
    space_cont='Press Space to Continue';
    Screen('TextSize',wind1,textsize)
    textbounds_rest_index=Screen('TextBounds',wind1,rest_index);
    textbounds_color_target_red=Screen('TextBounds',wind1,red_target_color);
    textbounds_color_target_green=Screen('TextBounds',wind1,green_target_color); %Set the color instruction text bound
    textbounds_thanks=Screen('TextBounds',wind1,thanks);
    textbounds_press_space=Screen('TextBounds',wind1,press_space);
    textbounds_block_end=Screen('TextBounds',wind1,block_end);
    textbounds_percentage=Screen('TextBounds',wind1,percentage);
    textbounds_mini_break=Screen('TextBounds',wind1,mini_break);
    textbounds_space_cont=Screen('TextBounds',wind1,space_cont);
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
    trial_type_store=[];
	color_set_store=[];
    test_set_stor=[]; 
    stimu_set_stor=[];
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
    loc=[x1 y1 x2 y2];  %location for test probes
    locL=[x1-dist y1 x2-dist y2]; %location for the left stimuli
    locR=[x1+dist y1 x2+dist y2]; %location for the left stimuli
    WaitSecs(1)
    %%
    %  start of actual experiment
    %
    instructions1_catvar(wind1);
    tot_trials=0;
    istim=0;

    for block=1:nblocks
        percent_correct=0;
        mean_rt=0;
        condition=condition_list{block};
        % get the stimuli for current block
        testCon=2;  %testCon 1: CM, 2 VM.
        for i=1:8
            block_stimuli{i}=fullset{session_order((block-1)*8+i)};
        end
        testCon=2;  %testCon 1: CM, 2 VM.
        % IF VM set in current block
        for i=1:8
            VM_simuli_set{i}=block_stimuli{i};
        end
        if condition<5   %CHANGE
         % If CM set in current block
            testCon=1;
            for i=1:4
                positive_set{i}=block_stimuli{i};
                negative_set{i}=block_stimuli{i+4};
            end
        end

        block_filler_order=(1+(block-1)*8):block*8; %%% set of position to draw from filler set
        %%
        block_start=[press_space num2str(block)];
  %      color_con=color_list{block}; Manually set color blocks
        if color_con == 1
            color_instr=red_target_color;
            textbounds_color_instr=textbounds_color_target_red;
            color_blockT=[160 0 0];  %when the target color is red
            color_blockF=[0 160 0];  %and the filler color  is green
        else
            color_instr=green_target_color;
            textbounds_color_instr=textbounds_color_target_red;
            color_blockT=[0 160 0];  %target color is green
            color_blockF=[160 0 0];
        end
        Screen('TextSize',wind1,textsize)
        Screen('DrawText',wind1,rest_index,rect(3)/2-textbounds_rest_index(3)/2,rect(4)/2-textbounds_rest_index(4)/2-200)
        Screen('DrawText',wind1,color_instr,rect(3)/2-textbounds_color_instr(3)/2,rect(4)/2-textbounds_color_instr(4)/2+200)
        Screen('DrawText',wind1,block_start,rect(3)/2-textbounds_press_space(3)/2,rect(4)/2-textbounds_press_space(4)/2+400)
        Screen('FrameRect',wind1,[color_blockT(1) color_blockT(2) color_blockT(3)],[loc(1) loc(2) loc(3) loc(4)],7);
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
        %%
        %      present ntrials memory-search trials
        %
        for trial=1:ntrials
            index=randi(nsets);
            setsize=set_sizes(index);
            trial_filler_order=randperm(8); % position to draw from the block_filler_order
            %%
            % Check if subject might need a break
            %
            if trial~=1 && rem(trial-1,minibreak)==0
                minibreaktrials=[mini_break num2str(trial-1)];
                Screen('TextSize',wind1,textsize);
                Screen('DrawText',wind1,minibreaktrials,rect(3)/2-textbounds_mini_break(3)/2,rect(4)/2-textbounds_mini_break(4)/2-200);
                Screen('DrawText',wind1,space_cont,rect(3)/2-textbounds_space_cont(3)/2,rect(4)/2-textbounds_space_cont(4)/2)
                Screen('Flip',wind1) 
                legal=0;
                while legal == 0
                    [keydown secs keycode]=KbCheck;
                    key=KbName(keycode);
                    if strcmp(key,'space')
                        legal=1;
                    end
                end
            end

            %%
            % 
            % Old vs. New trial            
            old=1;
            if rand<.5
                old=2;
                serpos=0;
                lag=0;
            end
            if old==1
                serpos=randi(setsize);
                lag=setsize-serpos+1;
            end
            % CM vs. VM trial            
            
            % Start the trial      
            if testCon==2
                trial_order=randperm(8);  % randomize the 8 stimuli for VM
                for k=1:setsize
                    memset{k}=block_stimuli{trial_order(k)};%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%CHANGE BACK
                    trial_filler{k}=filler_set{block_filler_order(trial_filler_order(k))};
                    if rand < .500001
                        sideset(k)=1; % side set determins the stimuli position: 1 (right); -1 (left)
                    else
                        sideset(k)=-1;  
                    end
                end
                if old==1
                    teststim=memset{serpos};
                    kstim=trial_order(serpos);
                else
                    teststim=block_stimuli{trial_order(setsize+1)};
                    kstim=trial_order(setsize+1);
                end
            else  % for CM condition
                trial_order=randperm(4); %randomize the positive set
                for k=1:setsize
                    memset{k}=positive_set{trial_order(k)};
                    trial_filler{k}=filler_set{block_filler_order(trial_filler_order(k))};
                    if rand < .500001
                        sideset(k)=1; % side set determins the stimuli position: 1 (right); -1 (left)
                    else
                        sideset(k)=-1; %change to test the filler  
                    end
                end
                if old==1
                    kstim=trial_order(serpos);
                    teststim=memset{serpos};
                else
                    k=randi(4);
                    kstim=k+4;
                    teststim=negative_set{k};
                end
            end
            % Reset the presentation side regardless if it is CM or VM
            if rand<p
            for k =2:setsize
                sideset(k)=sideset(k-1);
            end
            end
            %%
            %   present the fixation point
            %
            Screen('TextSize',wind1,fixation_size)
            Screen('DrawText',wind1,fixation,rect(3)/2-textbounds_fixation(3)/2,rect(4)/2-textbounds_fixation(4)/2)
            Screen('Flip',wind1)
            WaitSecs(.5)
            Screen('Flip',wind1)
            WaitSecs(.1)
            %%
            %          present the memory set
            %
            for i=1:setsize
                img_stimuli=imread([image_location memset{i}]);      % call the stimuli image  
                img_filler=imread([image_location trial_filler{i}]); % CHECK call the filler image
                picture_stimuli=Screen('MakeTexture',wind1,img_stimuli);
                picture_filler=Screen('MakeTexture',wind1,img_filler);
                if sideset(i)==1
                    locS=locR;  % when the stimuli is on the right side of fixation 
                    locF=locL;   % the filler is on the left side of fixation
                else
                    locS=locL;  
                    locF=locR;   
                end
                % DRAW the color frame with stimuli color and filler
                % color 
                Screen('DrawTexture',wind1,picture_stimuli,[],[locS(1) locS(2) locS(3) locS(4)]);  % Location for the stimuli
                Screen('DrawTexture',wind1,picture_filler,[],[locF(1) locF(2) locF(3) locF(4)]);  % Location of the filler
                Screen('FrameRect',wind1,[color_blockT(1) color_blockT(2) color_blockT(3)],[locS(1) locS(2) locS(3) locS(4)],7); % Draw a frame with target color for the stimuli
                Screen('FrameRect',wind1,[color_blockF(1) color_blockF(2) color_blockF(3)],[locF(1) locF(2) locF(3) locF(4)],7); % Draw a frame with filler color for the filler
                Screen('TextSize',wind1,fixation_size) %keep the fixation point
                Screen('DrawText',wind1,fixation,rect(3)/2-textbounds_fixation(3)/2,rect(4)/2-textbounds_fixation(4)/2)
                Screen('Close',picture_stimuli);
                Screen('Close',picture_filler);
                Screen('Flip',wind1);
                WaitSecs(present_rate);
                %Screen('Flip',wind1)
                Screen('TextSize',wind1,fixation_size) %keep the fixation point
                Screen('DrawText',wind1,fixation,rect(3)/2-textbounds_fixation(3)/2,rect(4)/2-textbounds_fixation(4)/2)
                Screen('Flip',wind1)
                WaitSecs(isi);
                Screen('Flip',wind1)
            end
            
            % For the last item presented, we need to ensure a full 0.9 sec
            % to collect CDA information
            %
            %  present the fixation point for CDA
            %
            Screen('TextSize',wind1,fixation_size); %keep the fixation point
            Screen('DrawText',wind1,fixation,rect(3)/2-textbounds_fixation(3)/2,rect(4)/2-textbounds_fixation(4)/2);
            WaitSecs(pre_test_presentation);
            Screen('Flip',wind1)
            %
            %   present the prompt
            %
            Screen('TextSize',wind1,fixation_size)
            Screen('DrawText',wind1,prompt,rect(3)/2-textbounds_fixation(3)/2,rect(4)/2-textbounds_fixation(4)/2)
            Screen('Flip',wind1)
            WaitSecs(.5)
            Screen('Flip',wind1)
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
                img_data=imread([image_location teststim]);
                picture=Screen('MakeTexture',wind1,img_data);
                Screen('DrawTexture',wind1,picture,[],[loc(1) loc(2) loc(3) loc(4)]);
                Screen('FrameRect',wind1,[color_blockT(1) color_blockT(2) color_blockT(3)],[loc(1) loc(2) loc(3) loc(4)],7); % Draw a frame with target color for the stimuli
                Screen('Flip',wind1);
                Screen('Close',picture);
                start=GetSecs;
                legal=0;
                while legal == 0
                    [keydown secs keycode] = KbCheck;
                    key=KbName(keycode);
                    if any(strcmp(key,legalkeys))
                        rt=secs-start;
                        legal=1;
                    end
                end
                mean_rt=mean_rt+rt;
                Screen('Flip',wind1)
                WaitSecs(.5)
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
                %%
                %             present feedback
                %
                Screen('TextSize',wind1,textsize)
                if corr == 1
                    percent_correct=percent_correct+1;
                    Screen('DrawText',wind1,text_correct,rect(3)/2-textbounds_correct(3)/2,rect(4)/2-textbounds_correct(4)/2);
                    Screen('Flip',wind1)
                else
                    Screen('DrawText',wind1,text_incorrect,rect(3)/2-textbounds_incorrect(3)/2,rect(4)/2-textbounds_incorrect(4)/2);
                    Screen('Flip',wind1)
                end
                WaitSecs(1);
                Screen('Flip',wind1)
                %%
                %               record results
                %
                tot_trials=tot_trials+1;
                block_store(tot_trials)=block;
                trial_store(tot_trials)=trial;
                serpos_store(tot_trials)=serpos;
                lag_store(tot_trials)=lag;
                trial_type_store(tot_trials)=old;
                color_set_store(tot_trials)=color_con;
                test_set_stor(tot_trials)=testCon;
                %stimu_set_stor(tot_trials)=stimu_type;
                setsize_store(tot_trials)=setsize;
                resp_store(tot_trials)=resp;
                rt_store(tot_trials)=rt;
                corr_store(tot_trials)=corr;
                probe_store(tot_trials)=kstim;
                for i=1:setsize
                    set_store(tot_trials,i)=trial_order(i)*sideset(i);
                end
                %%
                %               write to output text file
                %
                fprintf(fid,'%5d',block_store(tot_trials),trial_store(tot_trials),trial_type_store(tot_trials), color_set_store(tot_trials),...
                     test_set_stor(tot_trials),setsize_store(tot_trials),serpos_store(tot_trials),...
                    lag_store(tot_trials),resp_store(tot_trials),corr_store(tot_trials));
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
                Screen('DrawText',wind1,warning,rect(3)/2-textbounds_warning(3)/2,rect(4)/2-textbounds_warning(4)/2);
                Screen('Flip',wind1)
                WaitSecs(5)
            end
        end
        percent_correct=round(100*percent_correct/ntrials);
        percentage_correct=[percentage num2str(percent_correct)];
        mean_rt=mean_rt/ntrials;
        Screen('TextSize',wind1,textsize)
        Screen('DrawText',wind1,percentage_correct,rect(3)/2-textbounds_percentage(3)/2,rect(4)/2-textbounds_percentage(4)/2-200);
        Screen('DrawText',wind1,block_end,rect(3)/2-textbounds_block_end(3)/2,rect(4)/2-textbounds_block_end(4)/2)
        Screen('Flip',wind1)
        WaitSecs(4)
    end
    fclose(fid);
    save(variedvars);
    debriefing(wind1);
    Screen('DrawText',wind1,thanks,rect(3)/2-textbounds_thanks(3)/2,rect(4)/2-textbounds_thanks(4)/2);
    Screen('Flip',wind1);
    WaitSecs(.5);
    legal=0;
    while legal==0
        [keydown,secs,keycode]=KbCheck;
        if strcmp(KbName(keycode),'q')
            legal=1;
        end
    end
    clear screen;
    clear all;
else
    disp('Error: the file already exists!')
end
