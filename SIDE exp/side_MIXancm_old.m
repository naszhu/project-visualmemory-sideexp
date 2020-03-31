
%% ok. Lea modified Rui's at - 9/10/2019
debug=1;
debug_block=1;
synct=0;
fn=' MIX ';

image_location=[pwd '\images\']; 
data_location=[pwd '/data/'];
subid=input(' subject # ');

if debug==1
    filename=[data_location 'MixDebug' num2str(subid) '.txt'];
else
    filename=[data_location 'MIX' num2str(subid) '.txt'];
end
consistvars=[data_location 'Mix' num2str(subid)];


if ~exist(filename)
    fid=fopen(filename,'wt');
    s=RandStream('mt19937ar','Seed','shuffle');
    RandStream.setGlobalStream(s);
    if debug==0
        HideCursor;
    end
    if debug_block==0
        nblocks=9; %changed # of blocks
        ntrials=24;
    else 
        nblocks=1; %changed # of blocks
        ntrials=10;
    end
    %recsiz=input('size of color squares (40) ');
    recsiz=200;
    %textsize=input('size of text (30) ');
    textsize=30;
    %fixation_size=input('size of fixation ');
    fixation_size=20;
    %nsets=input('# of mem set sizes');
    nsets=3; %changed # of mem set size
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
    %
    %%
    % set all constants
    %
    %%%%%%%%%%%load all images into the variable block_stimuli
    image_location=[pwd '/images/'];
    images=dir([image_location '*.jpg']);
    fullset=[];%%%%have all names for the images
     for i=1:length(images)
         fullset{i}=images(i).name;
     end
 
 %randomized the order of the images for all blocks
     block_stimuli=[];
    leftside_stimuli=[];
    rightside_stimuli=[];
    %session_order=randperm(length(fullset)); %%%if I want to retrive the name for the image go fullset{[session_order(i)]}


    set_sizes=[2 4 8]; %delete set size 1 and 16
    images=dir([image_location '*.jpg']);

    for i=1:length(images)
        fullset0{i}=images(i).name;
    end
    order=randperm(length(fullset));
    for i=1:length(images)
        fullset{i}=fullset0{order(i)};
    end

    for i=1:16 %8 on each set
        pictures{i}=fullset{i};
    end
    for i=1:8
        CM_left{i}=pictures{i};
        CM_right{i}=pictures{i+8};
    end
    
    if synct==0
        Screen('Preference', 'SkipSyncTests', 0);
    else
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
    textbounds_percentage=Screen('TextBounds',wind1,percentage);
    textbounds_block_end=Screen('TextBounds',wind1,block_end);
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
	study_set_store=[];
    test_set_stor=[];    
    setsize_store=[];
    serpos_store=[];
    lag_store=[];
    probe_store=[];
    legalkeys={'f','j'};
    %%
    %  compute location coordinates
    %   these coordinates are determined by the size of the images
    %
    imageoffset = recsiz+50;
    windowwidth = rect(3);
    windowheight = rect(4);
    horizontalcenter = round(windowwidth/2);
    verticalcenter = round(windowheight/2);
    leftimagecoords = [horizontalcenter-recsiz/2-imageoffset,verticalcenter-recsiz/2,...
        horizontalcenter+recsiz/2-imageoffset,verticalcenter+recsiz/2];
    rightimagecoords = [horizontalcenter-recsiz/2+imageoffset,verticalcenter-recsiz/2,...
        horizontalcenter+recsiz/2+imageoffset,verticalcenter+recsiz/2];
       
    x1=centerx-recsiz/2;
    y1=centery-recsiz/2;
    x2=x1+recsiz;
    y2=y1+recsiz;
    loc=[x1 y1 x2 y2];
    if debug~=1
        WaitSecs(1);
    end
    %%
    %  start of actual experiment
    %
    instructions1_catvar(wind1);
    tot_trials=0;
    antot=0; %total number of an
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
            key='f';
        end
        %%
        %      present ntrials memory-search trials
        %
        

        is.cm=999; %clear up the assignment to see if there is an error 
        for trial=1:ntrials
            index=randi(nsets);
            setsize=set_sizes(index);
            
                        
            testCon =1;

            %% choose old or new
            old=1;%sameside
            if rand<.5
                old=2;
                serpos=0; 
                lag=0;
            end
            if old==1
                serpos=randi(setsize); %randomly choose a number from 1:setsize
                lag=setsize-serpos+1;
            end

            %% CM and AN mixed within trials. 
            memsetlefttemp={}; %store memset before shuffle
            memsetrighttemp={};
            memsetleft={};%empty memory set
            memsetright={}; 
            %memset={};

            Sidechoice=round(rand);
            %%choose memory set half AN, half CM
            ordertemp=randperm(8); %8 CM in the pool at most, so use 8 is good
            ms_stimnum=99999999999;
            for k=1:(setsize/2) %you know, half CM, half AN for no matter setsize

                if Sidechoice==0 %left, only add an ingrediences when test on that side
                    memsetlefttemp{k}=CM_left{ordertemp(k)}; %MEMleft_cm; previously: CM_left{i}=pictures{i};
                    antot=antot+1;
                    memsetlefttemp{k+setsize/2}=fullset{16+antot};

                    ms_stimnum(k)=ordertemp(k);
                    ms_stimnum(k+setsize/2)=16+antot;
                else %right
                    memsetrighttemp{k}=CM_right{ordertemp(k)}; %memright_cm
                    antot=antot+1;
                    memsetrighttemp{k+setsize/2}=fullset{16+antot};
                    
                    ms_stimnum(k)=ordertemp(k);
                    ms_stimnum(k+setsize/2)=16+antot;
                end
            end

            %%randomize the temp memory set

            %ordernow stores the index for
            ordernow=randperm(setsize); % memsetleft/right temp have length setsize
            for k=1:setsize

                if Sidechoice==0 %only assign memsetleft when test on that side
                    memsetleft{k}=memsetlefttemp{ordernow(k)};
                else
                    memsetright{k}=memsetrighttemp{ordernow(k)};
                end

            end

            %%choose test item
            kstim=999999; %just to make sure we are not messing up  
            if Sidechoice==0 %left
                if old==1 %The response should be learned
%                     disp('order:');disp(ordernow(serpos));disp(class(ordernow(serpos)));
%                     disp('serpos:');disp(serpos);disp(class(serpos));
                    %m=ordernow(serpos);
                    kstim=ordertemp(ordernow(serpos));%EX:=2;serpos is randomly chosen prviously, this argument give the stimulus num 
                    teststim=memsetleft{serpos}; %memsetleft include CN and AN
                
                else % test new item from other side
                    newca=rand; %the variable to decide whether test CM or AN
                    if newca<0.5 % test CM new
                        k=randi(8); %8 for cm one side
                        kstim=k+8;
                        teststim=CM_right{k};
                    else
                        antot=antot+1;
                        kstim=16+antot;
                        teststim=fullset{kstim};
                    end
                end

            else %right
                if old==1 %The response should be learned
                    kstim=ordertemp(ordernow(serpos));%EX:=2;serpos is randomly chosen prviously, this argument give the stimulus num 
                    teststim=memsetright{serpos};
                
                else % test new item from other side
                    newca=rand;
                    if newca<0.5
                        k=randi(8); %8 for cm one side
                        kstim=k+8;
                        teststim=CM_left{k};
                    else
                        antot=antot+1;
                        kstim=16+antot;
                        teststim=fullset{kstim};
                    end


                end
            end


            if kstim>8
            % Store if this trial is cm or an. 1:cm | 0:an
                is.cm=1;
            else
                is.cm=0;
            end


            %%
            %   present the fixation point
            %
            Screen('TextSize',wind1,fixation_size);
            Screen('DrawText',wind1,fixation,rect(3)/2-textbounds_fixation(3)/2,rect(4)/2-textbounds_fixation(4)/2);
            Screen('Flip',wind1);
            if debug~=1
                WaitSecs(.5 );
            end  
            Screen('Flip',wind1);
            if debug~=1
                WaitSecs(.1);
            end

            %%
            %   present the memory set
            %
            
            if Sidechoice<0.5 %%%<0.5==left
                for i=1:setsize

                    img_data=imread([image_location memsetleft{i}]);
                    picture=Screen('MakeTexture',wind1,img_data);
                    Screen('DrawTexture',wind1,picture,[],leftimagecoords);
                    Screen('Close',picture);
                    Screen('Flip',wind1);
                    WaitSecs(present_rate);
                    Screen('Flip',wind1);
                    WaitSecs(isi);
                end
            else
                for i=1:setsize

                        img_data=imread([image_location memsetright{i}]);


                picture=Screen('MakeTexture',wind1,img_data);
                Screen('DrawTexture',wind1,picture,[],rightimagecoords);
                Screen('Close',picture);
                Screen('Flip',wind1);
                WaitSecs(present_rate);
                Screen('Flip',wind1);
                WaitSecs(isi);
               end


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
                if debug==0
                    [keydown secs keycode]=KbCheck;
                    if keydown==1
                        flag=1;
                    end
                end
            end
            %%
            %          present probe and collect legal response
            %
            if flag == 0
                img_data=imread([image_location teststim]);
                picture=Screen('MakeTexture',wind1,img_data);

                if Sidechoice<0.5
                    Screen('DrawTexture',wind1,picture,[],leftimagecoords);
                else
                    Screen('DrawTexture',wind1,picture,[],rightimagecoords);
                end

                Screen('Flip',wind1);
                Screen('Close',picture);
                start=GetSecs;
                legal=0;
                if debug==0
                    while legal == 0
                        [keydown secs keycode] = KbCheck;
                        key=KbName(keycode);
                        if any(strcmp(key,legalkeys))
                            rt=secs-start;
                            legal=1;
                        end
                    end
                else 
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
                if old == 1 && resp == 2 && Sidechoice==0 %%%when press f; the item is old; and the side is left
                    corr=1;
                end
                if old==1 && resp==1 && Sidechoice==1 %%%when press j; the item is new; and the side is right
                    corr=1;
                end
                if old==2 && resp==1 && Sidechoice==0
                    corr=1;
                end
                if old == 2 && resp == 2 && Sidechoice==1
                    corr=1;
                end
                %%
                %             present feedback
                %
                Screen('TextSize',wind1,textsize)
                if corr == 1
                    percent_correct=percent_correct+1;
                    Screen('DrawText',wind1,text_correct,rect(3)/2-textbounds_correct(3)/2,rect(4)/2-textbounds_correct(4)/2);
                    Screen('Flip',wind1);
                else
                    Screen('DrawText',wind1,text_incorrect,rect(3)/2-textbounds_incorrect(3)/2,rect(4)/2-textbounds_incorrect(4)/2);
                    Screen('Flip',wind1);
                end
                if debug~=1
                    WaitSecs(1);
                end
                Screen('Flip',wind1);
                %%
                %               record results
                %
                tot_trials=tot_trials+1;
                block_store(tot_trials)=block;
                trial_store(tot_trials)=trial;
                serpos_store(tot_trials)=serpos;
                lag_store(tot_trials)=lag;
                trial_type_store(tot_trials)=old;
                %study_set_store(tot_trials)=studyCon;
                test_set_stor(tot_trials)=testCon;
                setsize_store(tot_trials)=setsize;
                resp_store(tot_trials)=resp;
                rt_store(tot_trials)=rt;
                corr_store(tot_trials)=corr;
                probe_store(tot_trials)=kstim;
                side_store(tot_trials)=Sidechoice;
                stimkind_store(tot_trials)=is.cm;  %1: cm; 0:an
                                                    %the current test prob or the current tiral is cm or an 

                set_store=9999999999;
                for i=1:setsize
                    set_store(tot_trials,i)=ms_stimnum(i); %changed
                end
                %%
                %               write to output text file
                %
                fprintf(fid,'%5d',block_store(tot_trials),trial_store(tot_trials),trial_type_store(tot_trials), test_set_stor(tot_trials),...
                    setsize_store(tot_trials),serpos_store(tot_trials),lag_store(tot_trials),resp_store(tot_trials),...
                    corr_store(tot_trials),side_store(tot_trials));
                fprintf(fid,'%10d',round(1000*rt_store(tot_trials)));
                fprintf(fid,'%5d',probe_store(tot_trials));
                fprintf(fid,'%5d',stimkind_store(tot_trials));
                fprintf(fid,fn);
                for i=1:setsize
                    fprintf(fid,'%5d',set_store(tot_trials,i));
                end
                fprintf(fid,'\n');


            else
                %%
                %               warn subject for premature key press
                %
                Screen('DrawText',wind1,warning,rect(3)/2-textbounds_warning(3)/2,rect(4)/2-textbounds_warning(4)/2);
                Screen('Flip',wind1);
                if debug~=1
                    WaitSecs(5);
                end
            end
        end %end trial
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
    save(consistvars);
    debriefing(wind1);
    Screen('DrawText',wind1,thanks,rect(3)/2-textbounds_thanks(3)/2,rect(4)/2-textbounds_thanks(4)/2);
    Screen('Flip',wind1);
    if debug~=1
        WaitSecs(.5);
    end
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

