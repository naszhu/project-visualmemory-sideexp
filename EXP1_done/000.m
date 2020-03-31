%%mixed cmvman experiment
% mixedcmvman.m
% replicate rui's experiment
% 3 CM+, 3CM-, 6VM, unlimited AN
% 01/02/19
%
clear all;
image_location=[pwd '\images\'];
data_location=[pwd '\data\'];
subid=input(' subject # ');
filename=[data_location 'cmvman' num2str(subid) '.txt'];
variedvars=[data_location 'cmvman' num2str(subid)];
if ~exist(filename)
    for block=1:nblocks

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
            setsize=set_sizes(index); %random select a set size

            ss=setsize;  % abbrev. for convenience; ss is the current set size

            %% generate memory set: give what sequence we want and what stimuli we want
            % memset is a number's set 
            ordcm=randperm(n/2); % perm(3) 3 cm target's random sequence
            for i=1:ss/2
                memlist(i)=cmpos(ordcm(i));
            end
            ordvm=randperm(nvm); % perm(4)
            for i=1:ss/2 %only choose 1/2/3 vm item
                memlist(i+ss/2)=vm(ordvm(i));
            end

            % Replace cmpos with probability 0.1 
            ifcms=randi(10); % 0.1 probability of being chosen
            isspc=0; % for later choosing test item
            if ifcms==1
                memlist(1) = cmspe;
                isspc=1;
            end

            % Randomize memset
            ordss=randperm(ss);
            memset=[]
            for i=1:ss
                memset(i)=memlist(ordss(i)); % permutation memset
            end

            %% Test item:
            % randomly choose old/new 
            %
            old=1; %old
            if rand<.5 % 0.5 probability of being new
                old=2;
                serpos=0;
                lag=0;
            end

            %% force to always test special if it appears
            if isspc==1
                serpos = find(memset==cmspe); %find current potition of cspc
                testitem =  memset(serpos); %should equal to 11
                cond = 1; % cmSpc

            else  % if no special item
                if old==1 %gives condition and choose current test item

                    %% give weight to each item in memeory set
                    w=[]
                    for i=1:ss 
                        if memset(i)<= n/2 %i<=3
                            if ss==2
                                w(i)=2/3;
                            elseif ss==4
                                w(i)=1/3;
                            else %ss==6
                                w(i)=2/9;
                            end
                        elseif memset(i)>n % out of 7,8,9,10
                            if ss==2
                                w(i)=1/3;
                            elseif ss==4
                                w(i)=1/6;
                            else %ss==6
                                w(i)=1/9;
                            end
                        end
                    end

                    %% assign random test item
                    testitem = randsample(memset,1,true,w);
                    serpos = find(memset==testitem);
                    lag=setsize-serpos+1;

                    if testitem <=n/2
                        cond = 2; %cm
                    elseif testitem >n
                        cond = 3; %vm
                    end

                else   %choose the new testitem
                    %% Choose cm foil or vm foil
                    if rand<0.5 %test cm foil
                        wnew=[0.5, 0.25, 0.25]; %this can move to the front
                        testitem=randsample(cmneg,1,true,wnew);
                        if testitem==7 %***, which item assigned with .5 prob
                            cond = 1; %special cm
                        else 
                            cond = 2; %cm
                        end
                    else %test vm foil
                        testitem=vm(ordvm(ss/2+1)); 
                        % use ss/2 in vm, so ss/2+1 is the next availiable vm
                        cond = 3; %vm
                    end
                end
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
                img_data=imread([image_location pictures{memset(i)}]);
                picture=Screen('MakeTexture',wind1,img_data);
                Screen('DrawTexture',wind1,picture,[],[loc(1) loc(2) loc(3) loc(4)]);
                Screen('Close',picture);
                Screen('Flip',wind1);
                WaitSecs(present_rate);
                Screen('Flip',wind1)
                WaitSecs(isi);
            end
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
                img_data=imread([image_location pictures{testitem}]);
                picture=Screen('MakeTexture',wind1,img_data);
                Screen('DrawTexture',wind1,picture,[],[loc(1) loc(2) loc(3) loc(4)]);
                Screen('Flip',wind1);
                Screen('Close',picture);
                start=GetSecs;
                legal=0;
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
                
                %
                % record response information
                %
                resp_store(tot_trials)=resp;
                rt_store(tot_trials)=rt;
                corr_store(tot_trials)=corr;
                
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
    Screen('DrawText',wind1,thanks,rect(3)/2-textbounds_thanks(3)/2,rect(4)/2-textbounds_thanks(4)/2);
    Screen('Flip',wind1);
    WaitSecs(.5);
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
else
    disp('Error: the file already exists!')
    end
end





