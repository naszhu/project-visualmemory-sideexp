        
    %blockwise:
    antot=0;
        is.cm=1; %whether this trial is cm or an
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

            %% CM or AN trials
            memsetleft={};%empty memory set
            memsetright={}; 
            memset={};

            %%choose memory set
            if  is.cm==1  %cm trials   

                %CM_left or CM_right stores CM pools in order, and so order represents which stimuli is drawn in there
                order=randperm(8); %ex: 2     3     1     7     5     8     4     6
                for k=1:setsize
                    memsetleft{k}=CM_left{order(k)}; %MEMleft_cm; previously: CM_lift{i}=pictures{i};
                    memsetright{k}=CM_right{order(k)}; %memright_cm
                end
                is.cm=0; %test an next time

            else % an trials
                for k=1:setsize
                    antot=antot+1;
                        %16 are #of stimuli taken by CM; 
                    memset{k}=pictures{16+antot}; %CHANGE if need; all AN
                end
                is.cm=1; %test cm next time

            end
            
            %%choose test item
            Sidechoice=round(rand);
            if Sidechoice==0 %left
                if old==1 %The response should be learned
                     
                    if is.cm==1
                        kstim=order(serpos);%EX:=2;serpos is randomly chosen prviously, this argument give the stimulus num 
                        teststim=memsetleft{serpos};
                    else
                        teststim=memset{serpos};
                        kstim=16+antot-lag+1; %ex: 16+8-7+1=18 
                    end
                
                else % test new item from other side

                    if is.cm==1    
                        k=randi(8);
                        kstim=k+8; %k+8 because the index for first 8 is taken by CM_left
                        teststim=CM_right{k};
                    else 
                        antot=antot+1; %ex: (16+8)+1
                        teststim=pictures{16+antot};
                    end

                end

            else %right
                if old==1 %The response should be learned
                     
                    if is.cm==1
                        kstim=order(serpos);%EX:=2;serpos is randomly chosen prviously, this argument give the stimulus num 
                        teststim=memsetright{serpos};
                    else
                        teststim=memset{serpos};
                        kstim=16+antot-lag+1; %ex: 16+8-7+1=18 
                    end
                
                else % test new item from other side

                    if is.cm==1    
                        k=randi(8);
                        kstim=k+8;
                        teststim=CM_left{k};
                    else 
                        antot=antot+1; %ex: (16+8)+1
                        teststim=pictures{16+antot};
                    end

                end

               