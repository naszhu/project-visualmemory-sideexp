            if ss==2 %ss 2 or 4
                whatcm=randsample(testkind,1,true,wss2);
                if whatcm=="A" 
                    memlist(1)=cmpos_s;
                    memlist(2)=randsample(vm,1);
                elseif whatcm=="B" 
                    memlist(1)=cmpos_r;
                    memlist(2)=randsample(vm,1);
                elseif whatcm=="nAB"
                    ordvm=randperm(nvm);
                    for i:ss
                        memlist(i)=vm(ordvm(i));
                    end
                else
                    disp("Error ss manipulation 1");
                end
            elseif ss==4
                whatcm=randsample(testkind,1,true,wss4);
                if whatcm=="A"
                    mmlist(1)=cmpos_s;
                    ordvm=randperm(nvm);
                    for i=1:(ss-1) %1:3
                        memlist(i+1)=vm(ordvm(i)); %i+1 because 1 with A
                    end
                elseif whatcm=="B"
                    mmlist(1)=cmpos_r;
                    ordvm=randperm(nvm);
                    for i=1:(ss-1)
                        memlist(i+1)=vm(ordvm(i));
                    end
                elseif whatcm=="nAB"
                    ordvm=randperm(nvm);
                    for i=1:ss %1:4
                        memlist(i)=vm(ordvm(i));
                    end
                else 
                    disp("Error ss manipulation 2");
                end
            else %ss==6
                whatcm=randsample(testkind,1,true,wss6);
                if whatcm=="A"
                    mmlist(1)=cmpos_s;
                    

            end