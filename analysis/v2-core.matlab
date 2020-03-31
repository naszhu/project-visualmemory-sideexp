


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

    %% Define some variables
    wss={[0.12346,0.4,0.47654],[0.12346,0.8,0.07654],[0.12346,0.87654,0]}
    probchart_old={[0.9,0.1],[0.9,repelem(1/30,3)],[0.9,repelem(1/50,5)];...
    [1/2,1/2],[repelem(1/4,4)],[0.228169,repelem(0.154366,4)];...
    [1/2,1/2],[repelem(1/3,3)],[,0]};
    prob_c_d=[1/3,1/6]
    probchart_new={[prob_c_d,repelem(1/10,5)],[prob_c_d,repelem(1/6,3)],[prob_c_d,1/2];...
    [prob_c_d,repelem(1/10,5)],[prob_c_d,repelem(1/6,3)],[prob_c_d,1/2];...
    [1/3,1/8,repelem(1/8,4)],[1/3,1/6,repelem(1/4,2)],[0]
    }



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
            foilist(1)=cmneg_l; %=3
            foilist(2)=cmneg_h; %=4
 

            %% Memory list           
            testkind=["A","B","nAB"];

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
                cmassign=1;
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
                    serpos=0;
                    lag=0
                end
            else % not with special CM+,then all 50%
                if whattest<0.5
                    old=2;
                end
            end

            %% Draw test item
            if old==1
                testitem = randsample(memlist,1,true,probchart_new{cmassign,ss}); %there is no 3,6 
            else % New
                testitem = randsample(foilist,1,true,probchart_new{cmassign,ss});
            end

            serpos = find(memset==testitem);
            lag=setsize-serpos+1;

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
        end

