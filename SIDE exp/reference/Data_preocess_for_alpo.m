files = dir('*.mat');
for i=1:length(files)
    eval(['load ' files(i).name ' -ascii']);
end

%[filename,foldername] = uigetfile('*.mat','Choose one subject:') 
 files=dir('*.mat');

filename = 'color_all.txt';
filepath = cd;
file = fullfile(filepath, filename);
fid = fopen(file, 'wt');

s=length(files);
for j=1:s
    filename=files(j).name
load(filename)
    %check with setsize 1 2 or 4
    correct_for_setsize_1=0;
    correct_for_setsize_2=0;
    correct_for_setsize_4=0;
    trials_for_setsize_1=0;
    trials_for_setsize_2=0;
    trials_for_setsize_4=0;
     serial_setsize_2=[0,0];
    serial_setsize_4=[0,0,0,0];
    

    for i=1:200
        if setsize_store(i)==1
            trials_for_setsize_1=trials_for_setsize_1+1;
            if memset_store(i,1)==response_store(i,1)
                correct_for_setsize_1=1+correct_for_setsize_1;
            end
        end
        if setsize_store(i)==2
            for n=1:2
                trials_for_setsize_2=trials_for_setsize_2+1;
                if memset_store(i,n)==response_store(i,1)|| memset_store(i,n)==response_store(i,2)
                correct_for_setsize_2=1+correct_for_setsize_2;
                serial_setsize_2(n)=serial_setsize_2(n)+1;
                end
            end
        end
        if setsize_store(i)==4
            for n=1:4
                trials_for_setsize_4=trials_for_setsize_4+1;
                if memset_store(i,n)==response_store(i,1)|| memset_store(i,n)==response_store(i,2)||memset_store(i,n)==response_store(i,3)||memset_store(i,n)==response_store(i,4)
                     correct_for_setsize_4=1+correct_for_setsize_4;
                serial_setsize_4(n)=serial_setsize_4(n)+1;
                end
            end
        end
    end
    %percentage_correct_setsize1=correct_for_setsize_1/trials_for_setsize_1;
    %percentage_correct_setsize2=correct_for_setsize_2/trials_for_setsize_2;
    %percentage_correct_setsize4=correct_for_setsize_4/correct_for_setsize_4;
    
   fprintf(fid,'%5d',subid,correct_for_setsize_1,trials_for_setsize_1);
    for m=1:2
        fprintf(fid,'%5d',serial_setsize_2(m));
    end
    fprintf(fid,'%5d',trials_for_setsize_2/2);
    for m=1:4
        fprintf(fid,'%5d',serial_setsize_4(m));
    end
    fprintf(fid,'%5d',trials_for_setsize_4/4);
    fprintf(fid,'\n');
end
    fclose(fid);


