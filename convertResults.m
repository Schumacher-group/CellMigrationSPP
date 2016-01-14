% convert some old saved results
clear
cd('results')
files = dir('*.mat');
for fileCtr= 1:length(files)
    load(files(fileCtr).name)
    cells = cells(:,1:3,:);
    if exist('selfAlign','var')
        save(files(fileCtr).name,'cells','T','N','L','alpha','beta','selfAlign','bcs')
    elseif exist('bcs','var')
        save(files(fileCtr).name,'cells','T','N','L','alpha','beta','bcs')
    else
        save(files(fileCtr).name,'cells','T','N','L','alpha','beta')
    end   
end
cd ..