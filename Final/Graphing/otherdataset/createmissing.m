clear;
clc;
load('TwoMoon.mat');

data = X;
num_view = length(data);

% For missing 0
percentage =10;
miss0={}; 
for i = 1:10
    temp =  [];
    for iv = 1:num_view
        view_data = data{iv}'; % get view data
        size_of_data = size(view_data); % get size of data
        rep_index = ones(size_of_data(2),1); % construct ones matrix;
        % inCol= randi(size(tmpda,2),1) % get random number ; 
        % number to remove percent calculation
        sample_number = size_of_data(2);
        number_to_remove = ceil((sample_number/100)*percentage);
        generate_index = randperm(sample_number,number_to_remove);
        rep_index(generate_index) = 1; % 400x1
        temp(:,iv) = rep_index;    
    end
    miss0{i} = temp;   
end



% For missing 10 
percentage =10;
miss10={}; 
for i = 1:10
    temp =  [];
    for iv = 1:num_view
        view_data = data{iv}'; % get view data
        size_of_data = size(view_data); % get size of data
        rep_index = ones(size_of_data(2),1); % construct ones matrix;
        % inCol= randi(size(tmpda,2),1) % get random number ; 
        % number to remove percent calculation
        sample_number = size_of_data(2);
        number_to_remove = ceil((sample_number/100)*percentage);
        generate_index = randperm(sample_number,number_to_remove);
        rep_index(generate_index) = 0; % 400x1
        temp(:,iv) = rep_index;    
    end
    miss10{i} = temp;   
end

%clearvars -except data miss10 miss20 miss30 miss40 miss50 miss60 miss70 miss80 truth;
clearvars -except X miss0 miss10 miss20 miss30 miss40 miss50 miss60 miss70 miss80 y0;
save('./otherdataset/TwoMoonMissing.mat');  % ******** CHANGE name Related to database;
