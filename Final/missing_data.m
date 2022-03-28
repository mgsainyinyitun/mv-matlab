clear;
clc;
load('./database/caltech7.mat');

num_view = length(data);

% For missing 60 
percentage = 60;
miss60={}; 
for i = 1:10
    temp =  [];
    for iv = 1:num_view
        view_data = data{iv}; % get view data
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
    miss60{i} = temp;   
end

% For missing 70
percentage = 70;
miss70={}; 
for i = 1:10
    temp =  [];
    for iv = 1:num_view
        view_data = data{iv}; % get view data
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
    miss70{i} = temp;   
end



% For missing 80
percentage = 80;
miss80={}; 
for i = 1:10
    temp =  [];
    for iv = 1:num_view
        view_data = data{iv}; % get view data
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
    miss80{i} = temp;   
end

clearvars -except data miss10 miss20 miss30 miss40 miss50 miss60 miss70 miss80 truth;

if ~exist('databaseallmissing','dir')
mkdir databaseallmissing;
end

save('databaseallmissing\caltech7.mat');  % ******** CHANGE name Related to database;
