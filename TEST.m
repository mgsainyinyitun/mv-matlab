% Load Data
clear
clc
load('./database/buaaRnSp.mat');
f = 8;
X = data; % complete data
folds = miss10;
ind_folds = folds{f};
truthF = truth;
numClust = length(unique(truthF));
num_view = length(X);

% Normalize Data
% make incomplete data by removing incomplete sample
% fill with average value
for iv = 1:num_view
    X1 = X{iv}';
    X1 = NormalizeFea(X1,1);
    ind_0 = find(ind_folds(:,iv) == 0);
    X1(ind_0,:) = missing;      % replace with NaN 
    Y{iv} = X1';                % transfer to Y   % m x n => m = feature, n = sample   
    avg{iv} = mean(Y{iv},2,'omitnan');
    Y{iv} = fillmissing(Y{iv}','constant',avg{iv});
    Y{iv} = Y{iv}';
end
clear X X1 W1
X = Y; % complete data with average value filled
clear Y 









