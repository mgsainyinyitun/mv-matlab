
clear
clc
% Load Data
%load('./database/bbcsport4vbigRnSp.mat');
%load('./database/100Leaves.mat');
%load('./database/ORL.mat');
load('./database/mfeatRnSp.mat');
%load('./database/WebKB.mat');
%load('./database/orlRnSp.mat');
%load('./database/caltech7.mat');
%load('./database/buaaRnSp.mat');
%load('./database/Mfeat.mat');
%load('./database/3sources.mat');

%load('./databaseallmissing/bbcsport4vbigRnSp.mat');

f = 1;
X = data; % complete data
folds = miss10;
ind_folds = folds{f};
truthF = truth;
numClust = length(unique(truthF));
num_view = length(X);
% construct incomplete data and index matrix
for iv = 1:num_view
    X1 = X{iv}';
    ind_0 = find(ind_folds(:,iv) == 0);
    X1(ind_0,:) = [];         % incomplete data  
    Y{iv} = X1';              % incomplete data           
    W1 = eye(size(ind_folds,1));
    W1(ind_0,:) = [];
    G{iv} = W1; % G1,G2,G3,G4 % four different view % ni x n      
end


c = length(unique(truthF));
y0 = truthF;

tic
[U,F,obj_value] = gmc_fusion(Y, c,G); % c: the # of clusters
toc
S = U;


%% cluster with kmean
disp('kmean clustering');
new_F = F;
norm_mat = repmat(sqrt(sum(new_F.*new_F,2)),1,size(new_F,2));
for i = 1:size(norm_mat,1)
    if (norm_mat(i,1)==0)
        norm_mat(i,:) = 1;
    end
end
new_F = new_F./norm_mat; 

repeat = 5;
for iter_c = 1:repeat
    pre_labels    = kmeans(real(new_F),numClust,'emptyaction','singleton','replicates',20,'display','off');
    result_LatLRR = ClusteringMeasure(truthF, pre_labels);       
    AC(iter_c)    = result_LatLRR(1)*100;
    MIhat(iter_c) = result_LatLRR(2)*100;
    Purity(iter_c)= result_LatLRR(3)*100;
end
mean_ACC = mean(AC)
mean_NMI = mean(MIhat)
mean_PUR = mean(Purity)
%% ----------------------
obj_value = obj_value';

