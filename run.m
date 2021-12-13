% Load Data
clear
clc

%load('./database/bbcsport4vbigRnSp.mat');
%load('./database/100Leaves.mat');
%load('./database/ORL.mat');
load('./database/mfeatRnSp.mat');
%load('./database/WebKB.mat');
f = 1;
X = data; % complete data
folds = miss10;
ind_folds = folds{f};
truthF = truth;  
numClust = length(unique(truthF));
num_view = length(X);
% Normalize Data
% make incomplete data
for iv = 1:num_view
    X1 = X{iv}';
    X1 = NormalizeFea(X1,1);
    ind_0 = find(ind_folds(:,iv) == 0);
    X1(ind_0,:) = [];    % incomplete data  
    Y{iv} = X1';  % incomplete data           
    W1 = eye(size(ind_folds,1));
    W1(ind_0,:) = [];
    G{iv} = W1; % G1,G2,G3,G4 % four different view                                         
end

clear X X1 W1
X = Y; % incomplete
clear Y 

% initial Z
for iv = 1:num_view
    options = [];
    options.NeighborMode = 'KNN';
    options.k = 3;
    options.WeightMode = 'Binary';
    Z1 = constructW(X{iv}',options);
    Z_ini{iv} = full(Z1); % Z_ini1,Z_ini2,Z_ini3,Z_ini4
    clear Z1;
end

% initialize H
for iv = 1:num_view
% ===Initial Vi=== k*ni
    ni_num(iv) = size(X{iv},2);% num of existing sample
    Htmp = litekmeans(X{iv}',numClust,'MaxIter',100);% Vtmp: ni*1
    tmp = zeros(ni_num(iv),numClust);
    tmp(sub2ind(size(tmp),[1:ni_num(iv)],Htmp'))=1;% ni*1->ni*k
    H{iv} = tmp';
end

% compute inverse GPU vers
for iv = 1:num_view
   Htemp = H{iv};
   Htemp = gpuArray(Htemp);
   invtemp = inv(Htemp*Htemp');
   invtemp =gather(invtemp);
   invH{iv} = invtemp;
end

max_iter=200;

% update 
%% update NMF
gamma = 0.001;
[HF,WF,Z,F,S]=algorithm(X,H,Z_ini,G,invH,num_view,truthF,gamma,max_iter);

% complete graph
% for iv=1:num_view
%     HF{iv} = HF{iv}*G{iv}; % 5x104 * 104x116
% end
% avgH=0;
% for iv=1:num_view
%     avgH = avgH + HF{iv};
% end
% avgH = avgH/num_view;




%% Graph learning
% complete graph
% average graph
% U = 0;
% for iv=1:num_view
%     Z{iv} = G{iv}'*Z{iv}*G{iv};
%     U = U + Z{iv};
% end
% U = U/num_view;

new_F = S;
norm_mat = repmat(sqrt(sum(new_F.*new_F,2)),1,size(new_F,2));
for i = 1:size(norm_mat,1)
    if (norm_mat(i,1)==0)
        norm_mat(i,:) = 1;
    end
end
new_F = new_F./norm_mat; 

%% fusion

% normalize Graph(optional)




%% cluster with kmean
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




%% cluster graph with spectral clustering
disp('Spectral clustering');
metric = spcclust(S', numClust, truth)








