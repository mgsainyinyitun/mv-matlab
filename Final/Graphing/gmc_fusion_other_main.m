%% Toy examples on Two-Moon data set and Three-Ring data set
% Graph-based Multi-view Clustering (GMC)
% 
%%
clc;  close all; clear all;
currentFolder = pwd;
addpath(genpath(currentFolder));

%% Load toy data set
datadir = 'otherdataset/';
m = 2; % number of views
X = cell(1,m);
dataname = input('Input the name of data set: (TwoMoonMissing or ThreeRingMissing)\n','s');
while(1)
    if strcmp(dataname,'TwoMoonMissing')
        dataf = [datadir, dataname];
        c = 2;
        load(dataf);
        break;
    elseif strcmp(dataname,'ThreeRingMissing')
        dataf = [datadir, dataname];
        c = 3;
        load(dataf);
        flag = 1;
        break;
    else
        dataType = input('Please only input TwoMoonMissing or ThreeRingMissing\n','s');
    end;
    
end

folds = miss0;
ind_folds = folds{1}; % 10% missing

%% Transport X
% for iv = 1:m
%     X{iv} = X{iv}';
% end

normData = 1;

%% Normalization: Z-score
if normData == 1
    for i = 1:m % 1-2
        dim = size(X{i});
        for  j = 1:dim(2) % 1-200
            normItem = std(  X{i}(:,j) );
            if (0 == normItem)
                normItem = eps;
            end
            X{i}(:,j) = (X{i}(:,j)-mean(X{i}(:,j)))/(normItem);
        end
    end
end  

%% Calculate For G
for iv = 1:m
    X1 = X{iv}';
    ind_0 = find(ind_folds(:,iv) == 0);
    X1(ind_0,:) = [];         % incomplete data  
    Y{iv} = X1';              % incomplete data           
    W1 = eye(size(ind_folds,1));
    W1(ind_0,:) = [];
    G{iv} = W1; % G1,G2,G3,G4 % four different view % ni x n         
end

%% Call gmc_fusion algorithm
num = size(X{1},1); % the number of samples
data = cell(1,m);
for i = 1:m
    data{i} = X{i}';
end
% [predY, U, S0, S0_initial, F, evs] = gmc_fusion(data, c, 1, 0);
[U,F,S0,S0_initial,obj_value] = gmc_fusion2(data,c,G);

% cluster with k-mean to get predY
S = U;
new_F = F;
norm_mat = repmat(sqrt(sum(new_F.*new_F,2)),1,size(new_F,2));
for i = 1:size(norm_mat,1)
    if (norm_mat(i,1)==0)
        norm_mat(i,:) = 1;
    end
end
new_F = new_F./norm_mat; 
numClust = length(unique(y0));
predY = kmeans(real(new_F),numClust,'emptyaction','singleton','replicates',20,'display','off');

metric = CalcMeasures(y0(:,1), predY);
fprintf('Data set %s-> ACC:%.4f\tNMI:%.4f\tARI:%.4f\terror_cnt:%d\n',dataname,metric(1),metric(2),metric(3),metric(4));

markerSize = 20;
%% Original data
for v = 1:m
    lab = y0(:,v);
    cLab = unique(lab);
    figure; 
    plot(X{v}(:,1),X{v}(:,2),'.k', 'MarkerSize', markerSize); hold on;
    plot(X{v}(lab==cLab(1),1),X{v}(lab==cLab(1),2),'.r', 'MarkerSize', 20); hold on;
    plot(X{v}(lab==cLab(2),1),X{v}(lab==cLab(2),2),'.', 'MarkerSize', markerSize); hold on;
    if flag
        plot(X{v}(lab==cLab(3),1),X{v}(lab==cLab(3),2),'.', 'Color', [79 79 79]/255, 'MarkerSize', markerSize); hold on;
    end;
%     set(gca,'xlim',[-1.7,1.7],'xtick',[-1.5:0.5:1.5]) % set x-axis
%     set(gca,'ylim',[-1.2,1.2],'ytick',[-1:0.5:1]) % set y-axix
    set(gca,'FontName','Times New Roman','FontSize',20,'LineWidth',1.2);
%     str = 'The orighnal data';
%     titlename = sprintf('%s: View-%d',str,v);
%     t = title(titlename);
%     set(t,'FontName','Times New Roman','FontSize',20.0);
    axis equal;
end;

%% Original connected graph with probabilistic neighbors, line width denotes similarity
S1 = cell(1,m);
for v = 1:m
    S1{v} = S0_initial{v};
    lab = y0(:,v);
    cLab = unique(lab);
    figure; 
    plot(X{v}(:,1),X{v}(:,2),'.k', 'MarkerSize', markerSize); hold on;
    plot(X{v}(lab==cLab(1),1),X{v}(lab==cLab(1),2),'.r', 'MarkerSize', markerSize); hold on;
    plot(X{v}(lab==cLab(2),1),X{v}(lab==cLab(2),2),'.', 'MarkerSize', markerSize); hold on;
    if flag
        plot(X{v}(lab==cLab(3),1),X{v}(lab==cLab(3),2),'.', 'Color', [79 79 79]/255, 'MarkerSize', markerSize); hold on;
    end;
    for ii = 1 : num
        for jj = 1 : ii
            weight = S1{v}(ii, jj);
            if weight > 0.5
                plot([X{v}(ii, 1), X{v}(jj, 1)], [X{v}(ii, 2), X{v}(jj, 2)], '-', 'Color', [0 197 205]/255, 'LineWidth', 5*weight), hold on;
            end;
        end;
    end;
%     set(gca,'xlim',[-1.7,1.7],'xtick',[-1.5:0.5:1.5]) % set x-axis
%     set(gca,'ylim',[-1.2,1.2],'ytick',[-1:0.5:1]) %set y-axis
    set(gca,'FontName','Times New Roman','FontSize',20,'LineWidth',1.2);
%     str = 'Connected graph with probabilistic neighbors';
%     titlename = sprintf('%s: View-%d',str,v);
%     t = title(titlename);
%     set(t,'FontName','Times New Roman','FontSize',20.0);
    axis equal;
end;

%% the learned graph of each view by GMC, line width denotes similarity
S2 = cell(1,m);
for v = 1:m
    S2{v} = S0{v};
    lab = y0(:,v); % label
    cLab = unique(lab); % 2 
    figure; 
    plot(X{v}(:,1),X{v}(:,2),'.k', 'MarkerSize', markerSize); hold on;
    plot(X{v}(lab==cLab(1),1),X{v}(lab==cLab(1),2),'.r', 'MarkerSize', markerSize); hold on;
    plot(X{v}(lab==cLab(2),1),X{v}(lab==cLab(2),2),'.', 'MarkerSize', markerSize); hold on;
    if flag
        plot(X{v}(lab==cLab(3),1),X{v}(lab==cLab(3),2),'.', 'Color', [79 79 79]/255, 'MarkerSize', markerSize); hold on;
    end;
    for ii = 1 : num
        for jj = 1 : ii
            weight = S2{v}(ii, jj);
            if weight > 0.05
                plot([X{v}(ii, 1), X{v}(jj, 1)], [X{v}(ii, 2), X{v}(jj, 2)], '-', 'Color', [0 197 205]/255, 'LineWidth', 5*weight), hold on;
            end;
        end;
    end;
%     set(gca,'xlim',[-1.7,1.7],'xtick',[-1.5:0.5:1.5]) % set x-axis
%     set(gca,'ylim',[-1.2,1.2],'ytick',[-1:0.5:1]) % set y-axis
    set(gca,'FontName','Times New Roman','FontSize',20,'LineWidth',1.2);
%     str = 'Connected graph with probabilistic neighbors';
%     titlename = sprintf('%s: View-%d',str,v);
%     t = title(titlename);
%     set(t,'FontName','Times New Roman','FontSize',20.0);
    axis equal;
end;

%% the learned unified graph by GMC, line width denotes similarity
U2 = U;
for v = 1:m
    lab = y0(:,v);
    cLab = unique(lab);
    figure; 
    plot(X{v}(:,1),X{v}(:,2),'.k', 'MarkerSize', markerSize); hold on;
    plot(X{v}(lab==cLab(1),1),X{v}(lab==cLab(1),2),'.r', 'MarkerSize', markerSize); hold on;
    plot(X{v}(lab==cLab(2),1),X{v}(lab==cLab(2),2),'.', 'MarkerSize', markerSize); hold on;
    if flag
        plot(X{v}(lab==cLab(3),1),X{v}(lab==cLab(3),2),'.', 'Color', [79 79 79]/255, 'MarkerSize', markerSize); hold on;
    end;
    for ii = 1 : num
        for jj = 1 : ii
            weight = U2(ii, jj);
            if weight > 0
                plot([X{v}(ii, 1), X{v}(jj, 1)], [X{v}(ii, 2), X{v}(jj, 2)], '-', 'Color', [0 197 205]/255, 'LineWidth', 5*weight), hold on;
            end;
        end;
    end;
%     set(gca,'xlim',[-1.7,1.7],'xtick',[-1.5:0.5:1.5]) % set x-axis
%     set(gca,'ylim',[-1.2,1.2],'ytick',[-1:0.5:1]) % set y-axix
    set(gca,'FontName','Times New Roman','FontSize',20,'LineWidth',1.2);
%     str = 'Learnt connected graph';
%     titlename = sprintf('%s: View-%d',str,v);
%     t = title(titlename);
%     set(t,'FontName','Times New Roman','FontSize',20.0);
    axis equal;
end;
