
function [U, F] = gmc_fusion(X, c,G, lambda, normData) 

%% input:
% X{}: multi-view dataset, each cell is a view, each column is a data point
% c: cluster number
% lambda: parameter (default 1)
%% output:
% U: the learned unified matrix
% F: the embedding representation
NITER = 20;
zr = 10e-11;
pn = 30;              % number of neighbours for constructS_PNG
islocal = 1;           % only update the similarities of neighbors if islocal=1
if nargin < 4
    lambda = 1;
end;
if nargin < 5
    normData = 1;
end;

num = size(X{1},2);         % number of incomplete instances
numC = size(G{1},2);        % number of complete instances
m = length(X);              % number of views
%% Normalization: Z-score
if normData == 1
    for i = 1:m
        for  j = 1:num
            normItem = std(X{i}(:,j));
            if (0 == normItem)
                normItem = eps;
            end
            X{i}(:,j) = (X{i}(:,j)-mean(X{i}(:,j)))/(normItem);
        end
    end
    
end  
%% initialize Z0: Constructing the SIG matrices
Z0 = cell(1,m);
for i = 1:m
    [Z0{i}, ~] = InitializeSIGs(X{i}, pn, 0);
end
Z0_initial = Z0;

%% initialize U, F and w
U = zeros(numC);                      % FIX % change matrix size to ni x ni => n x n
for i = 1:m
    U = U + G{i}'*Z0{i}*G{i};         % FIX % change Z0 to nxn by ( GT *  Z * G )
end
U = U/m;

for j = 1:numC      
    U(j,:) = U(j,:)/sum(U(j,:));
end

sU = (U+U')/2;
D = diag(sum(sU));
L = D - sU;
[F, ~, evs]=eig1(L, c, 0);

w = ones(1,m)/m;

idxx = cell(1,m);
ed = cell(1,m);

% ******
for v = 1:m
    ed{v} = L2_distance_1(X{v}, X{v});
    [~, idxx{v}] = sort(ed{v}, 2); % sort each row
end

%%  update ...
for iter = 1:NITER
    fprintf('Loop:%d \n',iter);
    
    % update Z^v
    parfor v = 1:m
        Z0{v} = zeros(num);
        for i = 1:num
            % TEMP For U
            U_star = G{v}*U*G{v}';
            
            id = idxx{v}(i,2:pn+2);
            di = ed{v}(i, id);
            numerator = di(pn+1)-di+2*w(v)*U_star(i,id(:))-2*w(v)*U_star(i,id(pn+1)); % FIX % U => G * U * GT
            denominator1 = pn*di(pn+1)-sum(di(1:pn));
            denominator2 = 2*w(v)*sum(U_star(i,id(1:pn)))-2*pn*w(v)*U_star(i,id(pn+1)); % FIX % U => G * U * GT
            Z0{v}(i,id) = max(numerator/(denominator1+denominator2+eps),0);
        end
    end
    
    % update w
    parfor v = 1:m
        US = U - G{v}'*Z0{v}*G{v}; % FIX % % change Z0 to nxn by ( GT *  Z * G )
        distUS = norm(US, 'fro')^2;
        if distUS == 0
            distUS = eps;
        end
        w(v) = 0.5/sqrt(distUS);
    end
    % disp(['weights: ',num2str(w)]);
    
    
    % update U
    dist = L2_distance_1(F',F');
    U = zeros(numC); % FIX % change matrix size to ni x ni => n x n
    for i=1:num
        idx = zeros();
        for v = 1:m
            s0 = Z0{v}(i,:);
            idx = [idx,find(s0>0)];
        end
        idxs = unique(idx(2:end));
        if islocal == 1
            idxs0 = idxs;
        else
            idxs0 = 1:num;
        end
        for v = 1:m
            s1 = Z0{v}(i,:);
            si = s1(idxs0);
            di = dist(i,idxs0);
            mw = m*w(v);
            lmw = lambda/mw;
            q(v,:) = si-0.5*lmw*di;
        end
        U(i,idxs0) = UpdateUStar(q,m);
        clear q;
    end

    sU = U;
    sU = (sU+sU')/2;
    D = diag(sum(sU));
    L = D-sU;
    F_old = F;
    [F, ~, ev]=eig1(L, c, 0, 0);
    evs(:,iter+1) = ev;
    % update lambda and the stopping criterion
    fn1 = sum(ev(1:c));
    fn2 = sum(ev(1:c+1));
    %fprintf('fn1=%d \n fn2=%d',fn1,fn2);
    if fn1 > zr
        lambda = 2*lambda;
        %lambda = 1.5*lambda;
    elseif fn2 < zr
        lambda = lambda/2;
        %lambda = lambda/1.5;
        F = F_old;
    else
        disp(['iter = ',num2str(iter),' lambda:',num2str(lambda)]);
        break;
    end
end

end





