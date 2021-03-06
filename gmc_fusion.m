
function [y, U, Z0, Z0_initial, F, evs] = gmc_fusion(X, c,G, lambda, normData) % gamma

%% input:
% X{}: multi-view dataset, each cell is a view, each column is a data point
% c: cluster number
% lambda: parameter (default 1)
%% output:
% S0: similarity-induced graph (SIG) matrix for each view
% y: the final clustering result, i.e., cluster indicator vector
% U: the learned unified matrix
% F: the embedding representation
% evs: eigenvalues of learned graph Laplacian matrix

NITER = 20;
zr = 10e-11;
pn = 15; % number of neighbours for constructS_PNG
islocal = 1; % only update the similarities of neighbors if islocal=1
if nargin < 3+1
    lambda = 1;
end;
if nargin < 4+1
    normData = 1;
end;

num = size(X{1},2); % number of instances
numC = size(G{1},2);
m = length(X); % number of views
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
U = zeros(numC);         % FIX % change matrix size to ni x ni => n x n
for i = 1:m
    U = U + G{i}'*Z0{i}*G{i};      % FIX % change Z0 to nxn by ( GT *  Z * G )
end
U = U/m;

for j = 1:numC % FIX
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
    %ed{v} = L2_distance_1(H{v},H{v});
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
%         for j = 1:num
%             normItem = sum(S0{v}(j,:));
%             if normItem == 0
%                 normItem = eps;
%             end;
%             S0{v}(j,:) = S0{v}(j,:)/normItem;
%         end;
    end
    
%     for i = 1:200
%     
%     % update W,
%     for iv = 1:m
%         %U{iv}=X{iv}*V{iv}'*invVV{iv};
%         Xtmp   = X{iv};
%         Htmp   = H{iv};
%         invtmp = invH{iv};
%         Xtmp    = gpuArray(Xtmp);
%         Htmp    = gpuArray(Htmp);
%         invtmp  = gpuArray(invtmp);
%         Wtmp= Xtmp * Htmp'*invtmp;
%         Wtmp   = gather(Wtmp);
%         WF{iv}  = Wtmp;
%            %coeW = getCoeW(X{iv},H{iv},W{iv});
%            %W{iv} = W{iv}.*coeW;
%     end
%     
%     
%     
%     % update H,
%     
%     for iv = 1:m
%         coeH = getCoeH(X{iv},WF{iv},H{iv},Z0{iv});
%         % coeH1 = getCoeH1(X{iv},WF{iv},H{iv},Z{iv});
%         
%         H{iv} = H{iv}.*coeH;
%         
%         % ---------- get inv of H again ---------
%         Htemp = H{iv};
%         Htemp = gpuArray(Htemp);
%         invtemp = inv(Htemp*Htemp');
%         invtemp =gather(invtemp);
%         invH{iv} = invtemp;
%     end
%     HF = H;
%     
%     objhistory = CalculateObj(X{1}, WF{1}, HF{1}');  
%     disp(objhistory);
%     
%     
%     end
    
    
    
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
        U(i,idxs0) = SloutionToP19(q,m);
        clear q;
    end
%         % choose the top-k neighbors
%         [~, ids] = sort(U,2,'descend');
%         ts = zeros(num);
%         for i =1:num
%             ts(i,ids(i,1:pn)) = U(i,ids(i,1:pn));
%         end
%         for j = 1:num
%             ts(j,:) = ts(j,:)/sum(ts(j,:));
%         end
%         sU = ts;
    % update F
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
    fprintf('fn1=%d \n fn2=%d',fn1,fn2);
    if fn1 > zr
        %lambda = 2*lambda;
        lambda = 1.5*lambda;
    elseif fn2 < zr
        %lambda = lambda/2;
        lambda = lambda/1.5;
        F = F_old;
    else
        disp(['iter = ',num2str(iter),' lambda:',num2str(lambda)]);
        break;
    end
end

%% generating the clustering result
[clusternum, y]=graphconncomp(sparse(sU)); y = y';
if clusternum ~= c
    fprintf('Can not find the correct cluster number: %d\n', c)
end; 


end


function coeH = getCoeH(X,W,H,Z)
% numerator = np.dot(W.T,X) + np.dot(H,Z) + np.dot(H,Z.T);
% denominator = np.dot( np.dot(W.T,W), H ) + H + np.dot( np.dot(H,Z),Z.T);
% return numerator/denominator; # (k,n) matrix

% W = ( m x 5)
% X = ( m x n)
% H = ( 5 x n)
% Z = ( n x n)

numerator = W'*X +(H*Z)+(H*Z');
denominator = (W'*W*H)+H+((H*Z)*Z');
coeH = numerator./denominator;
end




