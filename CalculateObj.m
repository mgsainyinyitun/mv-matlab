function [obj, dV] = CalculateObj(X, U, V, deltaVU, dVordU)

    if ~exist('deltaVU','var')
        deltaVU = 0;
    end
    if ~exist('dVordU','var')
        dVordU = 1;
    end
    
    
    dV = []; % empty matrix 
    
    maxM = 62500000;
    
    [mFea, nSmp] = size(X); % feature x sample (row, columns) (10 x 10) => [10,10]
    mn = numel(X);          % number of element in X => 100
    nBlock = floor(mn*3/maxM); % appro 0 % (100*3)/625000000;
    %e.g floor(3.5) => 3 

    if mn < maxM  % number of element less than maxM <  => (true) 100 < 62500000
        
        dX = U*V'-X;  % find different between dx = all zero matrix 
        
        obj_NMF = sum(sum(dX.^2)); % single maximun value
        
        if deltaVU % deltaVU = 0 % skip
            if dVordU % if dVordU == 1 
                dV = dX'*U;
            else
                dV = dX*V;
            end
        end
        % end of deltVU = 0
    else
        obj_NMF = 0; % reset to zero 
        if deltaVU  % 0
            if dVordU
                dV = zeros(size(V));
            else
                dV = zeros(size(U));
            end
        end
        % =================
        for i = 1:ceil(nSmp/nBlock) % 10
            if i == ceil(nSmp/nBlock)
                smpIdx = (i-1)*nBlock+1:nSmp;
            else
                smpIdx = (i-1)*nBlock+1:i*nBlock;
            end
            dX = U*V(smpIdx,:)'-X(:,smpIdx);
            obj_NMF = obj_NMF + sum(sum(dX.^2));
            if deltaVU
                if dVordU
                    dV(smpIdx,:) = dX'*U;
                else
                    dV = dU+dX*V(smpIdx,:);
                end
            end
        end
        if deltaVU
            if dVordU
                dV = dV ;
            end
        end
    end
   %obj_Lap = alpha*sum(sum((L*V).*V));
   
    obj = obj_NMF;