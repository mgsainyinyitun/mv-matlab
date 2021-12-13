function [HF,WF,Z,F,S]=algorithm(X,H,Z,G,invH,num_view,truthF,gamma,max_iter)
% X = list of all view {x1,x2,x3,...}
% W = list of all w {w1,w2,w3,...}
% H = list of all h {h1,h2,h3,...}

for iter = 1:max_iter
    % ---------- Update W ----------%
    for iv = 1:num_view
        %U{iv}=X{iv}*V{iv}'*invVV{iv};
        Xtmp   = X{iv};
        Htmp   = H{iv};
        invtmp = invH{iv};
        Xtmp    = gpuArray(Xtmp);
        Htmp    = gpuArray(Htmp);
        invtmp  = gpuArray(invtmp);
        Wtmp= Xtmp * Htmp'*invtmp;
        Wtmp   = gather(Wtmp);
        WF{iv}  = Wtmp;
    end
    % ---------- Update H ----------%
    for iv = 1:num_view
        coeH = getCoeH(X{iv},WF{iv},H{iv},Z{iv});
        H{iv} = H{iv}.*coeH;
        % ---------- get inv of H again ---------
        Htemp = H{iv};
        Htemp = gpuArray(Htemp);
        invtemp = inv(Htemp*Htemp');
        invtemp =gather(invtemp);
        invH{iv} = invtemp;
    end
    
    HF = H;
    
    
    % --------- update Z ----------------

    for iv = 1:num_view
        dim = size(HF{iv});
        I = eye(dim(2));
        to_inv = HF{iv}'*HF{iv}+gamma*I;
        invV = inv(to_inv);
        Z{iv} = invV*(HF{iv}'*HF{iv});
    end
    
%     for iv = 1:num_view
%         e = getError(X{iv},WF{iv},HF{iv});
%         %fprintf('Error of %i view is %i \n',iv,e);
%     end
    
    % objhistory = CalculateObj(X{1}, WF{1}, HF{1}');  
    % disp(objhistory);
    
end
% complete graph
% for iv=1:num_view
%     Z{iv} = G{iv}'*Z{iv}*G{iv};
% end
% graph fusion 
[F,S]= graphfusion(Z,HF,G,truthF,100,0.01);
    
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

function coeZ = getCoeZ(H,Z)
%      numerator = np.dot(H.T,H);
%      denominator = np.dot( np.dot(H.T,H),Z) ;
%      return numerator/denominator;  # (n,n) matrix
numerator = H'*H;
denominator = H'*H*Z;
coeZ = numerator./denominator;
end



function e = getError(X,W,H)
normX = norm(X);
H = H';
%sum_squared = norm_X * norm_X - 2 * np.trace(H.T.dot(X.T.dot(W)))+ np.trace((W.T.dot(W)).dot(H.T.dot(H)))
e = normX*normX - 2*trace(H'*(X'*W)) + trace(W'*W*(H'*H));
e = sqrt(max(e));
% disp(e);
end






    