

function [HF,WF]=test_alg(X,H,Z,G,invH,num_view,max_iter)

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
           %coeW = getCoeW(X{iv},H{iv},W{iv});
           %W{iv} = W{iv}.*coeW;
    end
    

    
    % ---------- Update H ----------%
    for iv = 1:num_view
        coeH = getCoeH(X{iv},WF{iv},H{iv},Z{iv});
        % coeH1 = getCoeH1(X{iv},WF{iv},H{iv},Z{iv});
        
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

%     for iv = 1:num_view
%         dim = size(HF{iv});
%         I = eye(dim(2));
%         to_inv = HF{iv}'*HF{iv}+gamma*I;
%         invV = inv(to_inv);
%         Z{iv} = invV*(HF{iv}'*HF{iv});
%     end
    
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
% [F,S]= graphfusion(Z,HF,G,truthF,lambda,beta,gamma);

% GMIC graph fusion


% complete H
 for iv = 1:num_view
     HF{iv} = HF{iv}*G{iv};
 end

end
