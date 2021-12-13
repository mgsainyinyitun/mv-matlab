function [F,S]=graphfusion(Zcomplete,HF,G,truth,beta,gamma)
% truth is the true class label.
% mv=size(X,1);
num_view=size(Zcomplete,2);
n=size( G{1}'*Zcomplete{1}*G{1} ,1);

Z=eye(n);

Zv=repmat(Z,[1,1,num_view]);
c=length(unique(truth));

wv=ones(num_view,1)/num_view;
%options = optimset( 'Algorithm','interior-point-convex','Display','off');
for ii=1:200
    fprintf('Number of Iteration:%i .\n',ii);
     % for F
    % use Z to achieve S(graph fusion) -> Z (n x n) 
    
    Z(Z<0)=0;      % Z less than 0 % fusion
    Z= (Z+Z')/2;         % Z+Z' /2   fusion
    Zold=Z;              % fusion
    S = Z;
    D = diag(sum(Z));    % fusion
    L = D-Z;
    
    [F, temp, ev]=eig1(L, c, 0); % solve eigen value decom
    
    for i=1:num_view
        % f=Z{i};
        % update Z
        % Zv(:,:,i)=(f*f'+alpha*eye(n)+beta*wv(i)*eye(n))\(beta*wv(i)*Z+f*f');
        
        %------------
         dim = size(HF{i});
         I = eye(dim(2));
         to_inv = HF{i}'*HF{i}+gamma*I;
         invV = inv(to_inv);
         Zcomplete{i} = invV*(HF{i}'*HF{i});
  
        %-----------
        % Z{iv} = G{iv}'*Z{iv}*G{iv};
        Zv(:,:,i)=G{i}'*Zcomplete{i}*G{i};
        
        T=Zv(:,:,i);
        T(T<0)=0;
        T=(T+T')/2;
        Zv(:,:,i)=T;
        wv(i)=1/2/norm(Zv(:,:,i)-Z,'fro');
        M(:,:,i)=wv(i)*Zv(:,:,i);
    end
    
    % for S
    parfor ij=1:n % 1 - 1200
        all=distance(F,n,ij);  % n x 1  
        Z(:,ij)=(sum(M(:,ij,:),3)   -gamma*all'/(4*beta))/    sum(wv);
    end
    
    % Z -> S (fusion)
    % ii greater 5 and 
    error = norm(Z-Zold,'fro')/norm(Zold,'fro');
    fprintf('Error value is %d \n',error);
    disp('Error is:');
    disp(error);
    if ii>5 &( (norm(Z-Zold,'fro')/norm(Zold,'fro') )<1e-3)
        break
    end
     
end

% res=zeros(10,3);
% for ij=1:10
% %actual_ids= kmeans(F, c, 'emptyaction', 'singleton', 'replicates', 1, 'display', 'off');
% actual_ids= kmeans(F, c);
% [res(ij,:)] = ClusteringMeasure( actual_ids,s);
% end
% result(1,1)=mean(res(:,1));result(1,2)=std(res(:,1));
% result(2,1)=mean(res(:,2));result(2,2)=std(res(:,2));
% result(3,1)=mean(res(:,3));result(3,2)=std(res(:,3));

end


function [all]=distance(F,n,ij);
  for ji=1:n
            all(ji)=(norm(F(ij,:)-F(ji,:)))^2;
  end
end










