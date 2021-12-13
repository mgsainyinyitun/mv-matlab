function [eigvec, eigval, eigval_full] = eig1(A, c, isMax, isSym) % L, c, 0,

if nargin < 2 % if argument is 1 or 0
    c = size(A,1);
    isMax = 1;
    isSym = 1;
elseif c > size(A,1) % if c > size(A,1);
    c = size(A,1);
end;

if nargin < 3 % arg => 2 or 1 
    isMax = 1;
    isSym = 1;
end;

if nargin < 4 % arg => 3 or 2 or 1
    isSym = 1;
end;

if isSym == 1
    A = max(A,A');
end;


[v d] = eig(A); % (nxn) , (n)
% V => vector ( full size matrix )  , D = value (diagonal vecor_ 

d = diag(d);
%d = real(d);

if isMax == 0
    [d1, idx] = sort(d); % sort eigen value ascending and get index
else
    [d1, idx] = sort(d,'descend'); % sort eigen value descending and get index
end;

idx1 = idx(1:c); % first c item 1,2,3,4,5,6

eigval = d(idx1);
eigvec = v(:,idx1);

eigval_full = d(idx);




