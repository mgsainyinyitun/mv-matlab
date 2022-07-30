load('./database/bbcsport4vbigRnSp.mat');

X = data{1}';
num = size(X);
num = num(1);

for  j = 1:num
  normItem = std(X(:,j));
  if (0 == normItem)
       normItem = eps;
  end
  X(:,j) = (X(:,j)-mean(X(:,j)))/(normItem);
end
Y = tsne(X,'Algorithm','barneshut','NumPCAComponents',50);

num = size(Y);
num = num(2);
for  j = 1:num
  normItem = std(Y(:,j));
  if (0 == normItem)
       normItem = eps;
  end
  Y(:,j) = (Y(:,j)-mean(Y(:,j)))/(normItem);
end

figure
numGroups = length(unique(truth));
clr = hsv(numGroups);
gscatter(Y(:,1),Y(:,2),truth,clr)