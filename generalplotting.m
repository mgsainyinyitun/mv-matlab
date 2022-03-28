% x-> missing percent > y=> accuracy 
% database -> bbcsport4vbigRnSp

X = [10 20 30 40 50];
% y1 = [80.1724 79.3103 81.0345 61.2069 41.8966];
% y2 = [83.1000 75.7875 63.6000 52.5875 37.2750];
% y3 = [69.4500,72.3000 61.7000 55.3500 41.8500];
% y4 = [82.6000,82.4700,83.1500,84.7500 46.0500];
% y5 = [67.6000 66.1500 58 49.2000 42.7];
% y6 = [74.2741 75.6581 72.6323 69.9457 37.0421];

% percent vs nmi 

y1 = [72.6406 75.1442 72.6566 58.0074 24.0983]; % bbcsport4vbigRnSp
y2 = []; % 100Leaves
y3 = []; % ORL
y4 = []; % mfeatRnSp
y5 = []; % caltech7


plot(X,y1,'-o',...
    X,y2,'-o',...
    X,y3,'-o',...
    X,y4,'-o',...
    X,y5,'-o',...
    'MarkerFaceColor',[0,0,1]);


title('Accuracy over missing rate of all dataset');
ylim([0 100])
ylabel('Accuracy (%)');
xticks([0 10 20 30 40 50])
xlabel('missing rate (%)');

legend('bbcsport4vbigRnSp','100Leaves','ORL',...
'mfeatRnSp','caltech7',...
'Location','southwest');
grid on
grid minor








