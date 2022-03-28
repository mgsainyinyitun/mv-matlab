% x-> missing percent > y=> accuracy 


X = [10 20 30 40 50];

% percent vs nmi 

y1 = [63.1213 66.3459 63.5873 56.0248 33.6355]; % OUR
y2 = [38.6383 45.3094 40.3537 39.767 28.1246]; % PMVC
y3 = [40.284362 38.123826 40.775 40.8318 3.4264]; % IMG
y4 = [36.4858 32.4897 23.7778 29.8042 31.6461]; % MIC
y5 = [46 45.3338 45.003 45.6132 32.1068]; % IMSC-AGL
y6 = [43.5064 40.818 41.6996 35.6921 22.1684] % DAIMC
y7 = [20.34 12.4263 12.134 11.4488 14.4296] ; % BSV


plot(X,y1,'-o',...
    X,y2,'-o',...
    X,y3,'-o',...
    X,y4,'-o',...
    X,y5,'-o',...
    X,y6,'-o',...
    X,y7,'-o',...
    'MarkerFaceColor',[1,1,1]);


title('Comparison of NMI for BBC sport');
ylim([0 100])
ylabel('NMI (%)');
xticks([0 10 20 30 40 50])
xlabel('missing rate (%)');

legend('OUR','PMVC','IMG',...
'MIC','IMSC-AGL','DAIMC','BSV',...
'Location','northeast');
grid on
grid minor








