% Clear Everything
clear;
clc;

% Data Retrieval
t = readtable('stockData.csv');
data = table2array(t(:,2:end));
names = t.Properties.VariableNames;
dates = table2array(t(:,1));

% Data Adjustment
invalidNames = {'Dow_DuPont','Visa'};
for i = 1:length(invalidNames)
    invalidIndices(i) = strmatch(cellstr(invalidNames(i)), names, 'exact');
end
data(:,invalidIndices-1) = [];
names([1,invalidIndices]) = [];

% Significant Dates
preStart = datenum(2007,1,3);
peakStart = datenum(2008,9,2);
postStart = datenum(2009,6,1);
preStartIndex = datefind(preStart,dates);
peakStartIndex = datefind(peakStart,dates);
postStartIndex = datefind(postStart,dates);

% Data Size
[dateCount,stockCount] = size(data);

%%% Question 1: Correlation Analysis
logData = log(data);
logReturn = logData(2:end,:) - logData(1:end-1,:);
R = corrcoef(logReturn);
R = R - 0.5*eye(stockCount);

% Maximum
[colMax, maxX] = max(R);
[wholeMax, maxY] = max(colMax);
%disp([wholeMax, names(maxX(maxY)), names(maxY)]);
% Minimum
[colMin, minX] = min(R);
[wholeMin, minY] = min(colMin);
% disp([wholeMin, names(minX(minY)), names(minY)]);

% Split Correlation Analysis
logReturnBefore = logData(2:preStartIndex-1,:) - logData(1:preStartIndex-2,:);
logReturnPre = logData(preStartIndex+1:peakStartIndex-1,:) - logData(preStartIndex:peakStartIndex-2,:);
logReturnPeak = logData(peakStartIndex+1:postStartIndex-1,:) - logData(peakStartIndex:postStartIndex-2,:);
logReturnPost = logData(postStartIndex+1:end,:) - logData(postStartIndex:end-1,:);
rBefore = corrcoef(logReturnBefore);
rPre = corrcoef(logReturnPre);
rPeak = corrcoef(logReturnPeak);
rPost = corrcoef(logReturnPost);

% Plot Heat Map
% imagesc(rPost);
% colorbar;
% xlabel('Stocks')
% ylabel('Stocks')

% Plot Histogram
pairCount = [];
for i = 1:stockCount
    for j = 1:i
        pairCount(length(pairCount)+1) = rPost(i,j);
    end
end
% histogram(pairCount);
% xlabel('Correlation Coefficient');
% ylabel('Total Pairs of Stocks');

%%% Question 2: Portfolio Short-Selling
% Initials
gfcData = data(preStartIndex:end,:);
gfcReturn = (gfcData(2:end,:) - gfcData(1:end-1,:))./gfcData(1:end-1,:);
S = cov(gfcReturn);
e = ones(28,1);
r = mean(gfcReturn)';

% Computing Critical Line
a = e'*inv(S)*e;
b = r'*inv(S)*e;
c = r'*inv(S)*r;
d = a*c-b^2;
alpha = (1/a)*inv(S)*e;
beta = inv(S) * (r - (b/a)*e);

% Displaying What t-values short sell each stock
alwaysShort = [];
neverShort = [];
for i = 1:stockCount
    indif = -alpha(i)/beta(i);
    if beta(i) > 0
        fullstr = strcat(names(i),' is short-sold if t < ',num2str(indif));
        if indif < 0
            neverShort(length(neverShort)+1) = i;
        end
    else
        fullstr = strcat(names(i),' is short-sold if t > ',num2str(indif));
        if indif < 0
            alwaysShort(length(alwaysShort)+1) = i;
        end
    end
    % Display When each stock is short sold
    %disp(fullstr);
end
% Display what stocks are always/never short sold
%disp(names(alwaysShort));
%disp(names(neverShort));

%%% Question 3: Portfolio Analysis
t = 0.15;
x = alpha + t*beta;
invest = 200000*x;
% Display Dollar investments
%disp(invest);

% Risk and Return
exp = (b+d*t)/a;
dev = sqrt((1+d*t^2)/a);
% disp(exp);
% disp(dev);

% hold on;
% 
% % MVF and EF
% t_plot = linspace(-0.35,0.35,1000);
% exp = (b+d*t_plot)/a;
% dev = sqrt((1+d*t_plot.^2)/a);
% plot(dev,exp,'DisplayName','Minimum Variance Frontier');
% plot(dev(end/2:end),exp(end/2:end),'DisplayName','Efficient Frontier');
% 
% % Plotting Random Portfolios
% portfolioCount = 0;
% sigmas = [];
% mus = [];
% while portfolioCount <= 1000
%     x_rand = 40*rand(28,1)-20;
%     x_sum = sum(x_rand);
%     x_rand_norm = x_rand/sum(x_rand);
%     mu_rand = x_rand_norm'*r;
%     sigma_rand = sqrt(x_rand_norm'*S*x_rand_norm);
%     if sigma_rand <= 0.05 && max(x_rand_norm) <= 20 && min(x_rand_norm) >= -20
%         portfolioCount = portfolioCount + 1;
%         sigmas(length(sigmas)+1) = sigma_rand;
%         mus(length(mus)+1) = mu_rand;
%     end
% end
% plot(sigmas,mus,'*','DisplayName','Random Portfolios');
% 
% % Plotting Stocks
% scatter(std(gfcReturn),mean(gfcReturn),'DisplayName','Individual Stocks');
% xlabel('Standard Deviation');
% ylabel('Expected Return');
% 
% % Indifference Curve
% z = -t*(r'*x) + 0.5*x'*S*x;
% sigma = linspace(0,0.045,1000);
% mu = (0.5*sigma.^2 - z)/t;
% plot(sigma,mu,'DisplayName','Indifference Curve');
% 
% legend('show');

%%% Question 8: EOFs

% Plotting Ordered Spectrum
%hold on;
A_full = ((data(2:end,:) - data(1:end-1,:))./data(1:end-1,:))';
A_before = ((data(2:preStartIndex-1,:) - data(1:preStartIndex-2,:))./data(1:preStartIndex-2,:))';
A_pre = ((data(preStartIndex+1:peakStartIndex-1,:) - data(preStartIndex:peakStartIndex-2,:))./data(preStartIndex:peakStartIndex-2,:))';
A_peak = ((data(peakStartIndex+1:postStartIndex-1,:) - data(peakStartIndex:postStartIndex-2,:))./data(peakStartIndex:postStartIndex-2,:))';
A_post = ((data(postStartIndex+1:end,:) - data(postStartIndex:end-1,:))./data(postStartIndex:end-1,:))';
% Finding EigenValues
eFull = eig(A_full*A_full');
eBefore = eig(A_before*A_before');
ePre = eig(A_pre*A_pre');
ePeak = eig(A_peak*A_peak');
ePost = eig(A_post*A_post');
% Plot Ordered Spectrum of Eigen Values
% plot(eFull,'DisplayName','All Time');
% plot(eBefore,'DisplayName','Before GFC');
% plot(ePre,'DisplayName','Pre GFC');
% plot(ePeak,'DisplayName','Peak GFC');
% plot(ePost,'DisplayName','Post GFC');
% ylabel('EigenValues');
% legend('show');

% Getting EigenVectors
[eVecFull,eD1] = eig(A_full*A_full');
[eVecBefore,eD2] = eig(A_before*A_before');
[eVecPre,eD3] = eig(A_pre*A_pre');
[eVecPeak,eD4] = eig(A_peak*A_peak');
[eVecPost,eD5] = eig(A_post*A_post');
% Plotting EigenVectors
% plot(eVecFull(:,1),'DisplayName','All Time 1');
% plot(eVecFull(:,2),'DisplayName','All Time 2');
% plot(eVecBefore(:,1),'DisplayName','Before GFC 1');
% plot(eVecBefore(:,2),'DisplayName','Before GFC 2');
% plot(eVecPre(:,1),'DisplayName','Pre GFC 1');
% plot(eVecPre(:,2),'DisplayName','Pre GFC 2');
% plot(eVecPeak(:,1),'DisplayName','Peak GFC 1');
% plot(eVecPeak(:,2),'DisplayName','Peak GFC 2');
% plot(eVecPost(:,1),'DisplayName','Post GFC 1');
% plot(eVecPost(:,2),'DisplayName','Post GFC 2');
% xlabel('Component in Eigenvector');
% ylabel('Component Value');
% legend('show');


