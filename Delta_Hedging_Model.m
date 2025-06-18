clear; clc; 

%% Set up display format for currency
format bank

%% Define portfolio and market parameters
% Strike prices for put and call options
Kp = 95; 
Kc = 110; 
% Volatility for Stock A and Index, correlation, and risk-free rate
s_A = 0.122; 
s_I = 0.18; 
rho = 0.911; 
r = 0.02; 
% Expected returns for Stock A and Index
returns = [0.083, 0.08]; 

%% Create a table summarizing the portfolio's asset classes
AssetClasses = {'Stock A (long)', 'Index Put (long, ATM)', 'Index Call (short, OTM)'};
Weights = [1, 1, -1]; % 1 for long, -1 for short
ExpReturn = [0.083, NaN, NaN]; % Only stock has explicit return
Volatility = [0.122, NaN, NaN]; % Only stock has explicit volatility
PortfolioTable = table(AssetClasses', Weights', ExpReturn', Volatility', ...
    'VariableNames', {'Asset_Class', 'Weight', 'Annualized_Return', 'Annualized_Volatility'})

%% Set up simulation parameters and time grid
C = [1 rho; rho 1]; % Correlation matrix
S = [s_A^2, s_A*s_I*rho; s_A*s_I*rho, s_I^2]; % Covariance matrix
V = chol(S); % Cholesky decomposition for correlated random numbers
Week = [0:1:52]'; % Weeks in a year
Time = (repmat(52,length(Week),1)-Week)/52; % Time to maturity in years
Time(end) = (Week(end)+0.5-52)/52; % Adjust last time step
N = 1000; % Number of Monte Carlo simulations

%% Initialize storage for simulation results
DynamicHedge = zeros(N, 6); % Stores results for each simulation
costHedge = zeros(N, 1);    % Stores hedging costs

for n = 1:N
    %% Simulate price paths for Stock A and Index
    M = [100 100; zeros(52,2)]; % Initial prices for Stock A and Index
    for t = 2:length(Time)
        % Simulate next price using geometric Brownian motion
        M(t,:) = M(t-1,:).*exp((returns-0.5*diag(S)')*1/52 + randn(1,2)*V*sqrt(1/52));
    end

    %% Calculate Black-Scholes option parameters and prices
    d1c = (log(M(:,2)/Kc) + (r + 0.5*s_I^2).*Time)./(s_I*sqrt(Time));
    d2c = d1c-s_I*sqrt(Time); 
    d1p = (log(M(:,2)/Kp) + (r + 0.5*s_I^2).*Time)./(s_I*sqrt(Time));
    d2p = d1p-s_I*sqrt(Time); 
    put = repmat(Kp,length(Time),1).*exp(-r.*Time).*normcdf(-d2p) - M(:,2).*normcdf(-d1p); 
    call = M(:,2).*normcdf(d1c) - repmat(Kc,length(Time),1).*exp(-r.*Time).*normcdf(d2c);
    delta = normcdf(d1c); % Delta of the call option

    %% Track changes in delta and associated hedging transactions
    X = [delta(1),M(1,2)*delta(1),r*M(1,2)*delta(1)/52;zeros(52,3)]; 
    for t = 2:length(Time) 
        X(t,1) = delta(t) - delta(t-1);  % Change in delta
        X(t,2) = M(t,2) * X(t,1);        % Value of shares bought/sold
        X(t,3) = r * sum(X(1:t-1,2)) / 52; % Interest cost on hedged position
    end

    %% Create a table summarizing the hedging process for each week
    DeltaHedgeAcct = table(Week,Time,M(:,2),M(:,1),d1c,d2c,d1p,d2p,put,call,delta,X(:,1),X(:,2),X(:,3),...
        'VariableNames',{'Week','Time','Index','StockA','d1c','d2c','d1p','d2p','Put','Call','Delta','Shares','Cost','Interest'}); 

    %% Display selected rows from the hedging table for review
    if n == 1
        FirstRows = DeltaHedgeAcct(1:10,:);
        LastRows = DeltaHedgeAcct(45:52,:);
        disp('First 10 rows for review:');
        disp(FirstRows);
        disp('Rows 45-52 for review:');
        disp(LastRows);
    end

    %% Calculate final values for put, call, and portfolio performance
    putValue = max(Kp-M(end,2),0); 
    callValue = max(M(end,2)-Kc,0);
    callValue_NH = call(1); % Call value at start (non-hedged)
    PnL = (putValue - put(1) + callValue - call(1) - sum(X(:,2)) - sum(X(:,3))) / 100; % Net profit/loss with hedging
    PnL_NH = (putValue - put(1) + callValue_NH - call(1)) / 100; % Net profit/loss without hedging
    BuyHold = (M(end,1)-M(1,1))/100; % Buy and hold return

    %% Store results for this simulation
    DynamicHedge(n,:) = [putValue, callValue, PnL, BuyHold, PnL_NH, callValue_NH]; 
    costHedge(n) = (M(end,2)>Kc)*Kc - sum(X(:,2)); 
end

%% Plot histogram of net profit/loss from all simulations
figure(1);
histogram(DynamicHedge(:,3), 21);
title('Net P/L Distribution (N=1000)');
xlabel('P/L (% of initial $100)');
ylabel('Frequency');

%% Calculate and display mean and standard deviation of net profit/loss
meanPnL = mean(DynamicHedge(:,3));
stdPnL = std(DynamicHedge(:,3));
fprintf('Net P/L Statistics:\n');
fprintf('Mean Net P/L: %.4f\n', meanPnL);
fprintf('Standard Deviation of Net P/L: %.4f\n', stdPnL);

%% Plot histograms for put value, call value, P&L, and buy & hold return
figure(2);
subplot(2,2,1); 
histogram(DynamicHedge(:,1), 21); 
title('Put Value');
xlabel('Value ($)');
subplot(2,2,2); 
histogram(DynamicHedge(:,2), 21); 
title('Call Value');
xlabel('Value ($)');
subplot(2,2,3); 
histogram(DynamicHedge(:,3), 21); 
title('P&L');
xlabel('P/L (%)');
subplot(2,2,4); 
histogram(DynamicHedge(:,4), 21); 
title('Buy & Hold');
xlabel('Return (%)');
sgtitle('Histogram Comparisons');

%% Calculate and display percentiles for key portfolio metrics
Centiles = [1, 25, 50, 75, 99];
CentileValues = prctile(DynamicHedge(:,1:4), Centiles);
PercentilesTable = table(Centiles', CentileValues(:,1), CentileValues(:,2), CentileValues(:,3), CentileValues(:,4),...
    'VariableNames', {'Percentile', 'PutValue', 'CallValue', 'PnL', 'BuyAndHold'});
disp(PercentilesTable);

%% Plot additional histograms for non-hedged results
figure(3);
subplot(2,1,1); 
histogram(DynamicHedge(:,5), 21);
title('P&L without Hedging');
xlabel('P/L (%)');
subplot(2,1,2); 
histogram(DynamicHedge(:,6), 21); 
title('Initial Call Value (Non-Hedged)');
xlabel('Value ($)');
sgtitle('Supplementary Analysis')
