clear, clc;
small_data_1=importdata('small_exp1.txt');

T = small_data_1.data(:,2); % Temperature
m = small_data_1.data(:,3); % Dependent variable, m, mass
t = small_data_1.data(:,1); % Time

%% CALCULATE REGRESSION PARAMETERS
% Step 1: Normalize variables****************** 
% This step minimizes numerical problems when fitting polynomials. 
% In general, it is a good idea to normalize variables, especially in the 
% case of nonlinear regression, it helps convergence.
Tn = (T-mean(T))./std(T);  % Normalized temperature (-)
mn = (m-mean(m))./std(m);  % Normalized mass (-)

% Step 2: Calculate regression parameters********
n = 61;          % number of data points
p = 59;          % number of regression parameters

v=n-p;           % Degrees of freedom

X = [ones(length(small_data_1.data(:,1)),1),small_data_1.data(:,1)];          % Define matrix X

beta = inv(X'*X)*X'*T;       % Calculate regression parameters

%% DETERMINE IF AT LEAST ONE REGRESSION PARAMETER IS STATISTICALLY SIGNIFICANT
% Hint: See MATLAB program in assignment 27
SSE = (T-X*beta)'*(T-X*beta); % See Appendix
MSE = SSE/v;   % Mean square error
T_bar = mean(T);  % Average 
SSr = sum((T-Tn).^2);           %Is this right?
MSr = SSr/(p-1);
Fobs = MSr/MSE;
alpha = 0.05;
Ftab = finv(1-alpha,p-1,n-p);

%% CALCULATE THE CORRELATION MATRIX
M = (X'*X^-1)/(sqrt((X'*X)^-1*(X'*X)^-1));    
Cm = zeros(p,p); %Preallocate
for i=1:p    
    for j=1:p    
        Cm(i,j)= M(i,j)./sqrt(M(i,i)*M(j,j)); % See definition of correlation matrix    
    end
end


