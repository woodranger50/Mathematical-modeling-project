clear all, clear, clc;
format long
t_span = 0:45:1800;
m_water = 0.1795;              %[kg]   Initial mass
T_water = 80+273.15;        %[K]   Initial temperature
y0=[T_water,m_water];


[t,dt]=ode45(@Mass_heat_flux,t_span,y0);

T=dt(:,1);
T_water=T-273.15;
m=dt(:,2);
m_water=m*10^3;

% figure
% hold on
% subplot(2,1,1)
% title('Temperature'); xlabel('Time (seconds)'); ylabel('Temperature [°C]')
% plot(t,T_water,'-o')
% subplot(2,1,2)
% title('Mass'); xlabel('Time (seconds)'); ylabel('Mass [g]')
% plot(t,m_water,'-o')
% hold off

%% Physical data plot
addpath('C:\Chalmers\Year 3\Matematisk modellering\Github\Mathematical-modeling-project');

small_data_1=importdata('small_exp1.txt');
small_data_2=importdata('small_exp2.txt');
small_data_3=importdata('small_exp3.txt');
small_data_4=importdata('small_exp4.txt');
small_data_5=importdata('small_exp5.txt');

large_data_1=importdata('large_exp1.txt');
large_data_2=importdata('large_exp2.txt');
large_data_3=importdata('large_exp3.txt');

figure
subplot(2,1,1)
title('Small beaker: Mass');  xlabel('Time (s)');  ylabel('Mass (g)');     legend;
hold on
% plot(small_data_1.data(1:61,1),small_data_1.data(1:61,3),'.')
% plot(small_data_2.data(1:61,1),small_data_2.data(1:61,3),'.')
% plot(small_data_3.data(1:61,1),small_data_3.data(1:61,3),'.')
% plot(small_data_4.data(1:61,1),small_data_4.data(1:61,3),'.')
plot(small_data_5.data(1:61,1),small_data_5.data(1:61,3),'.')
plot(t,m_water,'-')
legend('EX 5','Mod')
hold off

subplot(2,1,2)
title('Small beaker: Temperature (°C)');  xlabel('Time (s)');  ylabel('Temperature (°C)');     legend;
hold on
% plot(small_data_1.data(1:61,1),small_data_1.data(1:61,2),'.')
% plot(small_data_2.data(1:61,1),small_data_2.data(1:61,2),'.')
% plot(small_data_3.data(1:61,1),small_data_3.data(1:61,2),'.')
% plot(small_data_4.data(1:61,1),small_data_4.data(1:61,2),'.')
plot(small_data_5.data(1:61,1),small_data_5.data(1:61,2),'.')
plot(t,T_water,'-')
legend('EX 5','Mod')
hold off

% figure
% subplot(2,1,1)
% title('Large beaker: Mass');  xlabel('Time (s)');  ylabel('Mass (g)');     legend;
% hold on
% plot(large_data_1.data(1:61,1),large_data_1.data(1:61,3),'-o')
% plot(large_data_2.data(1:61,1),large_data_2.data(1:61,3),'-o')
% plot(large_data_3.data(1:61,1),large_data_3.data(1:61,3),'-o')
% plot(t,m_water,'-x')
% hold off
% 
% subplot(2,1,2)
% title('Large beaker: Temperature (°C)');  xlabel('Time (s)');  ylabel('Temperature (°C)');     legend;
% hold on
% plot(large_data_1.data(1:61,1),large_data_1.data(1:61,2),'-o')
% plot(large_data_2.data(1:61,1),large_data_2.data(1:61,2),'-o')
% plot(large_data_3.data(1:61,1),large_data_3.data(1:61,2),'-o')
% plot(t,T_water,'-x')
% hold off


exp=large_data_1.data(:,2);
model=T_water;

F_exp = griddedInterpolant(exp);
F_model = griddedInterpolant(model);

l_exp=length(large_data_1.data(:,2));
l_model=length(T_water);

out = [F_exp(linspace(1,numel(exp),l_model)'), F_model(linspace(1,numel(model),l_model)')];

diff=out(:,1)-out(:,2);
l=length(diff);
figure
hold on
plot(t,diff,'-o')
hold off
%% Renaming experimental data
T = small_data_1.data(:,2); % Temperature
m = small_data_1.data(:,3); % Dependent variable, m, mass
t = small_data_1.data(:,1); % Time

%% CALCULATE REGRESSION PARAMETERS
% Step 1: Normalize variables****************** 
% This step minimizes numerical problems when fitting polynomials. 
% In general, it is a good idea to normalize variables, especially in the 
% case of nonlinear regression, it helps convergence.

% Tn = (T-mean(T))./std(T);  % Normalized experimental temperature (-)
% mn = (m-mean(m))./std(m);  % Normalized experimental mass (-)

% Tm = (T-mean(T))./std(T);  % Normalized model temperature (-)
% mm = (m-mean(m))./std(m);  % Normalized model mass (-)

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
% SSr = sum((T-T).^2);           %Is this right?
MSr = SSr/(p-1);
Fobs = MSr/MSE;
alpha = 0.05;
Ftab = finv(1-alpha,p-1,n-p);




%% To do
%Check through all equations
% Better beta function?
%Check if Ra in T_surface_solve is less than 2*10^7
% Check different ODE?
%beta=1/T? For everything??
%R^2-value... the book was wrong, use the 1- version (wikipedia?). Use
%other indicators if value over 1 or under 0

% T_model=x1.*x.^2,x2.*x,x;
% tbl = table(T_model,out(:,1),'VariableNames',{'x','y'});
% tbl = table(exp);%out(:,1),out(:,2),'VariableNames',{'EXP Temp','MODEL Temp'});

% fitResults1 = polyfit(t,T_water,2)
% x1=fitResults1(1);
% x2=fitResults1(2);
% x3=fitResults1(3);


% mdl = fitlm(tbl,'y ~ x')
% mdl = fitlm(tbl,FORMULA)
% TBL = anova(mdl)
