clear all, clear, clc;
format long
t_span = [0 1800];
m_water = 0.1795;              %[kg]   Initial mass
T_water = 80+273.15;        %[K]   Initial temperature
y0=[T_water,m_water];


[t,dt]=ode45(@Mass_heat_flux,t_span,y0);

matrix=zeros(41,3);
for i=1:41
[~,Frac_evap,Frac_rad,Frac_conv]=Mass_heat_flux(0,[dt(i,1) dt(i,2)]);
matrix(i,1)=Frac_evap;
matrix(i,2)=Frac_rad;
matrix(i,3)=Frac_conv;
end

T_mod=dt(:,1);
T_water=T_mod-273.15;
m_mod=dt(:,2);
m_water=m_mod*10^3;

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
xlabel('Time (s)');  ylabel('Mass (g)');
hold on
% plot(small_data_1.data(1:61,1),small_data_1.data(1:61,3),'x')
% plot(small_data_2.data(1:61,1),small_data_2.data(1:61,3),'x')
% plot(small_data_3.data(1:61,1),small_data_3.data(1:61,3),'x')
% plot(small_data_4.data(1:61,1),small_data_4.data(1:61,3),'x')
plot(small_data_5.data(1:61,1),small_data_5.data(1:61,3),'x')
plot(t,m_water,'-')
legend('Experimental data','Model data')
hold off

subplot(2,1,2)
xlabel('Time (s)');  ylabel('Temperature (°C)');
hold on
% plot(small_data_1.data(1:61,1),small_data_1.data(1:61,2),'x')
% plot(small_data_2.data(1:61,1),small_data_2.data(1:61,2),'x')
% plot(small_data_3.data(1:61,1),small_data_3.data(1:61,2),'x')
% plot(small_data_4.data(1:61,1),small_data_4.data(1:61,2),'x')
plot(small_data_5.data(1:61,1),small_data_5.data(1:61,2),'x')
plot(t,T_water,'-')
legend('Experimental data','Model data')
hold off

% figure
% subplot(2,1,1)
% title('Large beaker: Mass');  xlabel('Time (s)');  ylabel('Mass (g)');     legend;
% hold on
% plot(large_data_1.data(1:61,1),large_data_1.data(1:61,3),'x')
% % plot(large_data_2.data(1:61,1),large_data_2.data(1:61,3),'x')
% % plot(large_data_3.data(1:61,1),large_data_3.data(1:61,3),'x')
% plot(t,m_water,'-')
% hold off
% 
% subplot(2,1,2)
% title('Large beaker: Temperature (°C)');  xlabel('Time (s)');  ylabel('Temperature (°C)');     legend;
% hold on
% plot(large_data_1.data(1:61,1),large_data_1.data(1:61,2),'x')
% % plot(large_data_2.data(1:61,1),large_data_2.data(1:61,2),'x')
% % plot(large_data_3.data(1:61,1),large_data_3.data(1:61,2),'x')
% plot(t,T_water,'-')
% hold off

%% Renaming experimental data
T = small_data_1.data(:,2); % Temperature
T_2 = small_data_2.data(:,2); % Temperature
T_3 = small_data_3.data(:,2); % Temperature
T_4 = small_data_4.data(:,2); % Temperature
T_5 = small_data_5.data(:,2); % Temperature

m = small_data_1.data(:,3); % Dependent variable, m, mass
m_2 = small_data_2.data(:,3); % Dependent variable, m, mass
m_3 = small_data_3.data(:,3); % Dependent variable, m, mass
m_4 = small_data_4.data(:,3); % Dependent variable, m, mass
m_5 = small_data_5.data(:,3); % Dependent variable, m, mass
time = small_data_1.data(:,1); % Time

%% Polyfit model temperature and mass
fitResults_T = polyfit(t,T_water,2);
fitResults_m = polyfit(t,m_water,2);

y_T=@(x)fitResults_T(1)*x.^2+fitResults_T(2)*x+fitResults_T(3);
y_m=@(x)fitResults_m(1)*x.^2+fitResults_m(2)*x+fitResults_m(3);

T_model=y_T(time);
m_model=y_m(time);

diff_T=(T-T_model)/std(T_model);
diff_m=(m-m_model)/std(m_model);



figure
subplot(2,1,1)
xlabel('Time (s)');   ylabel('Temperature (°C)');    ylim([-2.5 2.5]);
hold on
plot(time,1,'-')
plot(time,diff_T,'x')
plot([0,1800],[0,0],'k')
plot([0,1800],[1.96,1.96],'r')
plot([0,1800],[-1.96,-1.96],'r')
hold off
subplot(2,1,2)
xlabel('Time (s)');   ylabel('Mass (g)');   ylim([-2.5 2.5]);
hold on
plot(time,1,'-')
plot(time,diff_m,'x')
plot([0,1800],[0,0],'k')
plot([0,1800],[1.96,1.96],'r')
plot([0,1800],[-1.96,-1.96],'r')
hold off

%% CALCULATE REGRESSION PARAMETERS
% Step 1: Normalize variables****************** 
% This step minimizes numerical problems when fitting polynomials. 
% In general, it is a good idea to normalize variables, especially in the 
% case of nonlinear regression, it helps convergence.

% Tn_exp_1 = (T-mean(T))./std(T);  % Normalized experimental temperature (-)
% Tn_exp_2 = (T_2-mean(T_2))./std(T_2);  % Normalized experimental temperature (-)
% Tn_exp_3 = (T_3-mean(T_3))./std(T_3);  % Normalized experimental temperature (-)
% Tn_exp_4 = (T_4-mean(T_4))./std(T_4);  % Normalized experimental temperature (-)
% Tn_exp_5 = (T_5-mean(T_5))./std(T_5);  % Normalized experimental temperature (-)
% 
% mn_exp_1 = (m-mean(m))./std(m);  % Normalized experimental mass (-)
% mn_exp_2 = (m-mean(m_2))./std(m_2);  % Normalized experimental mass (-)
% mn_exp_3 = (m-mean(m_3))./std(m_3);  % Normalized experimental mass (-)
% mn_exp_4 = (m-mean(m_4))./std(m_4);  % Normalized experimental mass (-)
% mn_exp_5 = (m-mean(m_5))./std(m_5);  % Normalized experimental mass (-)
% 
% Tn_model = (T-mean(T))./std(T);  % Normalized model temperature (-)
% mn_model = (m-mean(m))./std(m);  % Normalized model mass (-)
% 
% T=Tn_exp_1;
% T_2=Tn_exp_2;
% T_3=Tn_exp_3;
% T_4=Tn_exp_4;
% T_5=Tn_exp_5;
% 
% m=mn_exp_1;
% m_2=mn_exp_2;
% m_3=mn_exp_3;
% m_4=mn_exp_4;
% m_5=mn_exp_5;

% Step 2: Calculate regression parameters********

p = 2;          % number of regression parameters
m_nr=5;         %nr of data ponts times number of similar experiments
n = 61*m_nr;          % number of data points
v=n-p;           % Degrees of freedom
v_reg=p-1; 
v_lof=m_nr-p;
v_pe=n-m_nr;
v_t=n-1;
% X = [ones(length(small_data_1.data(:,1)),1),small_data_1.data(:,1)];          % Define matrix X

% beta = inv(X'*X)*X'*T;       % Calculate regression parameters

%% DETERMINE IF AT LEAST ONE REGRESSION PARAMETER IS STATISTICALLY SIGNIFICANT
% Hint: See MATLAB program in assignment 27
alpha = 0.05;

%Tempterature
% T_bar = (T+T_2+T_3+T_4+T_5)./5;  % Average 
T_bar=(1/length(T_5))*(sum(T_5))
T_bar_1 = abs(mean(T));
%SS
SStot_T=sum((T_5-T_bar).^2)
SSres_T=sum((T_model-T_bar).^2)
SSreg_T = sum((T_model-T_bar_1).^2);
SSE_T = ((sum((T_model-T).^2))+(sum((T_model-T_2).^2))+(sum((T_model-T_3).^2))...
        +(sum((T_model-T_4).^2))+(sum((T_model-T_5).^2))); % See Appendix
SSPe_T=((sum((T-T_bar).^2))+(sum((T_2-T_bar).^2))+(sum((T_3-T_bar).^2))...
        +(sum((T_4-T_bar).^2))+(sum((T_5-T_bar).^2)));                             % Check later
SSlof_T=SSE_T-SSPe_T;
SSt_T=SSreg_T+SSE_T;

%MS
MSreg_T = SSreg_T/v_reg;
s_squared_T=SSE_T/v;
MSlof_T=SSlof_T/v_lof;
MSE_T = SSE_T/v;   % Mean square error
MSpe_T=SSPe_T/v_pe;

s_T=sqrt(s_squared_T);   %Standard deviation
Fobs_T = MSreg_T/s_squared_T;
Flof_T=MSlof_T/MSpe_T;
Ftab_T = finv(1-alpha,p-1,n-p);
R_sq_T=SStot_T/SSres_T

%Mass
m_bar = (m+m_2+m_3+m_4+m_5)./5;  % Average 
m_bar_1 = abs(mean(m));

%SS
SSreg_m = sum((m_model-m_bar_1).^2);
SSE_m = ((sum((m_model-m).^2))+(sum((m_model-m_2).^2))+(sum((m_model-m_3).^2))...
        +(sum((m_model-m_4).^2))+(sum((m_model-m_5).^2))); % See Appendix
SSPe_m=((sum((m-m_bar).^2))+(sum((m_2-m_bar).^2))+(sum((m_3-m_bar).^2))...
    +(sum((m_4-m_bar).^2))+(sum((m_5-m_bar).^2)));                             % Check later
SSlof_m=SSE_m-SSPe_m;
SSt_m=SSreg_m+SSE_m;

%MS
MSreg_m = SSreg_m/v_reg;
s_squared_m=SSE_m/v;
MSlof_m=SSlof_m/v_lof;
MSE_m = SSE_m/v;   % Mean square error
MSpe_m=SSPe_m/v_pe;

s_m=sqrt(s_squared_m);   %Standard deviation
Fobs_m = MSreg_m/s_squared_m;
Flof_m=MSlof_m/MSpe_m;
Ftab_m = finv(1-alpha,p-1,n-p);
R_sq_m=SSreg_m/(SSreg_m+SSE_m);

format short

Source = {'Regression'; 'Residual'};
SS = {SSreg_T; SSE_T};
df = {v_reg; v};
MS = {MSreg_T; s_squared_T};
Fobs = {Fobs_T; []};
Table_T = table(Source, SS, df, MS, Fobs);
% disp(Table_T)

Source = {'Regression'; 'Residual'};
SS = {SSreg_m; SSE_m};
df = {v_reg; v};
MS = {MSreg_m; s_squared_m};
Fobs = {Fobs_m; []};
Table_m = table(Source, SS, df, MS, Fobs);
% disp(Table_m)

%% Physical analysis
matrix_sum=matrix(:,1)+matrix(:,2)+matrix(:,3);

figure
hold on
xlabel('Time (s)');  ylabel('Fraction'); ylim([0 1.1]);
plot(t,matrix(:,1),'-')
plot(t,matrix(:,2),'-')
plot(t,matrix(:,3),'-')
plot(t,matrix_sum,'k')
legend('Evaporation','Radiation','Convection','Sum of Fractions')
hold off

%% To do
%Check through all equations
% Better beta function?
%Check if Ra in T_surface_solve is less than 2*10^7
% Check different ODE?
%beta=1/T? For everything??
%R^2-value... the book was wrong, use the 1- version (wikipedia?). Use
%other indicators if value over 1 or under 0

%Why? And also Fractions? Radiation going up?

%%
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

% exp=large_data_1.data(:,2);
% model=T_water;
% 
% F_exp = griddedInterpolant(exp);
% F_model = griddedInterpolant(model);
% 
% l_exp=length(large_data_1.data(:,2));
% l_model=length(T_water);
% 
% out = [F_exp(linspace(1,numel(exp),l_model)'), F_model(linspace(1,numel(model),l_model)')];
% 
% diff=out(:,1)-out(:,2);
% l=length(diff);
% figure
% hold on
% plot(t,diff,'-o')
% hold off
