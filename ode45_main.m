clear all, clear, clc;
format long
t_span = 0:45:1800;
m_water = 0.57538;              %[kg]   Initial mass
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

% figure
% subplot(2,1,1)
% title('Small beaker: Mass');  xlabel('Time (s)');  ylabel('Mass (g)');     legend;
% hold on
% plot(small_data_1.data(1:61,1),small_data_1.data(1:61,3))
% plot(small_data_2.data(1:61,1),small_data_2.data(1:61,3))
% plot(small_data_3.data(1:61,1),small_data_3.data(1:61,3))
% plot(small_data_4.data(1:61,1),small_data_4.data(1:61,3))
% plot(small_data_5.data(1:61,1),small_data_5.data(1:61,3))
% plot(t,m_water,'-o')
% legend('EX 1','EX 2','EX 3','EX 4','EX 5','Mod')
% hold off
% 
% subplot(2,1,2)
% title('Small beaker: Temperature (°C)');  xlabel('Time (s)');  ylabel('Temperature (°C)');     legend;
% hold on
% plot(small_data_1.data(1:61,1),small_data_1.data(1:61,2))
% plot(small_data_2.data(1:61,1),small_data_2.data(1:61,2))
% plot(small_data_3.data(1:61,1),small_data_3.data(1:61,2))
% plot(small_data_4.data(1:61,1),small_data_4.data(1:61,2))
% plot(small_data_5.data(1:61,1),small_data_5.data(1:61,2))
% plot(t,T_water,'-o')
% legend('EX 1','EX 2','EX 3','EX 4','EX 5','Mod')
% hold off

figure
subplot(2,1,1)
title('Large beaker: Mass');  xlabel('Time (s)');  ylabel('Mass (g)');     legend;
hold on
plot(large_data_1.data(1:61,1),large_data_1.data(1:61,3),'-o')
plot(large_data_2.data(1:61,1),large_data_2.data(1:61,3),'-o')
plot(large_data_3.data(1:61,1),large_data_3.data(1:61,3),'-o')
plot(t,m_water,'-x')
hold off

subplot(2,1,2)
title('Large beaker: Temperature (°C)');  xlabel('Time (s)');  ylabel('Temperature (°C)');     legend;
hold on
plot(large_data_1.data(1:61,1),large_data_1.data(1:61,2),'-o')
plot(large_data_2.data(1:61,1),large_data_2.data(1:61,2),'-o')
plot(large_data_3.data(1:61,1),large_data_3.data(1:61,2),'-o')
plot(t,T_water,'-x')
hold off

fitResults1 = polyfit(t,T_water,2)
x1=fitResults1(1);
x2=fitResults1(2);
x3=fitResults1(3);
% FORMULA='quadratic';
x=T_water;
exp=large_data_1.data(:,2);
model=T_water;
F_exp = griddedInterpolant(exp);
F_model = griddedInterpolant(model);
l_exp=length(large_data_1.data(:,2));
l_model=length(T_water);
out = [F_exp(linspace(1,numel(exp),l_model)'), F_model(linspace(1,numel(model),l_model)')];
T_model=x1.*x.^2,x2.*x,x;
tbl = table(T_model,out(:,1),'VariableNames',{'x','y'});
% tbl = table(exp);%out(:,1),out(:,2),'VariableNames',{'EXP Temp','MODEL Temp'});

mdl = fitlm(tbl,'y ~ x')
% mdl = fitlm(tbl,FORMULA)
% TBL = anova(mdl)
%% To do
%Check through all equations
% Better beta function?
%Check if Ra in T_surface_solve is less than 2*10^7
% Check different ODE?
%beta=1/T? For everything??
%R^2-value... the book was wrong, use the 1- version (wikipedia?). Use
%other indicators if value over 1 or under 0