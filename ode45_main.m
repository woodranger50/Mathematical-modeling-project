clear all, clear, clc;
format long
t_span = [0 1800];
m_water = 0.2;              %[kg]   Initial mass
T_water = 80+273.15;        %[K]   Initial temperature
y0=[T_water,m_water];


[t,dt]=ode45(@Mass_heat_flux,t_span,y0);

T=dt(:,1);
T_water=T-273.15;
m=dt(:,2);
m_water=m*10^3;

figure
hold on
subplot(1,2,1)
title('Temperature'); xlabel('Time (seconds)'); ylabel('Temperature [°C]')
plot(t,T_water,'-o')
subplot(1,2,2)
title('Mass'); xlabel('Time (seconds)'); ylabel('Mass [g]')
plot(t,m_water,'-o')
hold off

%% To do
%Check through all equations
% Better beta function?
%Check if Ra in T_surface_solve is less than 2*10^7
% Check different ODE?
%beta=1/T? For everything??