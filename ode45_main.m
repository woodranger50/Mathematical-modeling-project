clear all, clear, clc;


t_span = [0 1800];
% m_water = 0.2;       %[kg]   Initial mass
T_water = 80+273.15;        %[K]   Initial temperature
% y0=[T_water,m_water];

% [time,Temperature,mass] = ode45(@(T_water,m_water) Good_version(T_water,m_water),t_span,T_water,m_water);
[t,dTdt]=ode45(@(t,dTdt)Good_version(T_water),t_span,T_water);

dTdt
% steps=length(T_water);
% t=linspace(t_span(1),t_span(2),steps);
T=dTdt-273.15;
figure
hold on
title('Temperature')
plot(t,T,'-o')
hold off

