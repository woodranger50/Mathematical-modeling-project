clear all, clear, clc;

t_span = [0 1800];
% m_water = 0.2;              %[kg]   Initial mass
T_water = 80+273.15;        %[K]   Initial temperature
% y0=[T_water,m_water];

% [time,Temperature,mass] = ode45(@(T_water,m_water) Good_version(T_water,m_water),t_span,T_water,m_water);
[t,dTdt]=ode45(@(t,dTdt)Mass_heat_flux(T_water),t_span,T_water);

% steps=length(T_water);
% t=linspace(t_span(1),t_span(2),steps);
T=dTdt-273.15;
figure
hold on
title('Temperature'); xlabel('Time (seconds)'); ylabel('Temperature (°C)')
plot(t,T,'-o')
hold off

%% To do
% Make T_water uppdate through time
% Better beta function?
% Make m_water change
% Why is temp -500? Debugg...
% Check different ODE?