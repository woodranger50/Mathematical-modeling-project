function G = T_surface_solve(T_surface)
%% Calculating initial values
T_water=273.15+80;
%Calc radiation wall-air
%Calc radiation water-air
%Fix convection for water surface

%Data for glass beakers are in ml, mm, mm, mm 
%and are size, height, outer diameter, and inner diameter in order
beaker_large=[800,135,98,93]*10^-3;   
beaker_small=[250,95,70,67]*10^-3;

%Data for glass material
k_glass=0.9;                %W/m*K
epsilon_glass=[0.94, 0.95]; %emissivity, maybe do average instead

%Additional data
T_surr=20+273.15;           %K ambient temperature
%Stefan-Boltzmann Law
sigma=5.676*10^-8;          %W/m^2*K^4

%% DECIDING VARIABLES
m_water=0.2;                             %[kg] mass of water in beaker

d_inner=beaker_small(4);            %[m] diameter water
d_outer=beaker_small(3);            %[m]
%%
V_water=m_water/rho_water(T_water);         %[m^3] volume water 
h_water=V_water/(pi*(((d_inner/2)^2)));     %[m] height beaker
A_inner=pi*(d_inner/2)^2;                   %[m^2]
L_vert=h_water;                             %[m]
delta=d_outer-d_inner;                      %[m]
A_cyl_mantle=h_water*(pi*d_outer);          %[m^2]
L_surface=A_inner/(pi*d_inner);             %[m]

val=2.035*10^-9;    %taken at 298 K. Needs a function.)   

%Problem children
C=val;  %C is a value derived from appendix 1 equal to g*beta*(rho^2)/(my^2) with unit 1/K*m^3

%% fsolve of surface temperature
%Convection from water to water surface 
Gr_L_water_surface=C*(L_surface^3)*(T_water-T_surface);   %where b is fluid coefficient of thermal expansion, g is grav. const., L is the significant length, dT is temp. diff., mu is fluid viscosity
Ra_L_water_surface=Gr_L_water_surface*pr_water(T_water);
Nu_L_water_surface=((0.825+0.387*(Ra_L_water_surface^(1/6)))/((1+(0.492/pr_water(T_water))^(9/16))^(8/27)))^2;
h_water_surface=Nu_L_water_surface*k_water(T_water)/L_surface;            %heat transfer coefficient for water-wall_inner
dQdt_water_surface=h_water_surface*A_inner*(T_water-T_surface);             %[J] heat transfer from water to water surface

%Convection from water surface to air
Gr_L_surface_air=C*(L_surface^3)*(T_surface-T_surr);   %where b is fluid coefficient of thermal expansion, g is grav. const., L is the significant length, dT is temp. diff., mu is fluid viscosity
Ra_L_surface_air=Gr_L_surface_air*pr_water(T_water);
Nu_L_surface_air=((0.825+0.387*(Ra_L_surface_air^(1/6)))/((1+(0.492/pr_water(T_water))^(9/16))^(8/27)))^2;
h_surface_air=Nu_L_surface_air*k_water(T_water)/L_surface;            %heat transfer coefficient for wall_outer-air
dQdt_surface_air=h_surface_air*A_inner*(T_surface-T_surr);                  %[J] heat transfer from water surface to air

%Radiation from water surface to air
Radiation=epsilon_glass(1)*sigma*A_inner*(T_surface)^4;  

%Evaporation from water surface to air
%Coming soon

G=(dQdt_water_surface)-(dQdt_surface_air)-(Radiation);                     %[J/s^2]




end