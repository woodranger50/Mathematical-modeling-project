function [G,m_flow] = T_surface_solve(T_surface,T_water,m_water)
%% Calculating initial values
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
n_water=18.01528*10^-3;     %[kg/mol] Molar mass
g=9.81;                     %m/s^2

%Stefan-Boltzmann Law
sigma=5.676*10^-8;          %W/m^2*K^4
p_atm=101.3*10^3;              %[Pa] 
R=8.31446261815324;         %[m^3*Pa/(K*mol)]
%% DECIDING VARIABLES
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

%% fsolve of surface temperature
%Convection from water to water surface 
beta_water_surface=1/T_water;%(1-rho_water(T_water)/rho_water(T_surface))/(T_water-T_surface);           %Maybe find another way to calculate beta
C_water_surface=beta_water_surface*g*(rho_water(T_water)^2)/(my_water(T_water)^2);                      %[1/K*m^3]

% Gr=(beta*g*rho^2*L^3*dT)/(my^2)


Gr_L_water_surface=C_water_surface*(L_surface^3)*(T_water-T_surface);   %where b is fluid coefficient of thermal expansion, g is grav. const., L is the significant length, dT is temp. diff., mu is fluid viscosity
Ra_L_water_surface=Gr_L_water_surface*pr_water(T_water);
Nu_L_water_surface=0.54*Ra_L_water_surface^(1/4);
h_water_surface=Nu_L_water_surface*k_water(T_water)/L_surface;            %heat transfer coefficient for water-wall_inner
dQdt_water_surface=h_water_surface*A_inner*(T_water-T_surface);             %[J] heat transfer from water to water surface

%Convection from water surface to air
beta_surface_air=1/T_surr;%(1-rho_air(T_surr)/rho_water(T_surface))/(T_surface-T_surr);           %Maybe find another way to calculate beta
C_surface_air=beta_surface_air*g*(rho_air(T_surr)^2)/(my_air(T_surr)^2);                      %[1/K*m^3]

Gr_L_surface_air=C_surface_air*(L_surface^3)*(T_surface-T_surr);   %where b is fluid coefficient of thermal expansion, g is grav. const., L is the significant length, dT is temp. diff., mu is fluid viscosity
Ra_L_surface_air=Gr_L_surface_air*pr_water(T_surface);
Nu_L_surface_air=0.54*Ra_L_surface_air^(1/4);
h_surface_air=Nu_L_surface_air*k_air(T_surr)/L_surface;            %heat transfer coefficient for wall_outer-air
dQdt_surface_air=h_surface_air*A_inner*(T_surface-T_surr);                  %[J] heat transfer from water surface to air

%Radiation from water surface to air
Radiation=epsilon_glass(1)*sigma*A_inner*(T_surface)^4;  

%Evaporation from water surface to air
p_A=p_water(T_water);                            %Partial pressure of water in bulk air
p_Ai=p_water(T_surface);                        %Partial pressure at film                       ITS GREAT maybe?

my_a=my_air(T_surr);                         %[Pa*s]    
rho_a=rho_air(T_surr);                       %[kg/m^3]
D_AB=2.634/p_atm;                               %[m^2/s] Table J1
v=0.3;                                          %[m/s] THIS IS A GUESS
Re_L=rho_a*v*L_surface/my_a;                    %Dim less :)
Sc=my_a/(rho_a*D_AB);                           %Dim less
k_c=0.664*(D_AB/L_surface)*Re_L^(1/2)*Sc^(1/3); %[m/s]

k_G=k_c/(R*T_surface);                          %[mol/s*m^2*Pa]

N_A=k_G*A_inner*(p_A-p_Ai);                     %[mol/s]
m_flow=N_A*n_water;                             %[kg/s] Mass flow
Evap_energy=m_flow*dHvap_water(T_surface);            %[J/s]


G=(dQdt_water_surface)-dQdt_surface_air-Radiation-Evap_energy;                     %[J/s^2]




end