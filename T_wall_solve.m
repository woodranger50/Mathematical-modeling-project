function F = T_wall_solve(T_wall)
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

%Additional data
T_surr=20+273.15;           %K ambient temperature

%% DECIDING VARIABLES
m_water=0.2;                             %[kg] mass of water in beaker
% Temperature=Temperature;          %[K] water temperature

d_inner=beaker_small(4);            %[m] diameter water
d_outer=beaker_small(3);            %[m]
%%
V_water=m_water/rho_water(T_water);         %[m^3] volume water 
h_water=V_water/(pi*(((d_inner/2)^2)));     %[m] height beaker
L_vert=h_water;                             %[m]
delta=d_outer-d_inner;                      %[m]

val=2.035*10^-9;    %taken at 298 K. Needs a function.)   

%Problem children
C=val;  %C is a value derived from appendix 1 equal to g*beta*(rho^2)/(my^2) with unit 1/K*m^3

%% fsolve Temp


% syms(sym('T_wall', [1 2]))


%Convection from water to inner wall 
Gr_L_water_wall=C*(L_vert^3)*(T_water-T_wall(1));   %where b is fluid coefficient of thermal expansion, g is grav. const., L is the significant length, dT is temp. diff., mu is fluid viscosity
Ra_L_water_wall=Gr_L_water_wall*pr_water(T_water);
Nu_L_water_wall=((0.825+0.387*(Ra_L_water_wall^(1/6)))/((1+(0.492/pr_water(T_water))^(9/16))^(8/27)))^2;
h_water_wall=Nu_L_water_wall*k_water(T_water)/L_vert;            %heat transfer coefficient for water-wall_inner

%Convection from outer wall to air
Gr_L_wall_air=C*(L_vert^3)*(T_wall(2)-T_surr);   %where b is fluid coefficient of thermal expansion, g is grav. const., L is the significant length, dT is temp. diff., mu is fluid viscosity
Ra_L_wall_air=Gr_L_wall_air*pr_water(T_water);
Nu_L_wall_air=((0.825+0.387*(Ra_L_wall_air^(1/6)))/((1+(0.492/pr_water(T_water))^(9/16))^(8/27)))^2;
h_wall_air=Nu_L_wall_air*k_water(T_water)/L_vert;            %heat transfer coefficient for wall_outer-air


T_h=T_water;
T_c=T_surr;
% T_wall(1)=T_wall_inner, T_wall(2)=T_wall_outer

F(1)=h_water_wall*(T_h-T_wall(1))-k_glass/delta*(T_wall(1)-T_wall(2));
F(2)=h_wall_air*(T_wall(2)-T_c)-k_glass/delta*(T_wall(1)-T_wall(2));




end