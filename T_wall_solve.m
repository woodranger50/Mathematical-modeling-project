function F = T_wall_solve(T_wall,T_water,m_water)
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
g=9.81;                     %m/s^2
%Stefan-Boltzmann Law
sigma=5.676*10^-8;          %W/m^2*K^4

%% DECIDING VARIABLES
d_inner=beaker_small(4);            %[m] diameter water
d_outer=beaker_small(3);            %[m]
%%
V_water=m_water/rho_water(T_water);         %[m^3] volume water 
h_water=V_water/(pi*(((d_inner/2)^2)));     %[m] height beaker
L_vert=h_water;                             %[m]
delta=d_outer-d_inner;                      %[m]
A_cyl_mantle=h_water*(pi*d_outer);          %[m^2]


%% fsolve of inside/outside wall temperatures
%Convection from water to inner wall 
beta_water_wall=(1-rho_water(T_water)/rho_water(T_wall(1)))/(T_water-T_wall(1));           %Maybe find another way to calculate beta
C_water_wall=(beta_water_wall*g*(rho_water(T_water)^2))/(my_water(T_water)^2);                      %[1/K*m^3]

Gr_L_water_wall=C_water_wall*(L_vert^3)*(T_water-T_wall(1));   %where b is fluid coefficient of thermal expansion, g is grav. const., L is the significant length, dT is temp. diff., mu is fluid viscosity
Ra_L_water_wall=Gr_L_water_wall*pr_water(T_water);
Nu_L_water_wall=((0.825+0.387*(Ra_L_water_wall^(1/6)))/((1+(0.492/pr_water(T_water))^(9/16))^(8/27)))^2;
h_water_wall=Nu_L_water_wall*k_water(T_water)/L_vert;            %heat transfer coefficient for water-wall_inner

%Convection from outer wall to air
beta_wall_air=(1-rho_water(T_water)/rho_water(T_wall(2)))/(T_water-T_wall(2));           %Maybe find another way to calculate beta
C_wall_air=(beta_wall_air*g*(rho_air(T_surr))^2)/(my_air(T_surr)^2);                      %[1/K*m^3]

Gr_L_wall_air=C_wall_air*(L_vert^3)*(T_wall(2)-T_surr);   %where b is fluid coefficient of thermal expansion, g is grav. const., L is the significant length, dT is temp. diff., mu is fluid viscosity
Ra_L_wall_air=Gr_L_wall_air*pr_water(T_water);
Nu_L_wall_air=((0.825+0.387*(Ra_L_wall_air^(1/6)))/((1+(0.492/pr_water(T_water))^(9/16))^(8/27)))^2;
h_wall_air=Nu_L_wall_air*k_water(T_water)/L_vert;            %heat transfer coefficient for wall_outer-air


T_h=T_water;    %Defines T_h as the temperature of the water
T_c=T_surr;     %Defines T_c as the surrounding temperature

F(1)=h_water_wall*(T_h-T_wall(1))-k_glass/delta*(T_wall(1)-T_wall(2));  %Function 1 
F(2)=h_wall_air*(T_wall(2)-T_c)-k_glass/delta*(T_wall(1)-T_wall(2))-epsilon_glass(1)*sigma*A_cyl_mantle*(T_wall(2))^4;


end