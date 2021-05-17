function F = T_wall_solve(T_wall,T_water,m_water)
%% Calculating initial values

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
d_inner=beaker_large(4);            %[m] diameter water
d_outer=beaker_large(3);            %[m]

V_water=m_water/rho_water(T_water);         %[m^3] volume water 
h_water=V_water/(pi*(((d_inner/2)^2)));     %[m] height beaker
L_vert=h_water;                             %[m]
delta=d_outer-d_inner;                      %[m]
A_cyl_mantle_inner=h_water*(pi*d_inner);    %[m^2]
A_cyl_mantle_outer=h_water*(pi*d_outer);    %[m^2] 

%Film temperatures
T_film_water_wall=(T_water+T_wall(1))/2;
T_film_wall_surr=(T_wall(2)+T_surr)/2;

%% fsolve of inside/outside wall temperatures
%Convection from water to inner wall 
beta_water_wall=1/T_film_water_wall;%(1-rho_water(T_water)/rho_water(T_wall(1)))/(T_water-T_wall(1));           %Maybe find another way to calculate beta
C_water_wall=(beta_water_wall*g*(rho_water(T_film_water_wall)^2))/(my_water(T_film_water_wall)^2);                      %[1/K*m^3]

Gr_L_water_wall=C_water_wall*(L_vert^3)*(T_water-T_wall(1));   %where b is fluid coefficient of thermal expansion, g is grav. const., L is the significant length, dT is temp. diff., mu is fluid viscosity
Ra_L_water_wall=Gr_L_water_wall*pr_water(T_film_water_wall);
Nu_L_water_wall=((0.825+0.387*(Ra_L_water_wall^(1/6)))/((1+(0.492/pr_water(T_film_water_wall))^(9/16))^(8/27)))^2;
h_water_wall=Nu_L_water_wall*k_water(T_film_water_wall)/L_vert;            %[W/(m^2*K)]heat transfer coefficient for water-wall_inner

%Convection from outer wall to air
beta_wall_air=1/T_film_wall_surr;%(1-rho_water(T_water)/rho_water(T_wall(2)))/(T_water-T_wall(2));           %Maybe find another way to calculate beta
C_wall_air=(beta_wall_air*g*(rho_air(T_film_wall_surr))^2)/(my_air(T_film_wall_surr)^2);                      %[1/K*m^3]

Gr_L_wall_air=C_wall_air*(L_vert^3)*(T_wall(2)-T_surr);   %where b is fluid coefficient of thermal expansion, g is grav. const., L is the significant length, dT is temp. diff., mu is fluid viscosity
Ra_L_wall_air=Gr_L_wall_air*pr_water(T_film_wall_surr);
Nu_L_wall_air=((0.825+0.387*(Ra_L_wall_air^(1/6)))/((1+(0.492/pr_water(T_film_wall_surr))^(9/16))^(8/27)))^2;
h_wall_air=Nu_L_wall_air*k_air(T_film_wall_surr)/L_vert;            %heat transfer coefficient for wall_outer-air

%% dQdt

dQdt_water_inner_wall=h_water_wall*A_cyl_mantle_inner*(T_water-T_wall(1));
dQdt_outer_wall_air=h_wall_air*A_cyl_mantle_outer*(T_wall(2)-T_surr);
dQdt_radiation=A_cyl_mantle_outer*epsilon_glass(1)*sigma*(T_wall(2))^4;
dQdt_conduction=k_glass*A_cyl_mantle_inner/delta*(T_wall(1)-T_wall(2));


F(1)=dQdt_water_inner_wall-dQdt_conduction; 
F(2)=dQdt_conduction-dQdt_outer_wall_air-dQdt_radiation;


end