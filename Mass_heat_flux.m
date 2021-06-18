function [dt,Frac_evap,Frac_rad,Frac_conv]=Mass_heat_flux(t,ICs)
%% Calculating initial values
T_water=ICs(1);
m_water=ICs(2);
%Data for glass beakers are in ml, mm, mm, mm 
%and are size, height, outer diameter, and inner diameter in order
beaker_large=[800,135,98,93]*10^-3;    
beaker_small=[250,95,70,67]*10^-3;

%Data for glass material
rho_glass=2500;             %kg/m^3
cp_glass=840;               %J/kg*K
k_glass=0.9;                %W/m*K
% k_polystyrene=0.03;
% k_glass=k_polystyrene;
epsilon_glass=[0.94, 0.95]; %emissivity, maybe do average instead

%Additional data
RH=0.3;                     %relative humidity
T_surr=20+273.15;           %K ambient temperature
g=9.81;                     %m/s^2
%Stefan-Boltzmann Law
sigma=5.676*10^-8;          %W/m^2*K^4

pa=101.3*10^3;              %[Pa]

%% DECIDING VARIABLES
d_inner=beaker_small(4);            %[m] diameter water
d_outer=beaker_small(3);            %[m]

V_water=m_water/rho_water(T_water);         %[m^3] volume water 
h_water=V_water/(pi*(((d_inner/2)^2)));     %[m] height beaker
A_inner=pi*(d_inner/2)^2;                   %[m^2]
A_outer=pi*(d_outer/2)^2;                   %[m^2]
A_cyl_mantle=h_water*(pi*d_outer);          %[m^2]
A_cyl=A_cyl_mantle+2*(A_outer);             %[m^2]
A_cyl_mantle_inner=h_water*(pi*d_inner);    %[m^2]
A_cyl_mantle_outer=h_water*(pi*d_outer);    %[m^2]
L_hori=A_inner/(pi*d_inner);                %[m]
L_vert=h_water;                             %[m]
L_surface=A_inner/(pi*d_inner);             %[m]
delta=d_outer-d_inner;                      %[m]
m_glass_mantle=A_cyl_mantle*delta*rho_glass;%[kg]  

%% fsolve to find the inside and outside wall temperatures

OPTIONS=optimoptions(@fsolve,'Display','off');

x0=[T_water-3, T_water-6];
F=@(T_wall)T_wall_solve(T_wall,T_water,m_water);
T_wall=fsolve(F,x0,OPTIONS);

x1=T_water-10;                                              %Initial guess of temperatures
G=@(T_surface)T_surface_solve(T_surface,T_water,m_water);   %Declaring function
T_surface=fsolve(G,x1,OPTIONS);                             %Fsolve for temperature profile

T_wall_inner=T_wall(1);
T_wall_outer=T_wall(2);

%Film temperatures
T_film_water_wall=(T_water+T_wall_inner)/2;
T_film_wall_surr=(T_wall_outer+T_surr)/2;

%% Heat flux

%Convection from water to inner wall 
beta_water_wall=1/T_water;%(1-rho_water(T_water)/rho_water(T_wall_inner))/(T_water-T_wall_inner);           %Maybe find another way to calculate beta
C_water_wall=(beta_water_wall*g*(rho_water(T_water)^2))/(my_water(T_water)^2);                      %[1/K*m^3]

Gr_L_water_wall=C_water_wall*(L_vert^3)*(T_water-T_wall(1));   %where b is fluid coefficient of thermal expansion, g is grav. const., L is the significant length, dT is temp. diff., mu is fluid viscosity
Ra_L_water_wall=Gr_L_water_wall*pr_water(T_water);
Nu_L_water_wall=((0.825+0.387*(Ra_L_water_wall^(1/6)))/((1+(0.492/pr_water(T_water))^(9/16))^(8/27)))^2;
h_water_wall=Nu_L_water_wall*k_water(T_water)/L_vert;            %heat transfer coefficient for water-wall_inner

%Convection from water to water surface
beta_water_surface=1/T_water;%(1-rho_water(T_water)/rho_water(T_surface))/(T_water-T_surface);           %Maybe find another way to calculate beta
C_water_surface=beta_water_surface*g*(rho_water(T_water)^2)/(my_water(T_water)^2);                      %[1/K*m^3]

Gr_L_water_surface=C_water_surface*(L_surface^3)*(T_water-T_surface);   %where b is fluid coefficient of thermal expansion, g is grav. const., L is the significant length, dT is temp. diff., mu is fluid viscosity
Ra_L_water_surface=Gr_L_water_surface*pr_water(T_water);
Nu_L_water_surface=0.54*Ra_L_water_surface^(1/4);
h_water_surface=Nu_L_water_surface*k_water(T_water)/L_surface;            %heat transfer coefficient for water-wall_inner

%% Re-calculations for the physical comparrison 
%Convection from outer wall to air
beta_wall_air=1/T_film_wall_surr;%(1-rho_water(T_water)/rho_water(T_wall(2)))/(T_water-T_wall(2));           %Maybe find another way to calculate beta
C_wall_air=(beta_wall_air*g*(rho_air(T_film_wall_surr))^2)/(my_air(T_film_wall_surr)^2);                      %[1/K*m^3]

Gr_L_wall_air=C_wall_air*(L_vert^3)*(T_wall(2)-T_surr);   %where b is fluid coefficient of thermal expansion, g is grav. const., L is the significant length, dT is temp. diff., mu is fluid viscosity
Ra_L_wall_air=Gr_L_wall_air*pr_water(T_film_wall_surr);
Nu_L_wall_air=((0.825+0.387*(Ra_L_wall_air^(1/6)))/((1+(0.492/pr_water(T_film_wall_surr))^(9/16))^(8/27)))^2;
h_wall_air=Nu_L_wall_air*k_air(T_film_wall_surr)/L_vert;            %heat transfer coefficient for wall_outer-air

%Fluid thermal coefficient of expansion
beta_surface_air=1/T_surr;               %[1/K]

%Constant for Grashof
C_surface_air=beta_surface_air*g*(rho_air(T_surface)^2)/(my_air(T_surface)^2);                      %[1/K*m^3]

Gr_L_surface_air=C_surface_air*(L_surface^3)*(T_surface-T_surr);   %where b is fluid coefficient of thermal expansion, g is grav. const., L is the significant length, dT is temp. diff., mu is fluid viscosity
Ra_L_surface_air=Gr_L_surface_air*pr_water(T_surface);

%Different equation Nu for different Ra
if Ra_L_surface_air<2*10^7
    Nu_L_surface_air=0.54*Ra_L_surface_air^(1/4);
else
    Nu_L_surface_air=0.14*Ra_L_surface_air^(1/3);
end

h_surface_air=Nu_L_surface_air*k_air(T_surface)/L_surface;            %heat transfer coefficient for wall_outer-air

%%
dQdt_water_wall=h_water_wall*A_cyl_mantle*(T_water-T_wall_inner);    %[J/s] Heat transfer from water to inner wall
dQdt_water_surface=h_water_surface*A_inner*(T_water-T_surface);      %[J/s] Heat transfer from water to water surface

dQdt=dQdt_water_surface+dQdt_water_wall;        %[J/s]

%Differential equation T(t)
[G,m_flow]=T_surface_solve(T_surface,T_water,m_water);
dt= zeros(2,1);
dt(1)=-dQdt/(cp_water(T_water)*m_water);                       %[K/s]
dt(2)=-m_flow;                                                 %[kg/s]

%% Physical analysis
%Radiation
Radiation_surface=0.95*sigma*A_inner*((T_surface^4)-(T_surr^4));
Radiation_walls=A_cyl_mantle_outer*epsilon_glass(1)*sigma*((T_wall(2)^4)-(T_surr^4));

%Convection
Convection_surface=h_surface_air*A_inner*(T_surface-T_surr);
Convection_walls=h_wall_air*A_cyl_mantle_outer*(T_wall_outer-T_surr);

%Evaporation
Evaporation=m_flow*dHvap_water(T_surface);

% Fractions
Radiation=Radiation_surface+Radiation_walls;
Convection=Convection_surface+Convection_walls;
Sum=Radiation+Convection+Evaporation;

Frac_rad=Radiation/Sum;
Frac_conv=Convection/Sum;
Frac_evap=Evaporation/Sum;

end
