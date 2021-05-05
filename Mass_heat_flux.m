function dTdt=Mass_heat_flux(t,T_water)
%% Calculating initial values
%Calc radiation wall-air
%Calc radiation water-air
%Fix convection for water surface

%Data for glass beakers are in ml, mm, mm, mm 
%and are size, height, outer diameter, and inner diameter in order
beaker_large=[800,135,98,93]*10^-3;   
beaker_small=[250,95,70,67]*10^-3;

%Data for glass material
rho_glass=2500;             %kg/m^3
cp_glass=840;               %J/kg*K
k_glass=0.9;                %W/m*K
epsilon_glass=[0.94, 0.95]; %emissivity, maybe do average instead

%Additional data
RH=0.3;                     %relative humidity
T_surr=20+273.15;           %K ambient temperature
g=9.81;                     %m/s^2
%Stefan-Boltzmann Law
sigma=5.676*10^-8;          %W/m^2*K^4

pa=101.3*10^3;              %[Pa]

%% DECIDING VARIABLES
m_water=0.2;                             %[kg] mass of water in beaker

d_inner=beaker_small(4);            %[m] diameter water
d_outer=beaker_small(3);            %[m]
%%
V_water=m_water/rho_water(T_water);         %[m^3] volume water 
h_water=V_water/(pi*(((d_inner/2)^2)));     %[m] height beaker
A_inner=pi*(d_inner/2)^2;                   %[m^2]
A_outer=pi*(d_outer/2)^2;                   %[m^2]
A_cyl_mantle=h_water*(pi*d_outer);          %[m^2]
A_cyl=A_cyl_mantle+2*(A_outer);             %[m^2]
L_hori=A_inner/(pi*d_inner);                %[m]
L_vert=h_water;                             %[m]
L_surface=A_inner/(pi*d_inner);             %[m]
delta=d_outer-d_inner;                      %[m]
m_glass_mantle=A_cyl_mantle*delta*rho_glass;%[kg]    

%fsolve to find the inside and outside wall temperatures
F=@T_wall_solve;
x0=[330 320];
OPTIONS=optimoptions(@fsolve,'Display','off');
[T_wall,y]=fsolve(F,x0,OPTIONS);
T_wall_inner=T_wall(1);
T_wall_outer=T_wall(2);

G=@T_surface_solve;                                     %Declaring function
x1=320;                                                %Initial guess of temperatures
[T_surface,y]=fsolve(G,x1,OPTIONS);                     %Fsolve for temperature profile

%% Heat flux

%Convection from water to inner wall 
beta_water_wall=(1-rho_water(T_water)/rho_water(T_wall_inner))/(T_water-T_wall_inner);           %Maybe find another way to calculate beta
C_water_wall=(beta_water_wall*g*(rho_water(T_water)^2))/(my_water(T_water)^2);                      %[1/K*m^3]

Gr_L_water_wall=C_water_wall*(L_vert^3)*(T_water-T_wall(1));   %where b is fluid coefficient of thermal expansion, g is grav. const., L is the significant length, dT is temp. diff., mu is fluid viscosity
Ra_L_water_wall=Gr_L_water_wall*pr_water(T_water);
Nu_L_water_wall=((0.825+0.387*(Ra_L_water_wall^(1/6)))/((1+(0.492/pr_water(T_water))^(9/16))^(8/27)))^2;
h_water_wall=Nu_L_water_wall*k_water(T_water)/L_vert;            %heat transfer coefficient for water-wall_inner

%Horizontal plates (with sides?)
beta_water_surface=(1-rho_water(T_water)/rho_water(T_surface))/(T_water-T_surface);           %Maybe find another way to calculate beta
C_water_surface=beta_water_surface*g*(rho_water(T_water)^2)/(my_water(T_water)^2);                      %[1/K*m^3]

Gr_L_water_surface=C_water_surface*(L_surface^3)*(T_water-T_surface);   %where b is fluid coefficient of thermal expansion, g is grav. const., L is the significant length, dT is temp. diff., mu is fluid viscosity
Ra_L_water_surface=Gr_L_water_surface*pr_water(T_water);
Nu_L_water_surface=((0.825+0.387*(Ra_L_water_surface^(1/6)))/((1+(0.492/pr_water(T_water))^(9/16))^(8/27)))^2;
h_water_surface=Nu_L_water_surface*k_water(T_water)/L_surface;            %heat transfer coefficient for water-wall_inner


%%

dQdt_water_wall=h_water_wall*A_cyl*(T_water-T_wall_inner);                  %[J]
dQdt_water_surface=h_water_surface*A_inner*(T_water-T_surface);             %[J] heat transfer from water to water surface

dQdt=dQdt_water_surface+dQdt_water_wall;

%Differential equation T(t)
dTdt=-dQdt/(cp_water(T_water)*m_water);                       %[K/s]
end
