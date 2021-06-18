function [G,m_flow] = T_surface_solve(T_surface,T_water,m_water)
%% Calculating initial values

%Data for glass beakers are in ml, mm, mm, mm 
%and are size, height, outer diameter, and inner diameter in order
beaker_large=[800,135,98,93]*10^-3;   
beaker_small=[250,95,70,67]*10^-3;

%Heat conductivity GLASS
k_glass=0.9;                                %W/m*K
% k_polystyrene=0.03;
% k_glass=k_polystyrene;
%Emissivity GLASS
epsilon_glass=[0.94, 0.95];                 %[]
%Ambient temperature
T_surr=20+273.15;                           %[K]
%Molar mass WATER
n_water=18.01528*10^-3;                     %[kg/mol]
%Gravitational constant
g=9.81;                                     %[m/s^2]
%Stefan-Boltzmann constant
sigma=5.676*10^-8;                          %[W/(m^2*K^4)]
%Atmospheric pressure
p_atm=101.325*10^3;                         %[Pa] 
%Ideal gas constant
R=8.31446261815324;                         %[m^3*Pa/(K*mol)]

%% Geometric calculations

%Inner diameter beaker
d_inner=beaker_small(4);                    %[m]
%Outer diameter beaker
d_outer=beaker_small(3);                    %[m]
%Inner area surface of water (top of beaker)
A_inner=pi*(d_inner/2)^2;                   %[m^2]
%Volume of water
V_water=m_water/rho_water(T_water);         %[m^3]
%Height of water in beaker
h_water=V_water/A_inner;                    %[m]
%Characteristic length (for vertical wall)
L_vert=h_water;                             %[m]
%Beaker wall thickness
delta=d_outer-d_inner;                      %[m]
%Characteristic length (for horizontal plate)
L_surface=A_inner/(pi*d_inner);             %[m]

%% Convection

%The following calculations are for the convection from the bulk water to the water surface

%Fluid thermal coefficient of expansion
beta_water_surface=1/T_surface;             %[1/K]

%Possible alternatice equation for beta for vertical wall
%beta_water_surface=(1-rho_water(T_water)/rho_water(T_surface))/(T_water-T_surface);    %[1/K]

%Constant for Grashof
C_water_surface=beta_water_surface*g*(rho_water(T_surface)^2)/(my_water(T_surface)^2);  %[1/(m^3*K)]                      %[1/K*m^3]

%Grashof
Gr_L_water_surface=C_water_surface*(L_surface^3)*(T_water-T_surface);   %[]

%Rayleigh
Ra_L_water_surface=Gr_L_water_surface*pr_water(T_surface);              %[]

%Different equation Nussel for different values of Rayleigh
if Ra_L_water_surface<2*10^7
    Nu_L_water_surface=0.54*Ra_L_water_surface^(1/4);                   %[]
else
    Nu_L_water_surface=0.14*Ra_L_water_surface^(1/3);                   %[]
end

%Heat transfer coefficient
h_water_surface=Nu_L_water_surface*k_water(T_surface)/L_surface;        %[W/(m^2*K)]
%Heat transfer
dQdt_water_surface=h_water_surface*A_inner*(T_water-T_surface);         %[J]

%The following calculations are for the convection from the water surface
%to the bulk air

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
dQdt_surface_air=h_surface_air*A_inner*(T_surface-T_surr);                  %[J] heat transfer from water surface to air

%% Radiation
%Radiation from water surface to air
Radiation=0.95*sigma*A_inner*((T_surface)^4-(T_surr^4));  

%% Evaporation
%Evaporation from water surface to air
p_A=p_water(T_surface);                            %Partial pressure of water in bulk air
p_Ai=p_water(T_surr);                        %Partial pressure at film

D_AB=2.634/p_atm;                               %[m^2/s] Table J1

Sc=my_water(T_surr)/(rho_water(T_surr)*D_AB);
Pr=pr_water(T_surr);
% Pr_air=my_air(T_surface)*cp_air(T_surface)/k_air(T_surface);

%Chilton-Colburn analogy
k_c=h_water_surface*((Pr)^(2/3))/((Sc^(2/3)*rho_water(T_surr)*cp_water(T_surr))); %[m/s]
k_G=k_c/(R*T_surr);                          %[mol/s*m^2*Pa]

% Confirmation that Pr and Sc are calculated correctly
Le=Sc/Pr;
Le_conf=k_water(T_surr)/(rho_water(T_surr)*cp_water(T_surr)*D_AB);

N_A=k_G*(p_A-p_Ai);                     %[kg*mol/(m^2*s*Pa)]
m_flow=N_A*A_inner*n_water;                             %[kg/s] Mass flow
Evap_energy=m_flow*dHvap_water(T_surface);            %[J/s]

%% dQdt

G=(dQdt_water_surface)-dQdt_surface_air-Radiation-Evap_energy;                     %[J/s]

end