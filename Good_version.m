function T_water=Good_version(T_water)
%% Calculating initial values
%Data for glass beakers are in ml, mm, mm, mm 
%and are size, height, outer diameter, and inner diameter in order
beaker_large=[800,135,98,93*10^-3];   
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
% Temperature=Temperature;          %[K] water temperature

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

delta=d_inner-d_outer;                      %[m]

T_film=(T_water+T_surr)/2;  %[K]

val=2.035*10^-9;    %taken at 298 K. Needs a function.)   

%Problem children
C=val;  %C is a value derived from appendix 1 equal to g*beta*(rho^2)/(my^2) with unit 1/K*m^3
%Lots of good values can be found as a function of temperature in appendix 1
B=1; %evaporation coefficient

%Horizontal plates (with sides?)
%use horizontal plates
Gr_L_hori=C*(L_hori^3)*(T_film-T_surr);   %where b is fluid coefficient of thermal expansion, g is grav. const., L is the significant length, dT is temp. diff., mu is fluid viscosity
Ra_L_hori=Gr_L_hori*pr_water(T_water);
Nu_L_hori=0.54*Ra_L_hori^(1/4);
h_hori=Nu_L_hori*k_water(T_water)/L_hori; 

%Vertical plates
Gr_L_vert=C*(L_vert^3)*(T_film-T_surr);   %where b is fluid coefficient of thermal expansion, g is grav. const., L is the significant length, dT is temp. diff., mu is fluid viscosity
Ra_L_vert=Gr_L_vert*pr_water(T_water);
Nu_L_vert=((0.825+0.387*(Ra_L_vert^(1/6)))/((1+(0.492/pr_water(T_water))^(9/16))^(8/27)))^2;
h_vert=Nu_L_vert*k_water(T_water)/L_vert;

%Partial pressures
pws=exp((77.345+0.0057*T_surr-7235/T_surr)/(T_surr*8.2));
pw=pws*RH;      %water vapor saturation/partial pressure
%%
%Thermal convection
%Water surface
% dQ_h1=alfa*S_1*(T-T_w)*dt
dQdt_h1=h_hori*A_inner*(T_surr-T_film);                     %[J/s^2] T_w film temp?
%Wall
% dQ_h2=alfa_2*S_1*(T-T_w)*dt
dQdt_h2=h_vert*A_cyl_mantle*(T_surr-T_water);               %[J/s^2]T_w film temp?

%Thermal radiation
% dQ_c=epsilon*sigma*S_1*(T+273.15)^4*dt
dQdt_c=epsilon_glass(1)*sigma*A_cyl*(T_water)^4;            %[J/s^2]

%Evaporation
% dQ_w=beta*S_1*(P_s-P)*dt
dQdt_w=0;%B*(P_s-P);                                        %[J/s^2] Need to find out what B is. 

%Conduction
% dQ_r=-lamda*S_2*(T_w-T)/delta*dt
dQdt_r=-k_glass*(T_water-T_surr)/delta;                     %[J/s^2] delta wall thickness. lamda=k_glass

%Summ of heat loss
% dQ=dQ_h1+dQ_h2+dQ_c+dQ_w+dQ_r
dQdt=dQdt_h1*1+dQdt_h2*1+dQdt_c*1+dQdt_w*1+dQdt_r*1;                  %[J/s^2]

%Relationship between temperature and energy
% dQ=cp*m*T
dQdt=cp_water(T_water)*m_water*(T_water-T_surr);            %[J/s^2]

%Differential equation T(t)
dTdt=dQdt/(cp_water(T_water)*m_water);                       %[K/s]
T_water=T_water-dTdt

end
