function dTdt=Good_version(T_water)
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
delta=d_outer-d_inner;                      %[m]
m_glass_mantle=A_cyl_mantle*delta*rho_glass;%[kg]    
 
val=2.035*10^-9;    %taken at 298 K. Needs a function.)   

%Problem children
C=val;  %C is a value derived from appendix 1 equal to g*beta*(rho^2)/(my^2) with unit 1/K*m^3


%fsolve to find the inside and outside wall temperatures
F=@T_wall_solve;
x0=[330 320];
OPTIONS=optimoptions(@fsolve,'Display','off');
[T_wall,y]=fsolve(F,x0,OPTIONS);
T_wall_inner=T_wall(1);
T_wall_outer=T_wall(2);

G=@T_surface_solve;                                     %Declaring function
x1=280;                                                %Initial guess of temperatures
[T_surface,y]=fsolve(G,x1,OPTIONS);                     %Fsolve for temperature profile
T_surface=real(T_surface);                             %NOT RIGHT!!! FIX!!!
% T_film=(T_surr+T_water)/2;   %Temperature for call functions may be at film temperature




%% Convection 

%Convection from water to inner wall 
Gr_L_water_wall=C*(L_vert^3)*(T_water-T_wall_inner);   %where b is fluid coefficient of thermal expansion, g is grav. const., L is the significant length, dT is temp. diff., mu is fluid viscosity
Ra_L_water_wall=Gr_L_water_wall*pr_water(T_water);
Nu_L_water_wall=((0.825+0.387*(Ra_L_water_wall^(1/6)))/((1+(0.492/pr_water(T_water))^(9/16))^(8/27)))^2;
h_water_wall=Nu_L_water_wall*k_water(T_water)/L_vert;            %heat transfer coefficient for water-wall_inner

%Calc conv water-inner wall
% q=m_glass_mantle*cp*dT;

val=2.035*10^-9;    %taken at 298 K. Needs a function.)   

%Problem children
C=val;  %C is a value derived from appendix 1 equal to g*beta*(rho^2)/(my^2) with unit 1/K*m^3
%Lots of good values can be found as a function of temperature in appendix 1
B=1; %evaporation coefficient

%Horizontal plates (with sides?)
%use horizontal plates
Gr_L_hori=C*(A_inner^3)*(T_surface-T_surr);   %where b is fluid coefficient of thermal expansion, g is grav. const., L is the significant length, dT is temp. diff., mu is fluid viscosity
Ra_L_hori=Gr_L_hori*pr_water(T_surface);
Nu_L_hori=0.54*Ra_L_hori^(1/4);
h_surface=Nu_L_hori*k_water(T_surface)/A_inner; 

%Partial pressures
pws=exp((77.345+0.0057*T_surr-7235/T_surr)/(T_surr*8.2));
pw=pws*RH;      %water vapor saturation/partial pressure
%%
% %Thermal convection
% %Water surface
% % dQ_h1=alfa*S_1*(T-T_w)*dt
% dQdt_water_air=h_hori*A_inner*(T_water-T_surr);                     %[J/s^2] T_w film temp?
% %Water to inner wall
% % dQ_h2=alfa_2*S_1*(T-T_w)*dt
% dQdt_water_wall=h_water_wall*A_cyl_mantle*(T_water-T_wall_inner);               %[J/s^2]T_w film temp?
% %Outer wall to air
% % dQ_h2=alfa_2*S_1*(T-T_w)*dt
% dQdt_wall_air=h_wall_air*A_cyl_mantle*(T_wall_outer-T_surr);               %[J/s^2]T_w film temp?
% 
% %Thermal radiation
% % dQ_c=epsilon*sigma*S_1*(T+273.15)^4*dt
% dQdt_c=epsilon_glass(1)*sigma*A_cyl*(T_water)^4;            %[J/s^2]
% 
% %Evaporation
% % dQ_w=beta*S_1*(P_s-P)*dt
% dQdt_w=0;%B*(P_s-P);                                        %[J/s^2] Need to find out what B is. 
% 
% %Conduction
% % dQ_r=-lamda*S_2*(T_w-T)/delta*dt
% dQdt_r=0;%k_glass*(T_water-T_surr)/delta;                     %[J/s^2] delta wall thickness. lamda=k_glass
% 
% %Summ of heat loss
% % dQ=dQ_h1+dQ_h2+dQ_c+dQ_w+dQ_r
% dQdt=dQdt_water_air*1+dQdt_wall_air*1+dQdt_c*1+dQdt_w*1+dQdt_r*1;                  %[J/s^2]
% 
% %Relationship between temperature and energy
% % dQ=cp*m*T
% %dQdt=cp_water(T_water)*m_water*(T_water-T_surr);            %[J/s^2]

dQdt_water_surface=h_surface*A_inner*(T_water-T_surface);
dQdt_water_wall=h_water_wall*A_cyl*(T_water-T_wall_inner);

dQdt=dQdt_water_surface+dQdt_water_wall;

%Differential equation T(t)
dTdt=-dQdt/(cp_water(T_water)*m_water);                       %[K/s]
end
