%% Fetching experimental data
clear, clc;


addpath('C:\Chalmers\Year 3\Mathematical modeling\Project');


% for i=1:5
% num=i;
% number=num2str(num);
% small_data=importdata(['small_exp', number, '.txt']);
% end
% 
% for i=1:3
% num=i;
% number=num2str(num);
% data_large=importdata(['large_exp', number, '.txt']);
% end
    
small_data_1=importdata('small_exp1.txt');
small_data_2=importdata('small_exp2.txt');
small_data_3=importdata('small_exp3.txt');
small_data_4=importdata('small_exp4.txt');
small_data_5=importdata('small_exp5.txt');

large_data_1=importdata('large_exp1.txt');
large_data_2=importdata('large_exp2.txt');
large_data_3=importdata('large_exp3.txt');

%Data for glass beakers are in ml, mm, mm, mm 
%and are size, height, outer diameter, and inner diameter in order
beaker_large=[800,135,98,93];   
beaker_small=[250,95,70,67];

%Data for glass material
rho_glass=2500;             %kg/m^3
cp_glass=840;               %J/kg*K
k_glass=0.9;                %W/m*K
epsilon_glass=[0.94, 0.95]; %maybe do average instead

%Additional data
RH=0.3;                     %relative humidity
T_surr=20+273.15;           %K ambient temperature
g=9.81;                     %m/s^2

%Calculations with small_data_1
%Will be changed to for loop for all temp and beakers
T_water=small_data_1.data(1,2)+273.15;      %K
m_water=small_data_1.data(1,3)*10^-3;       %kg
rho_water=rho_water(T_water);               %kg/m^3
V_water=m_water/rho_water;                  %m^3
height_water=V_water/(pi*(((beaker_small(4)*(10^(-3)))/2)^2));
A_water=height_water*beaker_small(3)*pi
%% Convection
%Vertical wall (free convection vertical cylindar
%All properties taken at film temperature which is (T_s-T_infinity)/2
T=273.15;   %Temperature for call functions may be at film temperature
%Call to functions fetching variables dependent on temperature
p_water=p_water(T); Pr_water=pr_water(T);   dHvap_water=dHvap_water(T);   
cp_water=cp_water(T); my_water=my_water(T); k_water=k_water(T); rho_water=rho_water(T);
cp_air=cp_air(T);   my_air=my_air(T);   k_air=k_air(T); rho_air=rho_air(T);

val=2.035*10^-9;    %taken at 298 K. Needs a function to interpolate   

%Beaker characteristic heights in mm
L_large=beaker_large(4);    
L_small=height_water+(beaker_small(3)-beaker_small(4))*10^-3;   %Change since we need to include convection through glass

%D/L>=35/(Gr_L^(1/4));       %where d is the cylinder diameter, L is the boundry layer thickness
C=val;  %C is a value derived from appendix 1 equal to g*beta*(rho^2)/(my^2) with unit 1/K*m^3
%Lots of good values can be found as a function of temperature in appendix 1
Gr_L=C*(L_small^3)*dT;   %where b is fluid coefficient of thermal expansion, g is grav. const., L is the significant length, dT is temp. diff., mu is fluid viscosity
Ra_L=Gr_L*Pr_water;
Nu_L=((0.825+0.387*(Ra_L^(1/6)))/((1+(0.492/Pr_water)^(9/16))^(8/27)))^2;
h=Nu_L*k_water/L_small;
Q=h*A*(T-T_w);
%Horizontal plates (with sides?)

%% Radiation
%From sides (cylindrical walls)
%Water surface
%Beaker and water surface

%% Evaporation
%Water surface

%% Conduction
%Conduction from water, through glass, to air
%The bottom is probably assumed to be 0.





