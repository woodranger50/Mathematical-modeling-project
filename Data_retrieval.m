%% Fetching experimental data
clf, clear, clc;


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
rho_water_var=rho_water(T_water);               %kg/m^3
V_water=m_water/rho_water_var;                  %m^3
height_small_water=V_water/(pi*(((beaker_small(4)*(10^(-3)))/2)^2));
A_small_water=height_small_water*beaker_small(3)*pi;
%% Convection
%Vertical wall (free convection vertical cylindar
%All properties taken at film temperature which is (T_s-T_infinity)/2
T_film=(T_surr+T_water)/2;   %Temperature for call functions may be at film temperature
%Call to functions fetching variables dependent on temperature
p_water_var=p_water(T_film);        Pr_water_var=pr_water(T_film);      dHvap_water_var=dHvap_water(T_film);   
cp_water_var=cp_water(T_film);      my_water_var=my_water(T_film);      k_water_var=k_water(T_film);        
rho_water_var=rho_water(T_film);    cp_air_Var=cp_air(T_film);      
my_air_var=my_air(T_film);          k_air_var=k_air(T_film); rho_air_var=rho_air(T_film);

val=2.035*10^-9;    %taken at 298 K. Needs a function to interpolate (g*beta*(rho^2)/(my^2) with unit 1/K*m^3)   

%Beaker characteristic heights in mm
L_large=beaker_large(4);            %Get working for large later
L_small=height_small_water;

%D/L>=35/(Gr_L^(1/4));       %where d is the cylinder diameter, L is the boundry layer thickness
C=val;  %C is a value derived from appendix 1 equal to g*beta*(rho^2)/(my^2) with unit 1/K*m^3
%Lots of good values can be found as a function of temperature in appendix 1
Gr_L=C*(L_small^3)*(T_film-T_surr);   %where b is fluid coefficient of thermal expansion, g is grav. const., L is the significant length, dT is temp. diff., mu is fluid viscosity
Ra_L=Gr_L*Pr_water_var;
Nu_L=((0.825+0.387*(Ra_L^(1/6)))/((1+(0.492/Pr_water_var)^(9/16))^(8/27)))^2;
h=Nu_L*k_water_var/L_small;
Q_small=h*A_small_water*(T_film-T_surr)   %KJ

%Horizontal plates (with sides?)
%not sure if one can use horizontal plates

%% Radiation
%From sides (cylindrical walls)
%Water surface
%Beaker and water surface

%% Evaporation
%Water surface
m_flow_matrix=zeros(1,60);
i=0
for i=1:60
    m_flow_1=(small_data_1.data(i+1,3)-small_data_1.data(i,3))/(small_data_1.data(i+1,1)-small_data_1.data(i,1));
    m_flow_matrix(i)=m_flow_1;
    temp_matrix(i)=small_data_1.data(i,1);      %finish creating matrix and plot T vs m_flow
% m_1=small_data_1.data((1:61),3)
% t_1=small_data_1.data((1:61),2)
% 
% m_2=small_data_2.data((1:61),3)
% t_2=small_data_2.data((1:61),2)
end
% figure
% subplot(2,1,1)
% plot(m_1,t_1)
% subplot(2,1,2)
% plot(m_2,t_2)
%% Conduction
%Conduction from water, through glass, to air
%The bottom is probably assumed to be 0.





