% function out=Data_retrieval(T)


%We want to make a function that uses intial conditions to calculate heat
%loss. Use that heat loss to calculate a new temperature. Bind all of this
%together to an ode.

%% Fetching experimental data



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
epsilon_glass=[0.94, 0.95]; %emissivity, maybe do average instead

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
height_small_water=V_water/(pi*(((beaker_small(4)*(10^(-3)))/2)^2));    %m
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

val=2.035*10^-9;    %taken at 298 K. Needs a function.)   

%Beaker characteristic heights in mm
L=height_small_water;       %height [m]

%D/L>=35/(Gr_L^(1/4));       %where d is the cylinder diameter, L is the boundry layer thickness
C=val;  %C is a value derived from appendix 1 equal to g*beta*(rho^2)/(my^2) with unit 1/K*m^3
%Lots of good values can be found as a function of temperature in appendix 1
Gr_L=C*(L^3)*(T_film-T_surr);   %where b is fluid coefficient of thermal expansion, g is grav. const., L is the significant length, dT is temp. diff., mu is fluid viscosity
Ra_L=Gr_L*Pr_water_var;
Nu_L=((0.825+0.387*(Ra_L^(1/6)))/((1+(0.492/Pr_water_var)^(9/16))^(8/27)))^2;
h=Nu_L*k_water_var/L;
Q_vertical=h*A_small_water*(T_film-T_surr);   %J/s

%Horizontal plates (with sides?)
%use horizontal plates
Nu_L_horizontal=054*Ra_L^(1/4);
h_horizontal=Nu_L_horizontal*k_water_var/L;
Q_horizontal=h_horizontal*A_small_water*(T_film-T_surr); %J/s

Q_convection=Q_horizontal+Q_vertical;       %J/s

%% Radiation
%From sides (cylindrical walls)
%Water surface
%Beaker and water surface
%Assuming grey and opaque surface? -> rho=1-epsilon
%Assume water has same properties as glass.
%glass=_1, surroundings=_2
%beaker only sees surroundings?
%Assumes J_2 is much smaller than J_1 so J_2=0!?

F_12=1;
F_11=0;
F_21=0;
F_22=1;
% 
A_cyl=height_small_water*beaker_small(3)*10^-3*pi;
A_top=(pi*((beaker_small(3)*10^-3)/2)^2);
A_1=A_cyl+A_top;
% A_1*F_12=A_2*F_21;
rho_1=1-epsilon_glass(1);

G_12=A_1*F_12;

%Stefan-Boltzmann Law
sigma=5.676*10^-8;  %W/m^2*K^4
E_b1=sigma*T_water^4;          %Emissive power

J_1=rho_1*G_12+epsilon_glass(1)*E_b1;    %J/s*m^2
J_2=0;    %J/s*m^2
q_1=A_1*F_12*(J_1-J_2) %J/s Rate of heat exchange between 1 and 2

%% Evaporation
A_circle=(pi*((beaker_small(4)*10^-3)/2)^2);
%Antoine Equation

% A=8.140191; B=1810.94; C=244.485;
% pws=(10^(A-(B/(C+T))))*133,322;   %[Pa]

pa=101.3*10^3;  %[Pa]


pws=exp((77.345+0.0057*T_surr-7235/T_surr)/(T_surr*8.2));
pw=pws*RH;      %water vapor saturation/partial pressure

x=0.62198*pw/(pa-pw);
x_s=0.62198*pws/(pa-pws);

v=0;
Theta=(25+19*v);   % Theta = evaporation coefficient (kg/(m2 h)) % v = velocity of air above the water surface (m/s)

gh=Theta*A_circle*(x_s-x);    % gh = amount of evaporated water per hour (kg/h)   % A = water surface area (m2)
gs=gh*(3600/1000)       %gs = g/s
%This is an under estimate based on the experimental data.

%% Conduction
%Conduction from water, through glass, to air
%The bottom is probably assumed to be 0.

%% Plot of mass and temp
figure
subplot(2,1,1)
title('Small beaker: Mass');  xlabel('Time (s)');  ylabel('Mass (g)');     legend;
hold on
plot(small_data_1.data(1:61,1),small_data_1.data(1:61,3))
plot(small_data_2.data(1:61,1),small_data_2.data(1:61,3))
plot(small_data_3.data(1:61,1),small_data_3.data(1:61,3))
plot(small_data_4.data(1:61,1),small_data_4.data(1:61,3))
plot(small_data_5.data(1:61,1),small_data_5.data(1:61,3))
hold off

subplot(2,1,2)
title('Small beaker: Temperature (C°)');  xlabel('Time (s)');  ylabel('Temperature (C°)');     legend;
hold on
plot(small_data_1.data(1:61,1),small_data_1.data(1:61,2))
plot(small_data_2.data(1:61,1),small_data_2.data(1:61,2))
plot(small_data_3.data(1:61,1),small_data_3.data(1:61,2))
plot(small_data_4.data(1:61,1),small_data_4.data(1:61,2))
plot(small_data_5.data(1:61,1),small_data_5.data(1:61,2))
hold off

figure
subplot(2,1,1)
title('Large beaker: Mass');  xlabel('Time (s)');  ylabel('Mass (g)');     legend;
hold on
plot(large_data_1.data(1:61,1),large_data_1.data(1:61,3))
plot(large_data_2.data(1:61,1),large_data_2.data(1:61,3))
plot(large_data_3.data(1:61,1),large_data_3.data(1:61,3))
hold off

subplot(2,1,2)
title('Large beaker: Temperature (C°)');  xlabel('Time (s)');  ylabel('Temperature (C°)');     legend;
hold on
plot(large_data_1.data(1:61,1),large_data_1.data(1:61,2))
plot(large_data_2.data(1:61,1),large_data_2.data(1:61,2))
plot(large_data_3.data(1:61,1),large_data_3.data(1:61,2))
hold off

% Evaporation % THIS IS BAD
%Water surface
m_flow_small_1_matrix=zeros(60,1);
m_flow_small_2_matrix=zeros(60,1);
m_flow_small_3_matrix=zeros(60,1);
m_flow_small_4_matrix=zeros(60,1);
m_flow_small_5_matrix=zeros(60,1);

m_flow_large_1_matrix=zeros(60,1);
m_flow_large_2_matrix=zeros(60,1);
m_flow_large_3_matrix=zeros(60,1);

i=0
for i=1:60
    m_flow_1=-(small_data_1.data(i+1,3)-small_data_1.data(i,3))/(small_data_1.data(i+1,1)-small_data_1.data(i,1));
    m_flow_small_1_matrix(i)=m_flow_1;
    
    m_flow_2=-(small_data_2.data(i+1,3)-small_data_2.data(i,3))/(small_data_2.data(i+1,1)-small_data_2.data(i,1));
    m_flow_small_2_matrix(i)=m_flow_2;
    
    m_flow_3=-(small_data_3.data(i+1,3)-small_data_3.data(i,3))/(small_data_3.data(i+1,1)-small_data_3.data(i,1));
    m_flow_small_3_matrix(i)=m_flow_3;
    
    m_flow_4=-(small_data_4.data(i+1,3)-small_data_4.data(i,3))/(small_data_4.data(i+1,1)-small_data_4.data(i,1));
    m_flow_small_4_matrix(i)=m_flow_4;
    
    m_flow_5=-(small_data_5.data(i+1,3)-small_data_5.data(i,3))/(small_data_5.data(i+1,1)-small_data_5.data(i,1));
    m_flow_small_5_matrix(i)=m_flow_5;
    
    
    m_flow_1_large=-(large_data_1.data(i+1,3)-large_data_1.data(i,3))/(large_data_1.data(i+1,1)-large_data_1.data(i,1));
    m_flow_large_1_matrix(i)=m_flow_1_large;
    
    m_flow_2_large=-(large_data_2.data(i+1,3)-large_data_2.data(i,3))/(large_data_2.data(i+1,1)-large_data_2.data(i,1));
    m_flow_large_2_matrix(i)=m_flow_2_large;
    
    m_flow_3_large=-(large_data_3.data(i+1,3)-large_data_3.data(i,3))/(large_data_3.data(i+1,1)-large_data_3.data(i,1));
    m_flow_large_3_matrix(i)=m_flow_3_large;
    
    % m_1=small_data_1.data((1:61),3)
    % t_1=small_data_1.data((1:61),2)
    %
    % m_2=small_data_2.data((1:61),3)
    % t_2=small_data_2.data((1:61),2)
end


figure
subplot(2,1,1)
grid on
title('Small beaker: Mass flow');  xlabel('Temperature (C°)');  ylabel('Mass flow (g/s)');     legend;
hold on
plot(small_data_1.data(1:60,2),m_flow_small_1_matrix)
plot(small_data_2.data(1:60,2),m_flow_small_2_matrix)
plot(small_data_3.data(1:60,2),m_flow_small_3_matrix)
plot(small_data_4.data(1:60,2),m_flow_small_4_matrix)
plot(small_data_5.data(1:60,2),m_flow_small_5_matrix)
hold off

subplot(2,1,2)
title('Large beaker: Mass flow');  xlabel('Temperature (C°)');  ylabel('Mass flow (g/s)');     legend;
hold on
plot(small_data_1.data(1:60,2),m_flow_small_1_matrix)
plot(small_data_2.data(1:60,2),m_flow_small_2_matrix)
plot(small_data_3.data(1:60,2),m_flow_small_3_matrix)
hold off


m_flow_average=((m_flow_small_1_matrix+m_flow_small_2_matrix+m_flow_small_3_matrix+m_flow_small_4_matrix+m_flow_small_5_matrix)./5).*10.^-3;
m_flow_average_large=((m_flow_large_1_matrix+m_flow_large_2_matrix+m_flow_large_3_matrix)./3).*10.^-3;

T_average_small=(small_data_1.data(1:60,2)+small_data_2.data(1:60,2)+small_data_3.data(1:60,2)+small_data_4.data(1:60,2)+small_data_5.data(1:60,2))./5;
T_average_large=(large_data_1.data(1:60,2)+large_data_2.data(1:60,2)+large_data_3.data(1:60,2))./3;

T_film_average_small=((T_average_small+273.15+T_surr)./2);      %K
T_film_average_large=((T_average_large+273.15+T_surr)./2);      %K

Evap_energy_system_small=cp_water(T_film_average_small).*T_film_average_small.*m_flow_average;          %J/s
Evap_energy_system_large=cp_water(T_film_average_large).*T_film_average_large.*m_flow_average_large;    %J/s

% end