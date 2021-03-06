clear, clc;
m_water = 0.2;              %[kg]   Initial mass
T_water = 80+273.15;        %[K]   Initial temperature

%% SIDES

%Temperature profile
%Bulk - Inner wall - Outer wall - Surroundings

F=@T_wall_solve;                                    %Declaring function
x0=[330 320];                                       %Initial guess of temperatures
OPTIONS=optimoptions(@fsolve,'Display','off');      %Turning off diplay
[T_wall,y]=fsolve(F,x0,OPTIONS);                    %Fsolve for temperature profile
T_wall_inner=T_wall(1);                             %Temperature of inner wall
T_wall_outer=T_wall(2);                             %Temperature of outer wall

%Heat profile
%Convection - Conduction - Convection + Radiation


%% TOP


%Temp profile
%Bulk - Water surface - Surroundings

G=@T_surface_solve;                                    %Declaring function
x1=[330];                                       %Initial guess of temperatures
[T_surface,y]=fsolve(G,x1,OPTIONS);                    %Fsolve for temperature profile

%Heat profile
%Convection - Convection + Radiation + Evaporation



