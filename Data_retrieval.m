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

%% Convection
%Vertical wall (free convection vertical cylindar
%D/L>=35/(Gr_L^(1/4));       %where d is the cylinder diameter and L is the boundry layer thickness.

%Horizontal plates (with sides?)

%% Radiation
%From sides (cylindrical walls)
%Water surface
%Beaker and water surface

%% Evaporation
%Water surface

%% Conduction
%Probably assumed to be 0.





