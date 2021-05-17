clc
Ra_L_water_surface=2*10^8;
if Ra_L_water_surface<2*10^7 && Ra_L_water_surface>=10^5
        Nu_L_water_surface=0.54*Ra_L_water_surface^(1/4);
        elseif Ra_L_water_surface<=3*10^10 && Ra_L_water_surface>=2*10^7
        Nu_L_water_surface=0.14*Ra_L_water_surface^(1/3);
        else
        disp('Ra_L_water_surface out of bounds')
        Ra_L_water_surface
end