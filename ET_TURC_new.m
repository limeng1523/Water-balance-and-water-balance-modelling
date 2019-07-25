% ET_TURC calculates actual evaporation 
% Input Temperatur T in °C
% Global Radiation R in J/qcm 
% relative humidity in %

function ET=ET_TURC_new(a,b,f,T,R,SM,FC)

if T < 0
    ETP=0.;
else
    if f > 50 
        C_turc =1;
    elseif f <= 50
        C_turc=1+(50-f)/70;
    end
    ETP=C_turc*a*T/(T+15)*(R+b);
end
 if SM < FC
      ET=ETP*SM/FC;
 else
     ET=ETP;
 end
end


    
