function O = hbv_nl( ntime, P, Qobs, R_glo, T_air, rel_hum, v_wind, p_air, params_pso )

a=0.0031;
b=209.4;
Smax = 500; % maximum storage  as product of soil porosity and soil depth in mm
spin=3000; % Spin off phase in time steps restricted from calculation GOF 
SM_ini=250; % Fraction of total storage
SM=zeros(ntime,1); % soil moisture beta store [mm]
Qd=zeros(ntime,1); % direct runoff from beta store [mm/h]
ETP=zeros(ntime,1); %Potential ET [mm/h]
ET=zeros(ntime,1); % Actual ET [mm/h]
Sres_nl=zeros(ntime,1); % Storage in nonlinear reservoir [mm]
S1=zeros(ntime,1); % Storage in upper box [mm]
S2=zeros(ntime,1); % Storage in upper box [mm]
Q0=zeros(ntime,1); % Surface runoff
Q1=zeros(ntime,1); % Interflow
Q2=zeros(ntime,1); % Groundwater recharge
Qtot_nl=zeros(ntime,1); % Discharge height from non-linear reservoir

SM(1) = SM_ini;
Sres_nl(1) = params_pso(8);
S1(1) = params_pso(9); % Inital filling of the upper box in mm
S2(1) = params_pso(10); % Inital filling of the lower box in mm
for i=1:ntime %time loop
    % Percolation and ET are set as konstant

     ET(i) = ET_TURC_new(a,b,rel_hum(i),T_air(i),R_glo(i)*3600/10000*1/24,SM(i),params_pso(2)); % Evaporation (a,b,f,T,R)
    % Compute direct runoff, goes into linear reservoir 
    Qd(i)=betastore(P(i),1,SM(i),Smax,params_pso(1));
     % Update soil moisture using water balance eq.
    SM(i+1)=SM(i)+(P(i)-ET(i)-Qd(i));
    % Check constrains and repair over shoots
    if  SM(i+1) > Smax
         Qd(i)=Qd(i)+(SM(i+1)-Smax); %Storage overshoot goes to runoff
         SM(i+1)=Smax;
     end
    % Negative storage?
     if SM(i+1) < 0.
        ET(i)=ET(i) + SM(i+1) ; % Reduce ET by negtive storage
        if ET(i) < 0.
            ET(i)=0;
        end
        SM(i+1)=0.;      
     end
    %  Calculate discharge from non-linear store
    [Q0(i), Q1(i), Qperc(i), Q2(i)]=runoff(0,S1(i),S2(i),1,params_pso(3),params_pso(4),params_pso(5),params_pso(6),params_pso(7));
    Qtot_nl(i) = Q0(i) + Q1(i) + Q2(i);
    %Update reservoir storage using storage water balance
    S1(i+1)=S1(i)+(Qd(i)-Q0(i)-Q1(i));
    S2(i+1)=S2(i)+(Qperc(i)-Q2(i));
    Sres_nl(i+1)=Sres_nl(i)+(Qd(i)-Qtot_nl(i));
end
O = 1-(Qtot_nl(spin:ntime)'-Qobs(spin:ntime)')*(Qtot_nl(spin:ntime)-Qobs(spin:ntime))/((Qobs(spin:ntime)'-mean(Qobs(spin:ntime))*ones(length(Qobs(spin:ntime)),1)')*(Qobs(spin:ntime)-mean(Qobs(spin:ntime))*ones(length(Qobs(spin:ntime)),1)));

end

