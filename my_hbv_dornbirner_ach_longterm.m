% Predicts runoff based on the relative saturation of a reservoir
% Necessary Input: Precipitation read from file prec_file 
% Precipitation P
% Observed Discharge from input file

% Unit system that is used
% time in days
% lenghts in mm (more professional to use SI Units)
%variables
% Fluxes
%   Potential ET ETP
%   Real Evapotranspiration ET
%   Direct Runoff  Qd from betastore
%   Percolation G from soil store
%   Qtot Discharge from linear Reservoir
% States
%   Soil moisture SM in Betastore
%   Storage Sres in  linear reservoir

%parameters
% A_catch Catchment area to normalize observed discharge into spec. discharge
% beta store 
%   Maximum Soilmoisture Smax
%   beta value controls runoff productions beta
% linear reservoir 
%   resevoir konstant k_res
%%%%%%%%%%%%%%%%%%%%%%%%%% 
% for the komplete model percolation
%   hydraulic conductivity Ks
%   water content at field capacity FC
% Evaporation
%   water content at permanent wilting point PWP;
%   albedo 

clear all;
close all;
clc;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initialise input, model parameter and initial states
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
prec_file='Data_Dorbirnerach_1996_longterm_no_header.csv'; % name precipitation input 
% Set model parameters and inital state in code in code (more elegant to read them from input data
% files
fak=1;
Smax=fak*0.5*1000; % maximum storage  as product of soil porosity and soil depth in mm
beta=0.56;  % Beta parameter
Kres= 1/30 ;% Inverse of residence time in resvoir in days
spin=3000; % Spin off phase in time steps restricted from calculation GOF 
SM_ini=0.5*Smax; % Fraction of total storage
Sres_ini= 20; % Inital filling of the linear store in mm
A_catch=1.5*71*1e06; % Corrected Catchment area (Dornbirner Ach 71 km2)
FC=0.24*Smax; %Field capacity
%Turc parameters
a=0.0031;
b=209.4;
% Input
b_data=dlmread(prec_file,';');
% Precipitation in mm/h
P=fak*b_data(:,1);

Qobs=b_data(:,2)*(1000*3600)/A_catch; % observed discharge in mm/h
%%%%%%% ### New input variables by Erwin 11.3.
R_glo=b_data(:,3); % Global radiation
T_air=b_data(:,4); % Air Temperatur C
rel_hum=b_data(:,5); % relativ humidity %
v_wind=b_data(:,6); % wind speed m/s
p_air=b_data(:,7); % Air pressure in hpa not neede

ntime=length(b_data(:,1));
time=[1:1:ntime]';% array numbering timestepts  for plotting
dt=1; % time step one hour  

% Intialise variables as vector, same length as input time series
SM=zeros(ntime,1); % soil moisture beta store [mm]
Qd=zeros(ntime,1); % direct runoff from beta store [mm/h]
ETP=zeros(ntime,1); %Potential ET [mm/h]
ET=zeros(ntime,1); % Actual ET [mm/h]
G=zeros(ntime,1); % percolation [mm/h]
Sres=zeros(ntime,1); % Storage in linear reservoir [mm]
Qtot=zeros(ntime,1); % Discharge height from linear reservoir

% Initials state vectors 
SM(1)=SM_ini; % 
Sres(1)=Sres_ini; % Inital filling of the linear store in mm

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Processing and time loop
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
NSE_max = 0;
row = 1; 
for beta = 0.4:0.01:0.8
    col = 1;
    for FC = 0.1*Smax:0.01*Smax:0.5*Smax
        for i=1:ntime %time loop
            % Percolation and ET are set as konstant
            if (SM(i) > FC) 
                G(i)=0.0; %Percolation
            else 
                G(i)=0.0;
            end 
             ET(i)=1*ET_TURC_new(a,b,rel_hum(i),T_air(i),R_glo(i)*3600/10000*1/24,SM(i),FC); % Evaporation (a,b,f,T,R)
            % Compute direct runoff, goes into linear reservoir 
            Qd(i)=betastore(P(i),dt,SM(i),Smax,beta);
             % Update soil moisture using water balance eq.
            SM(i+1)=SM(i)+(P(i)-ET(i)-G(i)-Qd(i))*dt;
            bla=1;
            % Check constrains and repair over shoots
            if  SM(i+1) > Smax
                 Qd(i)=Qd(i)+(SM(i+1)-Smax); %Storage overshoot goes to runoff
                 SM(i+1)=Smax;
             end
            % Negative storage?
             if SM(i+1) < 0.
                ET(i)=ET(i)-(-SM(i+1)) ; % Reduce ET by negtive storage
                if ET(i) < 0.
                    ET(i)=0;
                end
                SM(i+1)=0.;      
             end
            %  Calculate discharge from linear store 
            Qtot(i)=linear_reservoir(Sres(i),dt,Kres);
            %Update reservoir storage using storage water balance
            Sres(i+1)=Sres(i)+(Qd(i)+G(i)-Qtot(i))*dt;
        end
        NSE(row, col)=1-(Qtot(spin:ntime)'-Qobs(spin:ntime)')*(Qtot(spin:ntime)-Qobs(spin:ntime))/((Qobs(spin:ntime)'-mean(Qobs(spin:ntime))*ones(length(Qobs(spin:ntime)),1)')*(Qobs(spin:ntime)-mean(Qobs(spin:ntime))*ones(length(Qobs(spin:ntime)),1)));
        if NSE(row, col) > NSE_max
            NSE_max = NSE(row, col);
            beta_best = beta;
            FC_best = FC;
        end
        col = col + 1;
    end
    row = row + 1;
end
FC = FC_best;
beta = beta_best;
figure(1);
subplot(1, 2, 1);
plot(0.4:0.01:0.8, NSE(:, 20));
xlabel('beta');
ylabel('NSE');
subplot(1, 2, 2);
plot(0.1*Smax:0.01*Smax:0.5*Smax, NSE(20, :));
xlabel('FC');
ylabel('NSE');
figure(2);
surf(0.4:0.01:0.8, 0.1*Smax:0.01*Smax:0.5*Smax, NSE);
xlabel('beta');
ylabel('FC');
zlabel('NSE');
figure(3);

[params_op,NSE_op] = pso( 100, 5, 1.1, 1.1, 0.5, ntime, P, Qobs, R_glo, T_air, rel_hum, v_wind, p_air );

beta = params_op(1);
FC = params_op(2);
k0 = params_op(3);
k1 = params_op(4);
k2 = params_op(5);
kperc = params_op(6);
L = params_op(7);
Sres_ini = params_op(8);
S1_ini = params_op(9);
S2_ini = params_op(10);

% Intialise variables as vector, same length as input time series
SM=zeros(ntime,1); % soil moisture beta store [mm]
Qd=zeros(ntime,1); % direct runoff from beta store [mm/h]
ETP=zeros(ntime,1); %Potential ET [mm/h]
ET=zeros(ntime,1); % Actual ET [mm/h]
G=zeros(ntime,1); % percolation [mm/h]
Sres_nl=zeros(ntime,1); % Storage in nonlinear reservoir [mm]
S1=zeros(ntime,1); % Storage in upper box [mm]
S2=zeros(ntime,1); % Storage in upper box [mm]
Q0=zeros(ntime,1); % Surface runoff
Q1=zeros(ntime,1); % Interflow
Q2=zeros(ntime,1); % Groundwater recharge
Qtot_nl=zeros(ntime,1); % Discharge height from non-linear reservoir

% Initials state vectors 
SM(1)=SM_ini; % 
Sres_nl(1)=Sres_ini; % Inital filling of the linear store in mm
S1(1)=S1_ini; % Inital filling of the upper box in mm
S2(1)=S2_ini; % Inital filling of the lower box in mm

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Processing and time loop
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i=1:ntime %time loop
    % Percolation and ET are set as konstant
    if (SM(i) > FC) 
        G(i)=0.0; %Percolation
    else 
        G(i)=0.0;
    end 
     ET(i)=1*ET_TURC_new(a,b,rel_hum(i),T_air(i),R_glo(i)*3600/10000*1/24,SM(i),FC); % Evaporation (a,b,f,T,R)
    % Compute direct runoff, goes into linear reservoir 
    Qd(i)=betastore(P(i),dt,SM(i),Smax,beta);
     % Update soil moisture using water balance eq.
    SM(i+1)=SM(i)+(P(i)-ET(i)-G(i)-Qd(i))*dt;
    % Check constrains and repair over shoots
    if  SM(i+1) > Smax
         Qd(i)=Qd(i)+(SM(i+1)-Smax); %Storage overshoot goes to runoff
         SM(i+1)=Smax;
     end
    % Negative storage?
     if SM(i+1) < 0.
        ET(i)=ET(i)-(-SM(i+1)) ; % Reduce ET by negtive storage
        if ET(i) < 0.
            ET(i)=0;
        end
        SM(i+1)=0.;      
     end
    %  Calculate discharge from non-linear store
    [Q0(i), Q1(i), Qperc(i), Q2(i)]=runoff(G(i),S1(i),S2(i),dt,k0,k1,kperc,k2,L);
    Qtot_nl(i) = Q0(i) + Q1(i) + Q2(i);
    %Update reservoir storage using storage water balance
    S1(i+1)=S1(i)+(Qd(i)-Q0(i)-Q1(i))*dt;
    S2(i+1)=S2(i)+(Qperc(i)-Q2(i))*dt;
    Sres_nl(i+1)=Sres_nl(i)+(Qd(i)+G(i)-Qtot_nl(i))*dt;
end

% Compute Objective function Sum((Qtot-Qobs)^2);
 O=1-(Qtot(spin:ntime)'-Qobs(spin:ntime)')*(Qtot(spin:ntime)-Qobs(spin:ntime))/((Qobs(spin:ntime)'-mean(Qobs(spin:ntime))*ones(length(Qobs(spin:ntime)),1)')*(Qobs(spin:ntime)-mean(Qobs(spin:ntime))*ones(length(Qobs(spin:ntime)),1)));
 bias = sum(Qtot(spin:ntime)*dt)/sum(Qobs(spin:ntime)*dt);
 bias_nl = sum(Qtot_nl(spin:ntime)*dt)/sum(Qobs(spin:ntime)*dt);

% 
%AFALHJH
% % Plot water balance time series
figure(4); % figure
subplot(3,3,1);
bar(time(1:ntime,1),P(1:ntime,1));
xlabel( 'time [h]','fontsize',16);
ylabel( 'P [mm/h]','fontsize',16);
xlim([1 ntime]); % sets axis limits 
set(gca,'fontsize',16);

subplot(3,3,2); % tells matlab that the figure consists of 3 plots in a column
plot(time(1:ntime,1),ET(1:ntime,1),'r--','linewidth',2); % plots sm against time and set linewidth to 2 ppt
xlabel( 'time [h]','fontsize',16); % xlabel with fontsize 16 ppt
ylabel( 'ET [mm/h]','fontsize',16); % ylabel with fontsize 16 ppt
set(gca,'fontsize',16); % sets fonstsize of the current axis to 16 ppt 
xlim([1 ntime]); % sets axis limits

subplot(3,3,3); % start of second plot
plot(time(1:ntime,1),Qd(1:ntime,1),'b-','linewidth',2);
xlabel( 'time [h]','fontsize',16);
ylabel( 'Qd [mm/h]','fontsize',16);
set(gca,'fontsize',16);
xlim([1 ntime]); % sets axis limits

subplot(3,3,4); % tells matlab that the figure consists of 3 plots in a column
plot(time(1:ntime,1),SM(1:ntime,1),'r-','linewidth',2); % plots sm against time and set linewidth to 2 ppt
xlabel( 'time [h]','fontsize',16); % xlabel with fontsize 16 ppt
ylabel( 'Soil moisture store [mm]','fontsize',16); % ylabel with fontsize 16 ppt
set(gca,'fontsize',16); % sets fonstsize of the current axis to 16 ppt 
axis([1 ntime 0 Smax]); % sets axis limits

subplot(3,3,5); % tells matlab that the figure consists of 3 plots in a column
plot(time(1:ntime,1),Qtot(1:ntime,1),'r-','linewidth',2); % plots sm against time and set linewidth to 2 ppt
hold on;
plot(time(1:ntime,1),Qobs(1:ntime,1),'b-','linewidth',2); % plots sm against time and set linewidth to 2 ppt
title(['NSE=' num2str(O),' bias=' num2str(bias)],'fontsize',16);
xlabel( 'time [h]','fontsize',16); % xlabel with fontsize 16 ppt
ylabel( 'Discharge [mm/h]','fontsize',16); % ylabel with fontsize 16 ppt
set(gca,'fontsize',16); % sets fonstsize of the current axis to 16 ppt 
xlim([1 ntime]); % sets axis limits 

subplot(3,3,6); % tells matlab that the figure consists of 3 plots in a column
plot(time(1:ntime,1),Qtot_nl(1:ntime,1),'r-','linewidth',2); % plots sm against time and set linewidth to 2 ppt
hold on;
plot(time(1:ntime,1),Qobs(1:ntime,1),'b-','linewidth',2); % plots sm against time and set linewidth to 2 ppt
title(['NSE=' num2str(NSE_op), ' bias=' num2str(bias_nl)],'fontsize',16);
xlabel( 'time [h]','fontsize',16); % xlabel with fontsize 16 ppt
ylabel( 'Discharge [mm/h]','fontsize',16); % ylabel with fontsize 16 ppt
set(gca,'fontsize',16); % sets fonstsize of the current axis to 16 ppt 
xlim([1 ntime]); % sets axis limits 

subplot(3,3,7); % tells matlab that the figure consists of 3 plots in a column
plot(time(1:ntime,1),Sres(1:ntime,1),'r-','linewidth',2); % plots sm against time and set linewidth to 2 ppt
xlabel( 'time [h]','fontsize',16); % xlabel with fontsize 16 ppt
ylabel( 'Reservoir storage [mm]','fontsize',16); % ylabel with fontsize 16 ppt
set(gca,'fontsize',16); % sets fonstsize of the current axis to 16 ppt 
xlim([1 ntime]); % sets axis limits 

figure(5); % figure
subplot(3,1,1);
plot(cumsum(P(1:ntime,1))/Smax,cumsum(Qd(1:ntime,1))/sum(P(1:ntime,1)),'b-','linewidth',2); % plots sm against time and set linewidth to 2 ppt
% 
xlabel( 'Sum P/Smax []','fontsize',16);
ylabel( 'Sum Q/P','fontsize',16);
set(gca,'fontsize',16); % sets fonstsize of the current axis to 16 ppt 
xlim([0 sum(P(1:ntime,1))/Smax]); % sets axis limits 
% 
subplot(3,1,2);
plot(cumsum(P(1:ntime,1)),cumsum(Qd(1:ntime,1))/sum(P(1:ntime,1)),'b-','linewidth',2); % plots sm against time and set linewidth to 2 ppt
% 
xlabel( 'Sum P  [mm]','fontsize',16);
ylabel( 'Sum Q/P ','fontsize',16);
set(gca,'fontsize',16); % sets fonstsize of the current axis to 16 ppt 
xlim([0 sum(P(1:ntime,1))]); % sets axis limits 
% 
subplot(3,1,3);
plot(cumsum(P(1:ntime,1))/sum(P(1:ntime,1)),cumsum(Qd(1:ntime,1))/sum(P(1:ntime,1)),'b-','linewidth',2); % plots sm against time and set linewidth to 2 ppt
% 
xlabel( 'Sum P/P  [mm]','fontsize',16);
ylabel( 'Sum Q/P ','fontsize',16);
set(gca,'fontsize',16); % sets fonstsize of the current axis to 16 ppt 
xlim([0 sum(P(1:ntime,1))/sum(P(1:ntime,1))]); % sets axis limits 

index = 1;
for i = 0:0.01:4
    rate_obs(index) = length(find(Qobs>=i)) / ntime * 100;
    rate_sim_nl(index) = length(find(Qtot_nl>=i)) / ntime * 100;
    rate_sim(index) = length(find(Qtot>=i)) / ntime * 100;
    index = index + 1;
end
error = sqrt(1 / length(rate_obs) * (rate_sim - rate_obs) * (rate_sim - rate_obs)');
error_nl = sqrt(1 / length(rate_obs) * (rate_sim_nl - rate_obs) * (rate_sim_nl - rate_obs)');

figure(6);
subplot(2, 1, 1);
plot(rate_obs, [0:0.01:4], 'b', rate_sim, [0:0.01:4], 'r');
hold on;
%plot(rate_obs(10:50:310), [0.1:0.5:3.1], '*', rate_sim(10:50:310), [0.1:0.5:3.1], '*');
xlabel('Flow exceed.probab. [%]');
ylabel('Discharge [mm/h]');
legend('observed', 'simulated');
title(['RMSE=' num2str(error)],'fontsize',16);
set(gca,'YTick',0:0.5:4);
set(gca,'fontsize',16);
subplot(2, 1, 2);
plot(time(1:ntime,1),Qtot(1:ntime,1),'r-','linewidth',2); % plots sm against time and set linewidth to 2 ppt
hold on;
plot(time(1:ntime,1),Qobs(1:ntime,1),'b-','linewidth',2); % plots sm against time and set linewidth to 2 ppt
title(['NSE=' num2str(O), ' bias=' num2str(bias)],'fontsize',16);
xlabel( 'time [h]','fontsize',16); % xlabel with fontsize 16 ppt
ylabel( 'Discharge [mm/h]','fontsize',16); % ylabel with fontsize 16 ppt
set(gca,'fontsize',16); % sets fonstsize of the current axis to 16 ppt 
xlim([1 ntime]); % sets axis limits

figure(7);
subplot(2, 1, 1);
plot(rate_obs, [0:0.01:4], 'b', rate_sim_nl, [0:0.01:4], 'r');
hold on;
%plot(rate_obs(10:50:310), [0.1:0.5:3.1], '*', rate_sim_nl(10:50:310), [0.1:0.5:3.1], '*');
xlabel('Flow exceed.probab. [%]');
ylabel('Discharge [mm/h]');
legend('observed', 'simulated');
title(['RMSE=' num2str(error_nl)],'fontsize',16);
set(gca,'YTick',0:0.5:4);
set(gca,'fontsize',16);
subplot(2, 1, 2);
plot(time(1:ntime,1),Qtot_nl(1:ntime,1),'r-','linewidth',2); % plots sm against time and set linewidth to 2 ppt
hold on;
plot(time(1:ntime,1),Qobs(1:ntime,1),'b-','linewidth',2); % plots sm against time and set linewidth to 2 ppt
title(['NSE=' num2str(NSE_op), ' bias=' num2str(bias_nl)],'fontsize',16);
xlabel( 'time [h]','fontsize',16); % xlabel with fontsize 16 ppt
ylabel( 'Discharge [mm/h]','fontsize',16); % ylabel with fontsize 16 ppt
set(gca,'fontsize',16); % sets fonstsize of the current axis to 16 ppt 
xlim([1 ntime]); % sets axis limits
