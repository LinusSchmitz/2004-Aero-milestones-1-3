% ASEN 2004 Lab 1 Group 013-09
% Maklen Estrada, AJ Lauffer, Fabrizio Roberts, Nate Sanchez,
% Linus Schmitz, and Nathan Tonella

% January 30, 2021

% Finding the coefficient of lift for a 3-D finite wing and comparing it to 
% 2-D coefficient of lift for an airfoil.
% Finding whole air craft and wing drag polar and comparing both values to
% given truth data
% Calculating and comparing various performance values between our
% aerodynamic model and the truth data provided


close all
clear
clc

%% 2-D Airfoil Data
% Note: all data came from the excel file provided on the ASEN 2004 Canvas
% page and was separated into their own files

% read in the airfoil data for both models
Temp = xlsread('Tempest_MH32_Airfoil_Data.xlsx'); % Temp = Tempest UAS
BO = xlsread('BOEING_BACJ_Airfoil_Data.xlsx'); % BO = BOEING 747

% extract coefficient of lift data (2-D)
CLTemp2D = Temp(:,2);
CLBO2D = BO(:,2);

% extract coeifficeint of drag data (2-D)
Cd_Temp = Temp(:,3);
Cd_BO = BO(:,3);

% extract angle of attack for each data point
AoATemp = Temp(:,1);
AoABO = BO(:,1);

%% Hard Coded Values

e = 0.9; % span efficiency factor

% aspect ratio
AR_Temp = 16.5; % found in "Tempest Characteristic"
AR_BO = 7; % found in "BOEING 747-200 Characteristics"

%% Find a0

% caluclate the slope of the 2-D airfoil using equation m=(y1-y2)/(x1-x2)
a0Temp = (CLTemp2D(7)-CLTemp2D(6))/(AoATemp(7)-AoATemp(6));
a0BO = (CLBO2D(7)-CLBO2D(6))/(AoABO(7)-AoABO(6));

%% Find AoA at L=0

% calculate the slope with same form as above using one index above the
% x axis and one index below the x axis
mTemp = (CLTemp2D(4)-CLTemp2D(3))/(AoATemp(4)-AoATemp(3));
mBO = (CLBO2D(5)-CLBO2D(4))/(AoABO(5)-AoABO(4));

% make a large vector of AoA values to use in estimating AoA at L=0
xTemp = linspace(AoATemp(1),AoATemp(end),1000);
xBO = linspace(AoABO(1),AoABO(end),1000);

% use point slope formula to find respective CL values
yTemp = CLTemp2D(4) + mTemp*(xTemp - AoATemp(4));
yBO = CLBO2D(4) + mBO*(xBO - AoABO(4));

% find the index closest to CL = 0
for i = 1:1000
   if yTemp(i) <= 0
       idxTemp = i;
   end
   
   if yBO(i) <= 0
      idxBO = i; 
   end
    
end

% use the index found to assign AoA when L=0
AoATemp_L0 = xTemp(idxTemp);
AoABO_L0 = xBO(idxBO); 

%% Find 3-D Coefficient of Lift
% Note: all equations came directly from the lab document given

% lift curve slope equation (3-D) (eq. 1)
aTemp = a0Temp/(1+((57.3*a0Temp)/(pi*e*AR_Temp)));
aBO = a0BO/(1+((57.3*a0BO)/(pi*e*AR_BO)));

% 3-D coefficent of lift equation (eq. 2)
CLTemp3D = aTemp*(AoATemp-AoATemp_L0);
CLBO3D = aBO*(AoABO-AoABO_L0);

%% Plot 2-D Airfoil vs Approximation 3-D Finite Wing

% plot for the Tempest UAS
figure(1)
hold on
plot(AoATemp,CLTemp2D,'o','Linewidth',1); % AoA vs 2-D airfoil
plot(AoATemp,CLTemp3D,'o','Linewidth',1); % AoA vs 3-D finite wing
xline(0); % y axis
yline(0); % x axis
% labels
title('\alpha vs C_L for Tempest UAS');
ylabel('C_L');
xlabel('\alpha [degrees]');
legend('2-D Airfoil','3-D Finite Wing','Location','southeast');
hold off

% plot for the BOEING 747
figure(2)
hold on
plot(AoABO,CLBO2D,'o','Linewidth',1); % AoA vs 2-D airfoil
plot(AoABO,CLBO3D,'o','Linewidth',1); % AoA vs 3-D finite wing
xline(0); % y axis
yline(0); % x axis
% labels
title('\alpha vs C_L for BOEING 747');
ylabel('C_L');
xlabel('\alpha [degrees]');
legend('2-D Airfoil','3-D Finite Wing','Location','southeast');
hold off

%% Find Wing Drag Polar and Coefficient of Lift at Minimum Drag

% 3-D wing drag polar (eq. 3)
CDwing_Temp = Cd_Temp + ((CLTemp3D.^2)./(pi*e*AR_Temp));
CDwing_BO = Cd_BO + ((CLBO3D.^2)./(pi*e*AR_BO));

% finding the index where wing drag polar is minimizes
[~,idxMinTemp] = min(CDwing_Temp);
[~,idxMinBO] = min(CDwing_BO);

% applying said index to find coefficient of lift at minimum drag
CLminD_Temp = CLTemp3D(idxMinTemp);
CLminD_BO = CLBO3D(idxMinBO);

% plot of wing drag polar vs AoA for Tempest UAS
figure(3)
hold on
plot(AoATemp,CDwing_Temp,'o','Linewidth',1);
xline(0); % y axis
yline(0); % x axis
% labels
title('\alpha vs C_D_W_i_n_g for Tempest UAS');
ylabel('C_D_W_i_n_g');
xlabel('\alpha [degrees]');
hold off

% plot of wing drag polar vs AoA for BOEING 747
figure(4)
hold on
plot(AoABO,CDwing_BO,'o','Linewidth',1);
xline(0); % y axis
yline(0); % x axis
% labels
title('\alpha vs C_D_W_i_n_g for BOEING 747');
ylabel('C_D_W_i_n_g');
xlabel('\alpha [degrees]');
hold off

% displaying coefficent of lift at minimum drag
%fprintf("Tempest: %f \nBOEING: %f \n",CLminD_Temp,CLminD_BO);

%% Drag vs Lift Plots
% BOEING 747 data

% read in truth data for both models
TempTruth = xlsread('Tempest_CFD_Drag_Polar.xlsx');
BOTruth = xlsread('Boeing_747_Drag_Polar_(Exp).xlsx');

% extract coefficient of lift (3-D)
CLtruthTemp = TempTruth(:,2);
CLtruthBO = BOTruth(:,1);

% extract coefficient of drag (3-D)
CDtruthTemp = TempTruth(:,3);
CDtruthBO = BOTruth(:,2);

% Hard coded Values
Cfe_Temp = 0.003; % equivalent skin friction coefficient
Swet_Temp = 2.027; % [m^2] wetted area (Tempest)
Sref_Temp = 0.63; % [m^2] projected area of the wing (Tempest)

Cfe_BO = 0.003; % equivalent skin friction coefficient
Swet_BO = 30000; % [m^2] wetted area (BOEING)
Sref_BO = 5500; % [m^2] projected area of the wing (BOEING

% minimum drag (eq. 10)
CDminTemp = Cfe_Temp*Swet_Temp/Sref_Temp;
CDminBO = Cfe_BO*Swet_BO/Sref_BO;

% Oswald's efficiency factor (eq. 12)
e0_Temp = 1.78*(1-0.045*(AR_Temp^0.68))-0.64; % method 1
e0_BO = 1.78*(1-0.045*(AR_BO^0.68))-0.64;

% method 2 for finding Oswald's Efficiency Factor
Q = 1.05;
P = 0.007;

e0_Temp_Method2 = 1/(Q+P*pi*AR_Temp);

% arbitrary constant (eq. 5)
k1_Temp = 1/(pi*e0_Temp*AR_Temp);
k1_BO = 1/(pi*e0_BO*AR_BO);

k1_Temp_Method2 = 1/(pi*e0_Temp_Method2*AR_Temp);

% whole aircraft drag (eq. 6)
CD_Temp = CDminTemp + k1_Temp*(CLTemp3D - CLminD_Temp).^2;
CD_BO = CDminBO + k1_BO*(CLBO3D - CLminD_BO).^2;

CD_Temp_Method2 = CDminTemp + k1_Temp_Method2*(CLTemp3D - CLminD_Temp).^2;

% plot for Tempest UAS Method 1
figure(5)
hold on
plot(CLTemp3D,CDwing_Temp,'ro','LineWidth',1); % 3-D finite wing drag polar
plot(CLTemp3D,CD_Temp,'bo','LineWidth',1); % whole aircraft drag polar
plot(CLtruthTemp,CDtruthTemp,'mo','LineWidth',1); % truth data drag polar provided
xline(0); % y axis
yline(0); % x axis
% labels
title('C_L vs C_D for Tempest UAS (Method 1)');
ylabel('C_D');
xlabel('C_L');
legend('Wing drag polar','Whole aircraft drag polar','Truth Data');
hold off

% plot for Tempest UAS Method 2
figure(6)
hold on
plot(CLTemp3D,CDwing_Temp,'ro','LineWidth',1); % 3-D finite wing drag polar
plot(CLTemp3D,CD_Temp_Method2,'bo','LineWidth',1); % whole aircraft drag polar
plot(CLtruthTemp,CDtruthTemp,'mo','LineWidth',1); % truth data drag polar provided
xline(0); % y axis
yline(0); % x axis
% labels
title('C_L vs C_D for Tempest UAS (Method 2)');
ylabel('C_D');
xlabel('C_L');
legend('Wing drag polar','Whole aircraft drag polar','Truth Data');
hold off

% plot for BOEING 747
figure(7)
hold on
plot(CLBO3D,CDwing_BO,'ro','LineWidth',1); % 3-D finite wing drag polar
plot(CLBO3D,CD_BO,'bo','LineWidth',1); % whole aircraft drag polar
plot(CLtruthBO,CDtruthBO,'mo','LineWidth',1); % truth data drag polar provided
xline(0); % y axis
yline(0); % x axis
% labels
title('C_L vs C_D for BOEING 747');
ylabel('C_D');
xlabel('C_L');
legend('Wing drag polar','Whole aircraft drag polar','Truth Data');
hold off

%% Comparisons

%Glide Range
% Solve Cl from drag
[ClCdBOMAX,iGRBO] = max(CLBO3D./CD_BO);
[ClCdTempMAX,iGRTemp] = max(CLTemp3D./CD_Temp);

[ClCdBOMAXTru,iGRBOTru] = max(CLtruthBO./CDtruthBO);
[ClCdTempMAXTru,iGRTempTru] = max(CLtruthTemp./CDtruthTemp);

hTemp = 1500; %m
hBO = 10668; %m

wBO = 3705368; %N
%wBO = 3559000; %N
wTemp = 62.78; %N 

pinfBO = 0.38034956; %kg/m^3
%pinfBO = 0.38857; %kg/m^3
pinfTemp = 1.0581; %kg/m^3

VinfBO = sqrt((2*wBO)/(CLBO3D(iGRBO)*Sref_BO*pinfBO));
VinfTemp = sqrt((2*wTemp)/(CLTemp3D(iGRTemp)*Sref_Temp*pinfTemp));

VinfBOTru = sqrt((2*wBO)/(CLtruthBO(iGRBOTru)*Sref_BO*pinfBO));
VinfTempTru = sqrt((2*wTemp)/(CLtruthTemp(iGRTempTru)*Sref_Temp*pinfTemp));

MaxGlideRange = [VinfBO VinfBOTru; VinfTemp VinfTempTru];

%Powered Endurance Propeller (Tempest)
%   largest possible propeller efficieency
%   lowest sfc
%   largest weight fraction of fuel
%   largest L^(3/2)/D
%   largest pinf (sea level)

[~,iETemp] = max((CLTemp3D.^(1.5))./CD_Temp);
[~,iETempTru] = max((CLtruthTemp.^(1.5))./CDtruthTemp);

VinfTemp = sqrt((2*wTemp)/(CLTemp3D(iETemp)*Sref_Temp*pinfTemp));
VinfTempTru = sqrt((2*wTemp)/(CLtruthTemp(iETempTru)*Sref_Temp*pinfTemp));

%Powered Endurance Jet (747)
%   lowest tsfc
%   largest weight fraction of fuel
%   largest L/D

VinfBO = sqrt((2*wBO)/(CLBO3D(iGRBO)*Sref_BO*pinfBO));
VinfBOTru = sqrt((2*wBO)/(CLtruthBO(iGRBOTru)*Sref_BO*pinfBO));

MaxPoweredEndurance = [VinfBO VinfBOTru; VinfTemp VinfTempTru];

%Powered Range Propeller (Tempest)
%   largest possible propeller efficieency
%   lowest sfc
%   largest weight fraction of fuel
%   largest L/D
[~,iPRTemp] = max(CLTemp3D./CD_Temp);
[~,iPRTempTru] = max(CLtruthTemp./CDtruthTemp);

VinfTemp = sqrt((2*wTemp)/(CLTemp3D(iPRTemp)*Sref_Temp*pinfTemp));
VinfTempTru = sqrt((2*wTemp)/(CLtruthTemp(iPRTempTru)*Sref_Temp*pinfTemp));

%Powered Range Jet (747)
%   lowest sfc
%   largest weight fraction of fuel
%   largest L^(1/2)/D
%   lowest pinf (high altittude)
[~,iPRBO] = max((CLBO3D.^(1/2))./CD_BO);
[~,iPRBOTru] = max((CLtruthBO.^(1/2))./CDtruthBO);
VinfBO = sqrt((2*wBO)/(CLBO3D(iPRBO)*Sref_BO*pinfBO));
VinfBOTru = sqrt((2*wBO)/(CLtruthBO(iPRBOTru)*Sref_BO*pinfBO));

MaxPoweredRange = [VinfBO VinfBOTru; VinfTemp VinfTempTru];




