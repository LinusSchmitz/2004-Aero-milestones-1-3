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


%%% INDIVIDUAL PORTION %%%

close all
clear
clc

%% 2-D Airfoil Data
% Note: all data came from the excel file AJ created using the graphs given
% on the second milestone document

% read in the airfoil data for both models
FlatPlate = xlsread('2D_Flat_Plate.xlsx'); % my best estimate of data for flat plate

% extract coefficient of lift data (2-D)
CL2D = FlatPlate(:,2);

% extract coeifficeint of drag data (2-D)
Cd = FlatPlate(:,3);

% extract angle of attack for each data point
AoA = FlatPlate(:,1);

%% Hard Coded Values

e = 0.9; % span efficiency factor
c = 0.23; %chord [m]
AR = 4;
% aspect ratio
b = c*AR;
%% Find a0

% caluclate the slope of the 2-D airfoil using equation m=(y1-y2)/(x1-x2)
a0 = (CL2D(7)-CL2D(6))/(AoA(7)-AoA(6));

%% Find AoA at L=0

% NOTE: below code is unnecessary unless you did not find a value for AOA
% where the cl was 0, then you have to interpolate it

% % calculate the slope with same form as above using one index above the
% % x axis and one index below the x axis
% m = (CL2D(NaN)-CL(NaN))/(AoA(NaN)-AoATemp(NaN));
% 
% % make a large vector of AoA values to use in estimating AoA at L=0
% x = linspace(AoA(1),AoA(end),1000);
% 
% % use point slope formula to find respective CL values
% y = CL2D(4) + m*(x - AoA(4));
% 
% % find the index closest to CL = 0
% for i = 1:1000
%    if y(i) <= 0
%        idx = i;
%    end     
% end
%
% % use the index found to assign AoA when L=0
% AoA_L0 = x(idx);
 

% if there was an index where Cl = 0 for some AoA value, use that index
idx = CL2D == 0;
AoA_L0 = AoA(idx);


%% Find 3-D Coefficient of Lift
% Note: all equations came directly from the lab document given

% lift curve slope equation (3-D) (eq. 1)
a = a0/(1+((57.3*a0)/(pi*e*AR)));

% 3-D coefficent of lift equation (eq. 2)
CL3D = a.*(AoA-AoA_L0);

CL3D(end) = NaN; % stop the data where it would no longer be linear

%% Plot 2-D Airfoil vs Approximation 3-D Finite Wing

figure(1)
hold on
plot(AoA,CL2D,'o','Linewidth',1); % AoA vs 2-D airfoil
plot(AoA,CL3D,'o','Linewidth',1); % AoA vs 3-D finite wing
xline(0); % y axis
yline(0); % x axis
% labels
title('\alpha vs C_L');
ylabel('C_L');
xlabel('\alpha [degrees]');
legend('2-D Airfoil','3-D Finite Wing','Location','southeast');
hold off

%% Find Wing Drag Polar and Coefficient of Lift at Minimum Drag

% 3-D wing drag polar (eq. 3)
CDwing = Cd + ((CL3D.^2)./(pi*e*AR));

% finding the index where wing drag polar is minimizes
[~,idxMin] = min(CDwing);
[~,idxCl] = min(CL3D);

% applying said index to find coefficient of lift at minimum drag
CLminD = CL3D(idxMin);

% plot of wing drag polar vs AoA
figure(3)
hold on
plot(AoA,CDwing,'o','Linewidth',1);
xline(0); % y axis
yline(0); % x axis
% labels
title('\alpha vs C_D_W_i_n_g');
ylabel('C_D_W_i_n_g');
xlabel('\alpha [degrees]');
hold off

%% Drag vs Lift Plots

% Hard coded Values
Cfe = 0.003; % equivalent skin friction coefficient
Sref = c*b; % [m^2] projected area of the wing
sWetFuselage = 0.534/2;
sWetNoTail = sWetFuselage + Sref * 2 + 0.014;
Stail = 0.06; % [m^2] tail wetted area
Swet = sWetNoTail + Stail;
Sother = sWetFuselage + 0.014 + Stail;

% minimum drag (eq. 10)
CD0 = Cfe*Swet/Sref;

% Oswald's efficiency factor (eq. 12) [different value?]
%e0 = 1.78*(1-0.045*(AR^0.68))-0.64; % method 1
e0 = 0.5; %Determined experimentally for a flat plate with AR =4, Lambda = 1
% method 2 for finding Oswald's Efficiency Factor
Q = 1.05;
P = 0.007;

e0_Method2 = 1/(Q+P*pi*AR);

% arbitrary constant (eq. 5)
k1 = 1/(pi*e0*AR);

k1_Method2 = 1/(pi*e0_Method2*AR);

%Use different equation Cdo +kCl^2 (eq 4.)

CD = CD0 + k1.*(CL3D).^2;
CD_Method2 = CD0 + k1_Method2.*(CL3D).^2;

% whole aircraft drag (eq. 6)
%{
CD = CDmin + k1.*(CL3D - CLminD).^2;

CD_Method2 = CDmin + k1_Method2.*(CL3D - CLminD).^2;
%}

% plot for Method 1
figure(5)
hold on
plot(CL3D,CD,'bo','LineWidth',1); % whole aircraft drag polar
xline(0); % y axis
yline(0); % x axis
yline(CD0, '-.b');% parasite drag
% labels
title('C_L vs C_D');
ylabel('C_D');
xlabel('C_L');
legend('C_D Total', 'C_{D0}','Location','SouthEast');
hold off

% plot for Method 2
figure(6)
hold on
plot(CL3D,CD_Method2,'bo','LineWidth',1); % whole aircraft drag polar
xline(0); % y axis
yline(0); % x axis
yline(CD0, '-.b');% parasite drag
% labels
title('C_L vs C_D(Method 2)');
ylabel('C_D');
xlabel('C_L');
legend('C_D Total', 'C_{D0}','Location','SouthEast');
hold off

%plot for CD vs alpha
figure (7)
hold on
plot(AoA,CD,'bo','LineWidth',1);
xline(0);
yline(0);
yline(CD0,'-.b');
%labels
title('\alpha vs C_D');
ylabel('C_D');
xlabel('\alpha');
legend('C_D', 'C_{D0}','Location','SouthEast');
%% Comparisons

%Glide Range
% Solve Cl from drag
[ClCdMAX,iGR] = max(CL3D./CD);

h = 7; % [m]
wConstant = 0.295*9.81; %[N/m^2]
wIndv = 3.078; % [N] 

%Individual Weight Components [N]
wFuselage = wConstant * 0.364;
wBallast = 0.4;
wWings = wConstant * Sref;
wPayload = 1.568;
wTail = Stail / 2 * wConstant;

%Total Weight [N]
w = wFuselage + wBallast + wWings + wPayload + wTail;

WSIndv = w/Sref;
%density
rho = 1.16; % [kg/m^3]
%find center of gravity

%Stability

tailangle = 35;

xcgwing = 0.0575;
xcgfuse = 0.15;
xcgtailH = 0.46;
xcgtailV = NaN;
SH = Stail*cosd(tailangle);
SV = Stail*sind(tailangle);
xcgtail = (xcgtailH*SH) + (xcgtailV*SV)/(SH + SV);
xcgballast = -0.349;
xcgpayload = 0.0825;
xcgTest = (wWings*xcgwing + wFuselage*xcgfuse + wBallast*xcgballast ...
    + wPayload*xcgpayload)/(w);
xcg = (wWings*xcgwing + wFuselage*xcgfuse + wTail*xcgtailH + wBallast*xcgballast ...
    + wPayload*xcgpayload)/(w);



SH = Stail*cosd(tailangle);
changeXH = 0.4;
wingChord = c;
changeXV = 0.25*sind(tailangle)*0.5+0.06;

%Have to shift the center of gravity 10 cm down from the top of the
%fuselage. 

wingspan = sqrt(AR * Sref);

VH = SH*changeXH/(Sref*c);
SV = Stail*sind(tailangle);
VV = (SV*changeXV/(Sref*wingspan));

angle = 15;
CL_desired = 0.4;
B = (changeXV*angle)/(wingspan*CL_desired);

fprintf("VH is %f, VV is %f, B is %f", VH, VV,B)

pinf = 1.0581; % [kg/m^3] 1.5km STD ATM

% Max Glide Range Velocity = Vinf
Vinf = sqrt((2*w)/(CL3D(iGR)*Sref*pinf)); % [m/s]

% Max Glide Range (lecture 4)
%   steady, unaccelerated, power off glide
MaxRange = h*(ClCdMAX); % [m]

%Velocity at Max Glide Range
vMaxR = sqrt(2*w/(rho*sqrt(CD0/k1)*Sref));

% Glide Endurance
% glide angle small angle approx. (lecture 4)
theta = 1/(ClCdMAX); % [rad]

% solving for sink rate vinf*D/W = PR/W (lecture 4)
SinkRate = Vinf*sin(theta); % [m/s]

% Max Glide Endurance (lecture 4)
MaxEndurance = h/SinkRate; % [s]

%Velocity at Max Glide Endurance
vMaxE = sqrt(2*w/(rho*sqrt(CD0/(3*k1))*Sref));

%% Sizing

%create W/S vector

%create Sref vec
sRef = 0.1:0.05:1.2;

%find a test weight
wWingTest = wConstant .* sRef;
w1 = wWingTest + wFuselage + wPayload + wBallast +wTail;

WS = w1./sRef;

Rmax = h./(2.*(k1^(1/2)).*(Cfe.*(2+Sother./(sRef))).^(1/2));

Cd0 = Cfe.*(2+(Sother./sRef));

Emax = h.*WS./(2.*Cd0.*rho.*(sqrt(2.*WS./(sqrt(3.*Cd0./k1).*rho))).^3);

ClmaxR = sqrt(Cfe.*(2+Sother./(w1.*1./(WS)))./k1);
ClmaxE = sqrt(3.*Cd0./k1);

%plot Rmax, ClmaxR
figure

hold on

plot(WS,Rmax,'r');
%plot min range
yline(70,'r--');
ylabel('Max Range');
xlabel('Wing Loading');

yyaxis right
plot(WS,ClmaxR,'b');
%plot Cl stall with it
yline(max(CL3D),'b--');
ylabel('Cl');
legend('Max Range','Minimum Required Max R','Cl required for R', 'Cl stall');

title('Maximum Range Versus Wing Loading');
hold off

%Plot Emax, ClmaxE, etc.
figure

hold on

plot(WS,Emax,'r');
%plot min endurance
yline(12,'r--');
ylabel('Max Endurance');
xlabel('Wing Loading');


yyaxis right
plot(WS,ClmaxE,'b');
%plot Cl stall
yline(max(CL3D),'b--');
ylabel('Cl');
legend('Max Endurance','Min Required Max E','Cl required for E','Cl stall');

title('Maximum Endurance versus Wing Loading');
hold off

%Plot the velocities compared to Vstall
vRmax = sqrt(2./(rho.*sqrt(Cd0./k1))).*sqrt(WS);
vEmax = sqrt(2./(rho.*sqrt(3.*Cd0./k1))).*sqrt(WS);
rMaxV = 12;
rMinV = 7;

%plot
figure
hold on
plot(WS,vRmax,'r');
plot(WS,vEmax,'g');
yline(rMaxV,'b--');
yline(rMinV,'b-.');

ylabel('Velocity');
xlabel('Wing Loading');

legend('Velocity for Max Range', 'Velocity for Max Endurance', ...
    'Maximum Velocity', 'Minimum Velocity');

title('Maximum Velocity versus Wing Loading');
hold off

%Variations in Wing Geometry
b = sqrt(sRef*AR);
bMax = 1.0;

figure
hold on
plot(WS,b,'r');
yline(bMax,'--');
legend('Wingspan', 'Max Wingspan');
ylabel('Wingspan [m]');
xlabel('Wing Loading');
hold off

%% Stability




