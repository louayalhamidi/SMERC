%% Global Inputs
clear all;
close all;

N_stations = 25; % Number of Stations

FCP = 30;

load('C:\Users\louay\Desktop\SMERC\r_fi_350.mat'); %Load Film Mixture Ratios
load('C:\Users\louay\Desktop\SMERC\film_status_350.mat'); %Load Film Status

%Cooling Inputs
m_dot_f = .53285;      % fuel mass flow rate (lb/s)
k_avg = 23;        % initial wall thermal conductivity (W/m.K)
P_co_in(N_stations) = 900; % coolant pressure entering station (psi)
T_co_in(N_stations) = 180; % coolant temperature entering station (R)
T_cc = 343.8; % methane critical temperature (R)
g = 386.088583;           % gravitational acceleration (in/s)
f = .02;            % friction coefficient
t = .0245;            % wall thickness (in)

%Structural Inputs
s_y_min = 46000;      % minimum material yield strength (psi) (1/2 of yield)
E = 20000000 ;            % Youngs modulus (psi)
nu = .284;              % Poissons Ratio
P_c(1,1:N_stations) = 670;

N_table = readtable('C:\Users\louay\Desktop\SMERC\N_vector_350.txt');
N_array = table2array(N_table);
N_array = N_array';
%% Creating Contour

contour_table = readtable('C:\Users\louay\Desktop\SMERC\contour_350.txt');
contour = table2array(contour_table);
x = contour(:,1);       % creating stream-wise distance array
r = contour(:,2);       % creating chamber radius array

%% Mass Flow Rate Calculations
m_dot_t = (1+(FCP/100))*m_dot_f;        % total mass flow rate w/film cooling (lb/s)
m_dot_fc = (FCP/100)*m_dot_f;           % film cooling mass flow rate (lb/s)
%% Initializing Variables

%Geometric Variables.
D_avg = []; % station diameter (in)
A = [];     % station surface area (in^2)
L = [];         % station length (in)

%Coolant Properties Variables
k_co = [];      % coolant thermal conductivity (BTU/in^2-s-F/in) 1 W/m-K = .00001337474527 BTU/in^2-s-F/in
mu_co = [];     % coolant dynamic viscosity (lb/in.s) 1 lb/ft-s = 1/12 lb/in-s
Cp_co = [];     % coolant specific heat (BTU/lb-R)  
rho_co = [];       % coolant density (lb/in^3) 1 lb/ft^3 = 1/1728 lb/in^3

%Coolant Temperature Variabes
dT_co = [];     % coolant temperature gain (R)
T_co_out = [];  % coolant temperature exiting station (R)
T_co_avg  = [];  % average coolant temperature in station (R)

%Coolant Flow Variables
V_co = [];      % coolant velocity (in/s)
dP = [];        % pressure drop (psi)
P_co_out  = [];  % coolant pressure exiting station (psi)
P_co_avg  = [];  % average coolant pressure in station (psi)

%Wall Temperature Variables
T_wc = [];  % cooling side wall temperature (R)
T_wa = [];  % average wall temperature (R)
dT_w = []; % temperature differential across the wall (R)
k_w = [];   % actual wall thermal conductivity ()1 W/m-K = .00001337474527 BTU/in^2-s-F/in

%Heat Transfer Variables
Q = [];         % heat absorbed by coolant (BTU/s)
h_c = [];       % coolant convective heat transfer coefficient (BTU/in^2-s-F)

%Structural Variables
P_crit = [];     % critical buckling pressure (psi)
s_h_co = [];    % coolant hoop stress (psi)

%Output Variables
C = [];
d = [];


%% Discretization

% Discretize the geometry based on the number of stations

for i = N_stations:-1:1
    gd = [];
    z=1-(200/N_stations);
    for k = 1:N_stations+1
        gd(k) = z+(200/N_stations);
        z = z + (200/N_stations);
    end

    D_avg(i) = r(gd(i))+r(gd(i+1));
    A(i) = pi*(r(gd(i))+r(gd(i+1)))*sqrt(((r(gd(i))-r(gd(i+1)))^2)+(x(gd(i))-x(gd(i+1)))^2);
    x_length(i) = abs(x(gd(i))-x(gd(i+1)));
    L(i) = sqrt(((r(gd(i))-r(gd(i+1)))^2)+(x(gd(i))-x(gd(i+1)))^2);
    
end

for i = 1:N_stations
    if i == 1
        x_length_rel(i) = x_length(i);
    end
    if i > 1
        x_length_rel(i) = x_length(i)+x_length_rel(i-1);
    end
end

for i = 1:N_stations
    if i == 1
        x_length_act(i) = x_length(i)/2;
    end
    if i > 1
        x_length_act(i) = sum(x_length(1,1:i-1)) + (x_length(i)/2);
    end
end

%% Section Calculations

% Discretize the geometry based on the number of stations
gd = [];
z=1-(200/N_stations);
for k = 1:N_stations+1
    gd(k) = z+(200/N_stations);
    z = z + (200/N_stations);
end

% Begin Cooling Calculations

for i = N_stations:-1:1
    
    % Find T_wg, q_conv, q_rad
    name_add = sprintf('s_%i',i);
    table_name = strcat('C:\Users\louay\Desktop\SMERC\tables\',name_add);
    heat_table = readtable(table_name);
    heat = table2array(heat_table);
    
    x_rpa = heat(:,1);
    q_conv_rpa = heat(:,2);
    q_rad_rpa = heat(:,3);
    T_wg_rpa = heat(:,4);
    
    [x_rpa, index] = unique(x_rpa);
    
    q_conv(i) = interp1(x_rpa,q_conv_rpa(index),x_length_act(i),'spline');
    q_rad(i) = interp1(x_rpa,q_rad_rpa(index),x_length_act(i),'spline');
    T_wg(i) = interp1(x_rpa,T_wg_rpa(index),x_length_act(i),'spline');
    
    % Apply Film Status
    
    if film_status(i) == 0 || film_status(i) ==1
        q(i) = q_rad(i);
    end
    
    if film_status(i) == 2
        q(i) = q_rad(i) + q_conv(i);
    end
    
    if i < N_stations
        P_co_in(i) = P_co_out(i+1);
        T_co_in(i) = T_co_out(i+1);
    end

    % Wall Transfer Calculations
    
    T_wc(i) = T_wg(i)-(q(i)*(t/(k_avg*.00001337474527)));
    T_wa(i) = (T_wc(i)+T_wg(i))/2;
    k_w = .015*T_wa + 11.002;
    T_wc(i) = T_wg(i)-(q(i)*(t/(k_w(i)*.00001337474527)));
    T_wa(i) = (T_wc(i)+T_wg(i))/2;
    dT_w = T_wg - T_wc;
    
    % Get Data from NIST Thermophysical Workbook website
    
    url_string1 = sprintf('https://webbook.nist.gov/cgi/fluid.cgi?Action=Data&Wide=on&ID=C74828&Type=IsoBar&Digits=5&P=%f&THigh=%f&TLow=%f',P_co_in(i),T_co_in(i),T_co_in(i));
    url_string2 = ('&TInc=.1&RefState=DEF&TUnit=R&PUnit=psia&DUnit=lbm%2Fft3&HUnit=Btu%2Flbm&WUnit=ft%2Fs&VisUnit=lbm%2Fft*s&STUnit=lb%2Fin');
    url = strcat(url_string1,url_string2);

    filename = 'adastra.txt';
    options.Timeout = 6000;
    outfilename = websave(filename,url);
    T = readtable(filename);
    delete 'adastra.txt';
    NIST_data = table2array(T(1,1:13));
    
    Cp_co(i)=NIST_data(1,9); %BTU/lb.R
    k_co(i)=NIST_data(1,13)*.00001337474527; % 1 W/m-K = .00001337474527 BTU/in^2-s-F/in
    mu_co(i)=NIST_data(1,12)*(1/12); %1 lb/ft-s = 1/12 lb/in-s
    rho_co(i)=NIST_data(1,3)*(1/1728); %1 lb/ft^3 = 1/1728 lb/in^3
    
    % Preform Heat Transfer Calculations
 
    Q(i) = q(i)*A(i);
    dT_co(i) = (A(i)*q(i))/(m_dot_t*Cp_co(i));
    T_co_out(i) = T_co_in(i) + dT_co(i);
    T_co_avg(i) = (T_co_in(i) + T_co_out(i))/2;
    
    h_c(i) = q(i)/(T_wc(i)-T_co_avg(i));
    
    C(i) = (.0185*(k_co(i)/h_c(i))*(((4*m_dot_t)/(pi*mu_co(i)))^.8)*(((mu_co(i)*Cp_co(i))/k_co(i))^.4)*((T_co_avg(i)/T_wc(i))^.1))^1.25;
    d(i) = (N_array(i)/C(i))^(1/(-2.25));
    
    % Fluid Analysis
    
    V_co(i) = (4*m_dot_t)/(pi*N_array(i)*(d(i)^2)*rho_co(i));
    dP_co(i) = f*(L(i)/d(i))*((rho_co(i)*(V_co(i)^2))/(2*g));
    P_co_out(i) = P_co_in(i) - dP_co(i);
    P_co_avg(i) = (P_co_in(i)+P_co_out(i))/2;
    
    %Structural Analysis
    
    s_h_co(i) = (P_co_avg(i)*(d(i)/2))/t;
    s_h_c(i) = ((P_c(i)*(D_avg(i)/2))/.12);
    display(i);
end
    
%% Structural Analysis

L_channel = sum(L);

for i = 1:N_stations
    P_crit_co(i) = (.855/(1-(nu^2))^(3/4))*(E/(((((D_avg(i)/2)/.1)^(5/2))*(L_channel/(D_avg(i)/2)))));
end

s_max_h = (2*max(s_h_co))+max(s_h_c);

RFoS = s_y_min/s_max_h;
strn = ['Realized FoS for Combined Hoop Stress is ' ,num2str(RFoS)];
disp(strn);

[P_crit_co_min ,I] = min(P_crit_co);

if P_co_avg(I) >= P_crit_co_min
    disp('Coolant Pressure Exceeds Critical Pressure - ENGINE WILL BUCKLE!!!\n');
    str = ['min Critical Pressure is ', num2str(P_crit_co_min),' at section ', num2str(I)];
    disp(str);
    str2 = ['Coolant Pressure is ' , num2str(P_co_avg(I))];
    disp(str2);
end

if s_max_h >= s_y_min
    disp('Minimum Yield Stress is Exceeded!!!');
    str3 = ['Maximum Hoop Stress is ',num2str(s_max_h)]; 
end


%% Global Calculations
                                                                                                                                                                                                                                                                
x_length_tot = sum(x_length);
A_engine = sum(A);
Cp_co_avg = sum(Cp_co)/length(Cp_co);
Q_c = m_dot_t*Cp_co_avg*(T_cc - T_co_in(N_stations));
Q_act = sum(Q);
Q_p = (Q_act/Q_c)*100;
dT_co_tot = sum(dT_co);
dP_co_tot = sum(dP_co);

%% Output

%Station
RESULTS = [];
RESULTS{1,1} = 'Station';
RESULTS{2,1} = 'x-length (in)';
RESULTS{3,1} = 'Average Diameter (in)';
RESULTS{4,1} = 'Station Area (in^2)';
RESULTS{5,1} = 'Channel Length (in)';
RESULTS{6,1} = 'Total Heat Flux (BTU/in^2.s)';
RESULTS{7,1} = 'Convective Heat Flux (BTU/in^2.s)';
RESULTS{8,1} = 'Radiation Heat Flux (BTU/in^2.s)';
RESULTS{9,1} = 'Gas-Side Wall Temp (R)';
RESULTS{10,1} = 'Cooling-Side Wall Temp (R)';
RESULTS{11,1} = 'Average Wall Temp (R)';
RESULTS{12,1} = 'Temp Diff Across Wall (R)';
RESULTS{13,1} = 'Wall Thermal Conducivity';
RESULTS{14,1} = 'Entering Coolant Temp (R)';
RESULTS{15,1} = 'Coolant Temp Gain (R)';
RESULTS{16,1} = 'Exiting Coolant Temp (R)';
RESULTS{17,1} = 'Average Coolant Temp (R)';
RESULTS{18,1} = 'Entering Coolant Pressure (psi)';
RESULTS{19,1} = 'Coolant Pressure Loss (psi)';
RESULTS{20,1} = 'Exiting Coolant Pressure (psi)';
RESULTS{21,1} = 'Average Coolant Pressure (psi)';
RESULTS{22,1} = 'Coolant Velocity (in/s)';
RESULTS{23,1} = 'Coolant Density (lb/in^3)';
RESULTS{24,1} = 'Coolant Thermal Conductivity (BTU.in/in^2.s.F)';
RESULTS{25,1} = 'Coolant Specific Heat (BTU/lb.R)';
RESULTS{26,1} = 'Coolant Viscosity (lb/in-s)';
RESULTS{27,1} = 'Coolant Heat Gain (BTU/s)';
RESULTS{28,1} = 'Coolant Convective Heat Transfer Coefficient (BTU/in^2.s.F)';
RESULTS{29,1} = 'Heat Transfer Constant';
RESULTS{30,1} = 'Channel Diameter (in)';

RESULTS(1,2:N_stations+1) = num2cell(1:N_stations);
RESULTS(2,2:N_stations+1) = num2cell(x_length);
RESULTS(3,2:N_stations+1) = num2cell(D_avg);
RESULTS(4,2:N_stations+1) = num2cell(A);
RESULTS(5,2:N_stations+1) = num2cell(L);
RESULTS(6,2:N_stations+1) = num2cell(q);
RESULTS(7,2:N_stations+1) = num2cell(q_conv);
RESULTS(8,2:N_stations+1) = num2cell(q_rad);
RESULTS(9,2:N_stations+1) = num2cell(T_wg);
RESULTS(10,2:N_stations+1) = num2cell(T_wc);
RESULTS(11,2:N_stations+1) = num2cell(T_wa);
RESULTS(12,2:N_stations+1) = num2cell(dT_w);
RESULTS(13,2:N_stations+1) = num2cell(k_w);
RESULTS(14,2:N_stations+1) = num2cell(T_co_in);
RESULTS(15,2:N_stations+1) = num2cell(dT_co);
RESULTS(16,2:N_stations+1) = num2cell(T_co_out);
RESULTS(17,2:N_stations+1) = num2cell(T_co_avg);
RESULTS(18,2:N_stations+1) = num2cell(P_co_in);
RESULTS(19,2:N_stations+1) = num2cell(dP_co);
RESULTS(20,2:N_stations+1) = num2cell(P_co_out);
RESULTS(21,2:N_stations+1) = num2cell(P_co_avg);
RESULTS(22,2:N_stations+1) = num2cell(V_co);
RESULTS(23,2:N_stations+1) = num2cell(rho_co);
RESULTS(24,2:N_stations+1) = num2cell(k_co);
RESULTS(25,2:N_stations+1) = num2cell(Cp_co);
RESULTS(26,2:N_stations+1) = num2cell(mu_co);
RESULTS(27,2:N_stations+1) = num2cell(Q);
RESULTS(28,2:N_stations+1) = num2cell(h_c);
RESULTS(29,2:N_stations+1) = num2cell(C);
RESULTS(30,2:N_stations+1) = num2cell(d);

%Global

RESULTS{31,1} = 'GLOBAL RESULTS';
RESULTS{32,1} = 'Total x length (in)';
RESULTS{33,1} = 'Total Channel Length (in)';
RESULTS{34,1} = 'Total Engine Area (in^2)';
RESULTS{35,1} = 'Total Average Specific Heat (BTU/lb.R)';
RESULTS{36,1} = 'Total Coolant Temperature Gain (R)';
RESULTS{37,1} = 'Total Coolant Pressure Drop (psi)';
RESULTS{38,1} = 'Cooling Capacity (BTU/s)';
RESULTS{39,1} = 'Total Heat Absorbed (BTU/s)';
RESULTS{40,1} = 'Critical Heat Gain Percentage (%)';
RESULTS{41,1} = 'Number of Cooling Channels';

RESULTS(32,2) = num2cell(x_length_tot);
RESULTS(33,2) = num2cell(L_channel);
RESULTS(34,2) = num2cell(A_engine);
RESULTS(35,2) = num2cell(Cp_co_avg);
RESULTS(36,2) = num2cell(dT_co_tot);
RESULTS(37,2) = num2cell(dP_co_tot);
RESULTS(38,2) = num2cell(Q_c);
RESULTS(39,2) = num2cell(Q_act);
RESULTS(40,2) = num2cell(Q_p);
RESULTS(41,2:N_stations+1) = num2cell(N_array);

filename = sprintf('TEST_RESULTS_P=%f_T=%f_FCP=%f_N=%f-%f_t=%f.xlsx',P_co_in(N_stations),T_co_in(N_stations),FCP,N_array(N_stations),N_array(((N_stations/2)+.5)+1),t);
xlswrite(filename,RESULTS);

%% Plots

%Heat Flux
left_color = [0 0 0];
right_color = [1 0 0];
set(figure,'defaultAxesColorOrder',[left_color; right_color]);
yyaxis left
plot(x,r,'k','linewidth',3);
xlabel('x - distance [in]');
ylabel('Radius [in]');
yyaxis right
plot(x_length_rel,q,'-.r','linewidth',2);
ylabel('q_{total} [BTU/in^{2}-s]');
set(gcf,'color','w');

%Gas side Wall temp
left_color = [0 0 0];
right_color = [1 0 0];
set(figure,'defaultAxesColorOrder',[left_color; right_color]);
yyaxis left
plot(x,r,'k','linewidth',3);
xlabel('x - distance [in]');
ylabel('Radius [in]');
yyaxis right
plot(x_length_rel,T_wg,'-.r','linewidth',2);
ylabel('T_{wg} [R]');
set(gcf,'color','w');

%Convective Heat Transfer Coefficient (Cooling)
left_color = [0 0 0];
right_color = [1 0 0];
set(figure,'defaultAxesColorOrder',[left_color; right_color]);
yyaxis left
plot(x,r,'k','linewidth',3);
xlabel('x - distance [in]');
ylabel('Radius [in]');
yyaxis right
plot(x_length_rel,h_c,'-.r','linewidth',2);
ylabel('h_{c} [BTU/in^{2}-s-F]');
set(gcf,'color','w');

%cooling channel diameter
left_color = [0 0 0];
right_color = [1 0 1];
set(figure,'defaultAxesColorOrder',[left_color; right_color]);
yyaxis left
plot(x,r,'k','linewidth',3);
xlabel('x - distance [in]');
ylabel('Radius [in]');
yyaxis right
plot(x_length_rel,d,'-.m','linewidth',2);
ylabel('Cooling Channel Diameter [in]');
set(gcf,'color','w');

%Coolant pressure
left_color = [0 0 0];
right_color = [0 0 1];
set(figure,'defaultAxesColorOrder',[left_color; right_color]);
yyaxis left
plot(x,r,'k','linewidth',3);
xlabel('x - distance [in]');
ylabel('Radius [in]');
yyaxis right
plot(x_length_rel,P_co_avg,'-.b','linewidth',2);
ylabel('P_{co} [PSI]');
set(gcf,'color','w');

%coolant Temperature
left_color = [0 0 0];
right_color = [0 0 1];
set(figure,'defaultAxesColorOrder',[left_color; right_color]);
yyaxis left
plot(x,r,'k','linewidth',3);
xlabel('x - distance [in]');
ylabel('Radius [in]');
ax.XAxis.Color = 'k';
yyaxis right
plot(x_length_rel,T_co_avg,'-.b','linewidth',2);
ylabel('T_{co} [R]');
set(gcf,'color','w');

%coolant velocity
left_color = [0 0 0];
right_color = [0 0 1];
set(figure,'defaultAxesColorOrder',[left_color; right_color]);
yyaxis left
plot(x,r,'k','linewidth',3);
xlabel('x - distance [in]');
ylabel('Radius [in]');
ax.XAxis.Color = 'k';
yyaxis right
plot(x_length_rel,V_co./12,'-.b','linewidth',2);
ylabel('v_{co} [ft/s]');
set(gcf,'color','w');

%coolant density
left_color = [0 0 0];
right_color = [0 0 1];
set(figure,'defaultAxesColorOrder',[left_color; right_color]);
yyaxis left
plot(x,r,'k','linewidth',3);
xlabel('x - distance [in]');
ylabel('Radius [in]');
ax.XAxis.Color = 'k';
yyaxis right
plot(x_length_rel,rho_co.*1728,'-.b','linewidth',2);
ylabel('\rho_{co} [lb/ft^3]');
set(gcf,'color','w');

%cooling hoop stress
left_color = [0 0 0];
right_color = [0 0 0];
set(figure,'defaultAxesColorOrder',[left_color; right_color]);
yyaxis left
plot(x,r,'k','linewidth',3);
xlabel('x - distance [in]');
ylabel('Radius [in]');
ax.XAxis.Color = 'k';
yyaxis right
plot(x_length_rel,s_h_co,'-.k','linewidth',2);
ylabel('s_{h,co} [PSI]');
set(gcf,'color','w');

%chamber hoop stress
left_color = [0 0 0];
right_color = [0 0 0];
set(figure,'defaultAxesColorOrder',[left_color; right_color]);
yyaxis left
plot(x,r,'k','linewidth',3);
xlabel('x - distance [in]');
ylabel('Radius [in]');
ax.XAxis.Color = 'k';
yyaxis right
plot(x_length_rel,s_h_c,'-.k','linewidth',2);
ylabel('s_{h,c} [PSI]');
set(gcf,'color','w');

close all;