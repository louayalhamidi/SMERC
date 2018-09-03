%% Inputs
clear all
% Global
dc = 53.848; % Combustion Chamber Diameter (mm)
dcdp = 6; % CC diameter/Pintle Diameter Ratio
Lsdp = 1; % Pintle Skip Length/Pintle Diameter Ratio
N_or = 20; % Total number of pintle orifices
P_or = 10; % Percentage of propellant going to secondary orifice rows

% Fuel
rho_f = 274.03; % Fuel Density (kg/m^3)
m_f = 0.24169; % Fuel Mass Flow Rate (kg/s)
Cd_f = 0.7; % Fuel Discharge Coefficient
dP_f = 130; % Fuel Pressure Drop (PSI)

% Oxidizer
rho_o = 1200; % Oxidizer Density
m_o = 0.410887; % Oxidizer Mass Flow Rate (kg/s)
Cd_o = 0.7; % Oxidizer Discharge Coefficient
dP_o = 130; % Oxidizer Pressure Drop (PSI)

% Film Cooling
FCP = 30; % Film Cooling Percentage of Fuel
Cd_fc = 0.7; % Film Cooling Discharge Coefficient
rho_fc = 274.03; % Film Cooling density (kg/m^3)
dP_fc = 130; % Film cooling Pressure Drop (PSI)
N_fc = 50; % Number of film cooling orifices

%% Flow Calculations

% Fuel
Q_f = m_f/rho_f; % Fuel Volumetric Flow rate (m^3/s)
A_total_f = m_f/(Cd_f*(sqrt(2*rho_f*dP_f*6894.76))); % Annular Total Injection Area (mm^2)
V_inj_f = Q_f/A_total_f; % Fuel Injection Velocity (m/s)


% Primary/Secondary Flow Calcualtions
m_o_s = m_o*(P_or/100);
m_o_p = m_o - m_o_s;

% Oxidizer (Primary Row)
Q_o_p = m_o_p/rho_o;
A_total_o_p = m_o_p/(Cd_o*(sqrt(2*rho_o*dP_o*6894.76))); % Center Total Area
V_inj_o_p = Q_o_p/A_total_o_p; 

% Oxidizer (Secondary Row)
Q_o_s = m_o_s/rho_o;
A_total_o_s = m_o_s/(Cd_o*(sqrt(2*rho_o*dP_o*6894.76))); % Center Total Area
V_inj_o_s = Q_o_s/A_total_o_s; 

% Film Cooling
m_fc = (FCP/100)*m_f; % Film Cooling Mass Flow Rate (kg/s)
Q_fc = m_fc/rho_fc; 
A_total_fc = m_fc/(Cd_fc*(sqrt(2*rho_fc*dP_fc*6894.76)));
V_inj_fc = Q_fc/A_total_fc;
%% Geometry Calculations

% Fuel
R_inner = ((dc/dcdp)/2)*0.001; % (m)
R_outer = sqrt((A_total_f/pi)+((R_inner)^2)); % Outer radius of annular channel (m)
Annular_Gap = R_outer - R_inner; % Annular Gap (m)

% Oxidizer
A_primary_c = A_total_o_p/(N_or/2); % (m^2)
D_primary_c = sqrt((4*A_primary_c)/pi); % (m)
A_secondary_c = A_total_o_s/(N_or/2); % (m^2)
D_secondary_c = sqrt((4*A_secondary_c)/pi); % (m)

% Film Cooling
A_fc = A_total_fc/N_fc; % (m^2)
D_fc = sqrt((4*A_fc)/pi); % (m)

%% Pintle Performance Calculations

dp = (dc/dcdp)*0.001; % Pintle Diameter (m)
Ls = (dp*Lsdp)*0.001; % Skip Length (m)
BF_p = (D_primary_c*(N_or/2))/(pi*dp); % Blockage Factor
BF_s = (D_secondary_c*(N_or/2))/(pi*dp); % Blockage Factor
TMR = ((m_o_p*V_inj_o_p)+(m_o_s*V_inj_o_s))/(m_f*V_inj_f); % Total Momentum Ratio (centered/annular)

%% Report

RESULTS = [];
RESULTS{1,1} = 'FUEL (ANNULAR)';
RESULTS{1,2} = 'METRIC (mm, m/s)';
RESULTS{1,3} = 'IMPERIAL (in, ft/s)';
RESULTS{2,1} = 'V_inj';
RESULTS(2,2) = num2cell(V_inj_f);
RESULTS(2,3) = num2cell(V_inj_f*3.28084);
RESULTS{3,1} = 'Annular Gap';
RESULTS(3,2) = num2cell(Annular_Gap*1000);
RESULTS(3,3) = num2cell(Annular_Gap*39.3701);
RESULTS{4,1} = 'OXIDIZER (RADIAL)';
RESULTS{5,1} = 'V_inj';
RESULTS(5,2) = num2cell(V_inj_o_p);
RESULTS(5,3) = num2cell(V_inj_o_p*3.28084);
RESULTS{6,1} = 'D_primary';
RESULTS(6,2) = num2cell(D_primary_c*1000);
RESULTS(6,3) = num2cell(D_primary_c*39.3701);
RESULTS{7,1} = 'D_secondary';
RESULTS(7,2) = num2cell(D_secondary_c*1000);
RESULTS(7,3) = num2cell(D_secondary_c*39.3701);
RESULTS{8,1} = 'FILM COOLING';
RESULTS{9,1} = 'V_inj';
RESULTS(9,2) = num2cell(V_inj_fc);
RESULTS(9,3) = num2cell(V_inj_fc*3.28084);
RESULTS{10,1} = 'D_fc';
RESULTS(10,2) = num2cell(D_fc*1000);
RESULTS(10,3) = num2cell(D_fc*39.3701);
RESULTS{11,1} = 'PERFORMANCE/PINTLE';
RESULTS{12,1} = 'd_p';
RESULTS(12,2) = num2cell(dp*1000);
RESULTS(12,3) = num2cell(dp*39.3701);
RESULTS{13,1} = 'L_s';
RESULTS(13,2) = num2cell(Ls*1000);
RESULTS(13,3) = num2cell(Ls*39.3701);
RESULTS{14,1} = 'BF_p';
RESULTS(14,2) = num2cell(BF_p);
RESULTS{15,1} = 'BF_s';
RESULTS(15,2) = num2cell(BF_s);
RESULTS{16,1} = 'TMR';
RESULTS(16,2) = num2cell(TMR);

filename = 'PINTLE_RESULTS_KRAKEN_350';
xlswrite(filename,RESULTS);
