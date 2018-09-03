%% Film Cooling Inputs
clear all
close all

N_stations = 25;

m_dot_f = .53285;
m_dot_o = .90585;

FCP = 30;          % film cooling percentage (%)
K_t = .003;         % Turbulent mixing intensity
H_s = .0775; %surface layer thickness (in)

T_f_inj=270; % (R)
T_vap= 341.52; % Vaporization temp at 650 psi (R)
Q_vap = 105.355/2; %0.07076526225279 Heat/Enthalpy of Vaporization (BTU/lb) = 2.6404 kJ/kmol
e = exp(1);
%% Creating Contour and Other Functions

%Contour
% contour_table = readtable('C:\Users\louay\Desktop\SMERC\contour_1005.txt');
contour_table = readtable('C:\Users\louay\Desktop\SMERC\contour_350.txt');

contour = table2array(contour_table);
x = contour(:,1);       % creating stream-wise distance array
r = contour(:,2);       % creating chamber radius array

%Pressure Function
% pressure_table = readtable('C:\Users\louay\Desktop\SMERC\Pressurefunction.txt');
pressure_table = readtable('C:\Users\louay\Desktop\SMERC\Pressurefunction_350.txt');

pressure = table2array(pressure_table);
[x_p I] = unique(pressure(:,1));
p = pressure(:,2);
p = p(I);

%Heat Flux at Vaporization temperature function
q_conv_T_vap_table = readtable('C:\Users\louay\Desktop\SMERC\q_T_vap_350.txt');
q_conv_T_vap_table_num = table2array(q_conv_T_vap_table);
[x_q Iq] = unique(q_conv_T_vap_table_num(:,1));
q_vap = q_conv_T_vap_table_num(:,4);
q_vap = q_vap(Iq);

%% Initialize

%Mixing Criteria
m_dot_fc = ((FCP/100)*m_dot_f);  % film cooling mass flow rate (lb/s)
m_dot_fi_rel = zeros(1,N_stations);
m_dot_s_rel = zeros(1,N_stations);
m_dot_tot(1,1:N_stations) = m_dot_fc+m_dot_f+m_dot_o;


%Heat fluxes
q_conv_T_f = []; %Convective heat flux to liquid coolant film at coolant film temp (BTU/in^2.s)
q_conv_T_vap = []; %Convective heat flux to liquid coolant film at coolant vaporization temp (BTU/in^2.s)
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

%Get relative x vector
for i = 1:N_stations
    if i == 1
        x_length_rel(i) = x_length(i);
    end
    if i > 1
        x_length_rel(i) = x_length(i)+x_length_rel(i-1);
    end
end
%% Film Cooling Calculations

for i = 1:N_stations
    
    if i == 1
        T_f(i) = T_f_inj;
    end
    
%% Heating Section
    if T_f(end) <= T_vap
        
        P_f(i) = interp1(x_p,p,x_length_rel(i));
        
%         url_string1 = sprintf('http://webbook.nist.gov/cgi/fluid.cgi?Action=Load&ID=C74828&Type=IsoBar&Digits=5&P=%f&THigh=%f&TLow=%f',P_f(i),T_f(i),T_f(i));
%         url_string2 = ('&TInc=.1&RefState=DEF&TUnit=R&PUnit=psia&DUnit=lbm%2Fft3&HUnit=Btu%2Flbm&WUnit=ft%2Fs&VisUnit=lbm%2Fft*s&STUnit=lb%2Fin');
%         url = strcat(url_string1,url_string2);
        
        url_string1 = sprintf('https://webbook.nist.gov/cgi/fluid.cgi?Action=Data&Wide=on&ID=C74828&Type=IsoBar&Digits=5&P=%f&THigh=%f&TLow=%f',P_f(i),T_f(i),T_f(i));
        url_string2 = ('&TInc=.1&RefState=DEF&TUnit=R&PUnit=psia&DUnit=lbm%2Fft3&HUnit=Btu%2Flbm&WUnit=ft%2Fs&VisUnit=lbm%2Fft*s&STUnit=lb%2Fin');
        url = strcat(url_string1,url_string2);

        
        filename = 'adastra.txt';
        options.Timeout = 60;
        outfilename = websave(filename,url);
        T = readtable(filename);
        delete 'adastra.txt';
        NIST_data = table2array(T(1,1:13));
    
        Cp_f(i)=NIST_data(1,9); %BTU/lb.R
        mu_f(i)=NIST_data(1,12)*(1/12); %1 lb/ft-s = 1/12 lb/in-s
        q_conv_T_f(i) = input([' Enter the convective heat transfer for film temperature of ' num2str(T_f(i)) ':']);
        Re_f(i) = (m_dot_fc)/(2*pi*(D_avg(i)/2)*mu_f(i));
        nu(i) = input([' Enter film effectivness for ' num2str(Re_f(i)) ':']);
        dT_f(i) = (2*pi*(D_avg(i)/2)*q_conv_T_f(i)*L(i))/(m_dot_fc*Cp_f(i)*nu(i));
        T_f(i+1) = T_f(i) + dT_f(i);
        m_dot_f_tot = zeros(1,i);
        
        film_status(i) = 0;

    end
    
%% Vaporization Section
    if (T_f(end) > T_vap) & (m_dot_f_tot(end) < m_dot_fc)
        
        
        q_conv_T_vap(i) = interp1(x_q,q_vap,x_length_rel(i));
        dm_dot_f(i) = (2*pi*(D_avg(i)/2)*q_conv_T_vap(i)*L(i))/(Q_vap);
        m_dot_f_tot(i+1) = m_dot_f_tot(i)+dm_dot_f(i);
        film_status(i) = 1;

        x_mix = zeros(1,i);
        zeta = zeros(1,i);

    end
    
%% Mixing Section
    if (m_dot_f_tot(end) >= m_dot_fc) & (zeta(end) <=1)
        

            
        if x_mix(end) == 0
            m_dot_fi(i) = m_dot_fc;
            m_dot_s(i) = m_dot_f+m_dot_o;
            m_dot_fi_rel(i) = m_dot_fi(i)/m_dot_tot(i);
            m_dot_s_rel(i) = m_dot_s(i)/m_dot_tot(i);
            
            x_mix = zeros(1,i);
            zeta = zeros(1,i);
            
            x_mix(i+1) = L(i);
            
            r_s(1,1:i) = m_dot_o/m_dot_f;
            r_fi(1,1:i) = 0;

            continue
    
        end
        
        if x_mix(end) > 0
            
            m_dot_s_rel(i) = (m_dot_s_rel(i-1)*(1-(zeta(i-1)/2)))+(m_dot_fi_rel(i-1)*(zeta(i-1)/2));
            m_dot_fi_rel(i) = (m_dot_fi_rel(i-1)*(1-(zeta(i-1)/2)))+(m_dot_s_rel(i-1)*(zeta(i-1)/2));
            
            k_ss(i) = (m_dot_s_rel(i-1)/m_dot_s_rel(i))*(1-(zeta(i-1)/2));
            k_sf(i) = (m_dot_fi_rel(i-1)/m_dot_s_rel(i))*(zeta(i-1)/2);
            k_ff(i) = (m_dot_fi_rel(i-1)/m_dot_fi_rel(i))*(1-(zeta(i-1)/2));
            k_fs(i) = (m_dot_s_rel(i-1)/m_dot_fi_rel(i))*(zeta(i-1)/2);
            
            r_s(i) = (k_ss(i)*r_s(i-1))+(k_sf(i)*r_fi(i-1));
            r_fi(i) = (k_fs(i)*r_s(i-1))+(k_ff(i)*r_fi(i-1));
            
            M(i) = K_t*(m_dot_s_rel(i)/m_dot_fi_rel(i));
            x_mix(i) = x_mix(i-1)+L(i-1);
            x_bar(i) = (x_mix(i)/H_s);
            zeta(i) = 1-(e^(-M(i)*x_bar(i)));
            
        end
        
    film_status(i) = 2;    
    
    end
    
    if (m_dot_f_tot(end) >= m_dot_fc) & (zeta(i) > 1)
        film_status(i) = 3 ;
    end
    
    
end

%% Output

film_RESULTS = [];
film_RESULTS{1,1} = 'Station';
film_RESULTS{2,1} = 'x - location (in)';
film_RESULTS{3,1} = 'Average Diameter (in)';
film_RESULTS{4,1} = 'Chamber Pressure (psi)';
film_RESULTS{5,1} = 'Liquid Film Temp. (R)';
film_RESULTS{6,1} = 'Convective Heat Flux to Liquid Film (BTU/in^s.s)';
film_RESULTS{7,1} = 'Liquid film Reynolds Number';
film_RESULTS{8,1} = 'Liquid film viscosity (lb/in-s)';
film_RESULTS{9,1} = 'Liquid film Specific Heat (BTU/lb.R)';
film_RESULTS{10,1} = 'Vaporized Film (lb/s)';
film_RESULTS{11,1} = 'Convective Heat Flux to mix film (BTU/in^2.s)';
film_RESULTS{12,1} = 'Film Status';
film_RESULTS{13,1} = 'Film Mixture Ratio';

film_RESULTS(1,2:N_stations+1) = num2cell(1:N_stations);
film_RESULTS(2,2:N_stations+1) = num2cell(x_length_rel);
film_RESULTS(3,2:N_stations+1) = num2cell(D_avg);
film_RESULTS(4,2:size(P_f,2)+1) = num2cell(P_f);
film_RESULTS(5,2:size(T_f,2)+1) = num2cell(T_f);
film_RESULTS(6,2:size(q_conv_T_f,2)+1) = num2cell(q_conv_T_f);
film_RESULTS(7,2:size(Re_f,2)+1) = num2cell(Re_f);
film_RESULTS(8,2:size(mu_f,2)+1) = num2cell(mu_f);
film_RESULTS(9,2:size(Cp_f,2)+1) = num2cell(Cp_f);
film_RESULTS(10,2:size(m_dot_f_tot,2)+1) = num2cell(m_dot_f_tot);
film_RESULTS(11,2:size(q_conv_T_vap,2)+1) = num2cell(q_conv_T_vap);
film_RESULTS(12,2:N_stations+1) = num2cell(film_status);
film_RESULTS(13,2:N_stations+1) = num2cell(r_fi);

filename = sprintf('film_RESULTS_305');
xlswrite(filename,film_RESULTS);


save('C:\Users\louay\Desktop\SMERC\r_fi_350.mat','r_fi');
save('C:\Users\louay\Desktop\SMERC\film_status_350.mat','film_status');
