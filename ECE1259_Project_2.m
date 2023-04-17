% This program is to simulate a coax cable design

format long e
clear all;
clc

%% opening statements
disp('This program is for a coxial cable design')
disp('------------------------------------------')
disp('ECE1249 Design Project 2')
disp('Authors: Rick, Seth, Andrew, Kamden')
fprintf('\n')

%% constants
h = 6.626e-34; % planks constant in Js
u_0 = pi*4e-7; % permeability constant in H/m
e_0 = 8.854e-12; % permitivitty constant F/m
c = 3e8; % speed of light in m/s
sweep = 1:1000000;

% grab these from the material data
sigma_d = 0;
sigma_c = 0;
e_r = 0;
u_r = 1;
Ebr = 0;
tan_d = 0;

%% inputs
%Listed conductors with sigma_c
Al = 3.8e7;
Car = 3e4;
Cu = 5.8e7;
Au = 4.1e7;
Gr = 7e4;
Fe = 1e7;
Pb = 5e6;
Nich = 1e6;
Ni = 1.5e7;
Ag = 6.2e7;
Sol = 7e6;
Stain = 1.1e6;
Sn = 8.8e6;
W = 1.8e7;

%Listed dielectrics with e_r, Ebr, tand, and sigma_d
Air = [1.0005, 3e6, 0, 0];
Alum = [9.9, 0, 0.0001, 0];
Bar = [1200, 7.5e6, 0, 0];
Glass = [10, 30e6, 0.004, 10^-12];
Ice = [4.2, 0, 0.12, 0];
Mica = [5.4, 200e6, 0.0003, 10^-15];
Poleth = [2.26, 47e6, 0, 10^-16];
Polsty = [2.56, 20e6, 0, 10^-17];
Quartz = [3.8, 30e6, 0.0002, 10^-17];
Si = [11.8, 0, 0, 4.4e-4];
Soil = [3, 0, 0.017, 2e-3]; %e_r between 3-4 so I put 3
Teflon = [2.1, 60e6, 0.0001, 10^-15]; %tand < 0.0002 so I put 0.0001
Water = [81, 0, 0.04, 10^-4];
Seawater = [72, 0, 0.9, 5];

%defining checking variable
correct = 'N';

%creating a menu for the user to give conductor 
conductor = menu('What conductor are you using?', 'Aluminum', 'Carbon', 'Copper', ...
    'Gold', 'Graphite', 'Iron', 'Lead', 'Nichrome', 'Nickel', 'Silver', 'Solder', ...
    'Stainless Steel', 'Tin', 'Tungsten', 'Other');

switch conductor
    case 1
        sigma_c = Al;
        
    case 2
        sigma_c = Car;
        
    case 3
        sigma_c = Cu;
       
    case 4
        sigma_c = Au;
        
    case 5
        sigma_c = Gr;
        
    case 6
        sigma_c = Fe;
       
    case 7
        sigma_c = Pb;
        
    case 8
        sigma_c = Nich;
        
    case 9 
        sigma_c = Ni;
        
    case 10
        sigma_c = Ag;
        
    case 11
        sigma_c = Sol;
        
    case 12
        sigma_c = Stain;
        
    case 13
        sigma_c = Sn;
        
    case 14 
        sigma_c = W;
        
    case 15
        sigma_c = input('Please enter a value for the conductor permiability in S/m: ');
        disp(['You chose Other with a permiability of ', num2str(sigma_c)]);
        
end

%Menu for the dielectric from list or given values.
dielectric = menu('What dielectric are you using?', 'Air', 'Alumina', 'Barium Titanate', ...
    'Glass', 'Ice', 'Mica', 'Polyethylene', 'Polystyrene', 'Quartz (fused)', 'Silicon (pure)', ...
    'Soil (dry)', 'Teflon', 'Water (distilled)', 'Seawater', 'Other');

switch dielectric
    case 1
        e_r = Air(1);
        Ebr = Air(2);
        tan_d = Air(3);
        sigma_d = Air(4);
        
    case 2
        e_r = Alum(1);
        Ebr = Alum(2);
        tan_d = Alum(3);
        sigma_d = Alum(4);
    case 3
        e_r = Bar(1);
        Ebr = Bar(2);
        tan_d = Bar(3);
        sigma_d = Bar(4);
        
    case 4
        e_r = Glass(1);
        Ebr = Glass(2);
        tan_d = Glass(3);
        sigma_d = Glass(4); 
       
    case 5
        e_r = Ice(1);
        Ebr = Ice(2);
        tan_d = Ice(3);
        sigma_d = Ice(4); 
        
    case 6
        e_r = Mica(1);
        Ebr = Mica(2);
        tan_d = Mica(3);
        sigma_d = Mica(4); 
        
    case 7
        e_r = Poleth(1);
        Ebr = Poleth(2);
        tan_d = Poleth(3);
        sigma_d = Poleth(4); 
        
    case 8 
        e_r = Polsty(1);
        Ebr = Polsty(2);
        tan_d = Polsty(3);
        sigma_d = Polsty(4); 
        
    case 9
        e_r = Quartz(1);
        Ebr = Quartz(2);
        tan_d = Quartz(3);
        sigma_d = Quartz(4); 
        
    case 10
        e_r = Si(1);
        Ebr = Si(2);
        tan_d = Si(3);
        sigma_d = Si(4); 
        
    case 11
        e_r = Soil(1);
        Ebr = Soil(2);
        tan_d = Soil(3);
        sigma_d = Soil(4); 
        
    case 12
        e_r = Teflon(1);
        Ebr = Teflon(2);
        tan_d = Teflon(3);
        sigma_d = Teflon(4); 
        
    case 13
        e_r = Water(1);
        Ebr = Water(2);
        tan_d = Water(3);
        sigma_d = Water(4); 
        
    case 14
        e_r = Seawater(1);
        Ebr = Seawater(2);
        tan_d = Seawater(3);
        sigma_d = Seawater(4); 
        
    case 15
        e_r = input('Please enter a value for the relative permitivity of the dielectric: ');
        Ebr = input('Please enter a value for the Electric field breakdown in V/m: ');
        tan_d = input('Please enter a value for the tangent of the dielectric phase at 1MHz: ');
        sigma_d = input('Please enter a value for the permiability of the dielectric in S/m: ');
        fprintf(['You chose Other with a relative permitivity of ', num2str(e_r), ', an Electric field breakdown of ',num2str(Ebr), ...
            ',\na tangent of the dielectric phase of ', num2str(tan_d), ', and a permiability of ', num2str(sigma_d), '.\n']);
        
end
      
%% asks the user to enter several design requirements
incorrect = true;
fprintf('----------------------------------------------------\n\n')
a = input('Please enter a value for the inner radius in mm: ');

%Replace with Warning instead 
%Make outer Diameter to Dielectric thickness and calcualte outer Diameter
%Have fixed values for Tickness, Length, and Voltage and give that info
%%when used


while incorrect
    b = input('Please enter a value for the outer radius in mm: ');
    if b > a
        incorrect = false;
       break
    end 
     disp('Error, outer radius must be bigger than the inner radius.');
end
f = input('Please enter a value for the operating frequency in Hz: ');
z_s = 0; % source impedance is assumed to be 0
while true
    l = input('Please enter the length in m: ');
    if l > 0
        break
    end
    disp('Error, Lenght needs to be greater than 0.')
end
while true
    v_s = input('Please enter the operating voltage in V: ');
    if v_s ~= 0
        break
    end
    disp('Error, input voltage can not be 0.')
end

fprintf('----------------------------------------------------\n\n')

%% lossy or lossless system
msg_loss = "Is your system lossless or lossy?"; % ask the user if there is loss
opts_loss = ["Lossy" "Lossless"]; % options
choice_loss = menu(msg_loss, opts_loss); % create interactive menu

% ***LOSSY LOSSY LOSSY LOSSY LOSSY LOSSY LOSSY***
if choice_loss == 1 % system is LOSSY
    disp('The system is now lossy!')
    check_lossy = input('Is this correct? [Y or N]: ', 's'); % confirm with user of choice

    if check_lossy == "y" || check_lossy == "Y"
        % load check
        msg_load_check = "Do you plan to connect a load?"; % ask the user if there is a load 
        opts_load_check = ["Yes" "No"]; % options
        choice_load_check = menu(msg_load_check, opts_load_check); % create interactive menu
        
        if choice_load_check == 1 % there IS a load 
            fprintf('----------------------------------------------------\n\n')
            disp('The system has a load!')
            z_l = input('Please enter a value for the load impedance in Ohms: ');
            check_load = input('Is this correct? [Y or N]: ', "s"); % confirm with user of choice
            fprintf('----------------------------------------------------\n\n')
        
            if check_load == "y" || check_load == "Y"
                w = 2 * pi * f; % angular freq
                omega = 2 * pi .* sweep;
                j = sqrt(-1); % imaginary number
                G = 2 * pi * sigma_d / log(b / a); % conductance
                C = 2 * pi * e_r * e_0 / log(b / a); % capacitance
                L = u_r * u_0 * log(b / a) / (2 * pi); % inductance
                R = (1 / (2 * pi)) * ((1 / a) +(1 / b)) * sqrt(pi * f * u_0 / sigma_c); % resistance
                R_sweep = (1 / (2 * pi)) * ((1 / a) +(1 / b)) * sqrt(pi .* sweep * u_0 / sigma_c);
                gamma = sqrt((R + j * w * L) * (G + j * w * C)); % propagation constant
                gamma_sweep = sqrt((R_sweep + j .* omega .* L) .* (G + j .* omega .* C));
                z_0 = (R + j * w * L) / gamma; % characteristic impednace 
                alpha = real(gamma); % attenuation constant
                Asweep = real(gamma_sweep);
                beta = imag(gamma); % phase constant
                wl = 2 * pi / beta; % wavelength
                u_p = w / beta; % propagation velocity
                G_dB = 10 * log10(exp(-2 * alpha * l)); % gain attenuation in dB
                reflect = (z_l - z_0) / (z_l + z_0); % reflection coefficient
                z_in = z_0 * (z_l + z_0 * tanh(gamma * l)) / (z_0 + z_l * tanh(gamma * l)); % input impedance
                v_in = v_s * (z_in / (z_in + z_s)); % input voltage
                v_0p = v_in / (exp(gamma * l) + reflect*exp(-gamma * l)); % received voltage 
                v_0n = reflect * v_0p; % reflected voltage
                v_l = v_0p + v_0n; % load voltage
                VSWR = (1 + abs(reflect)) / (1 - abs(reflect)); % voltage standing wave ratio
                
                % display
                disp(['alpha = ', num2str(alpha), ' Np/m'])
                disp(['beta = ', num2str(beta), ' rad/m'])
                disp(['Z_0 = ', num2str(z_0), ' Ohms'])
                disp(['VSWR = ', num2str(VSWR)])
                disp(['G_db = ', num2str(G_dB), ' dB'])
                disp(['Z_in = ', num2str(z_in), ' Ohms'])
                disp(['V_in = ', num2str(v_in), ' V'])
                disp(['V_L = ', num2str(v_l), ' V'])
                disp(['G_db = ', num2str(G_dB), ' dB'])

                %display the frequency response
                sys=tf([C*l 1], [R*C*l 1]);
                bode(sys)
                
                %display the frequency sweep
                figure(2);
                plot(sweep, Asweep);
                title('Attenuation vs. Frequency');
                xlabel('Frequency (Hz)');
                ylabel('Attenuation (Np/m)');

            else
                disp('You will have to restart the program.');
                return;
               
            end
        end

        if choice_load_check == 2 % there ISNT a load 
            fprintf('----------------------------------------------------\n\n')
            disp('The system does not have a load!')
            check_load = input('Is this correct? [Y or N]: ', "s"); % confirm with user of choice
            fprintf('----------------------------------------------------\n\n')

            if check_load == "y" || check_load == "Y"
                w = 2 * pi * f; % angular freq
                j = sqrt(-1); % imaginary number
                G = 2 * pi * sigma_d / log(b / a); % conductance
                C = 2 * pi * e_r * e_0 / log(b / a); % capacitance
                L = u_r * u_0 * log(b / a) / (2 * pi); % inductance
                R = (1 / (2 * pi)) * ((1 / a) +(1 / b)) * sqrt(pi * f * u_0 / sigma_c); % resistance
                gamma = sqrt((R + j * w * L) * (G + j * w * C)); % propagation constant
                z_0 = (R + j * w * L) / gamma; % characteristic impednace 
                alpha = real(gamma); % attenuation constant
                omega = 2 * pi .* sweep;
                R_sweep = (1 / (2 * pi)) * ((1 / a) +(1 / b)) * sqrt(pi .* sweep * u_0 / sigma_c);
                gamma_sweep = sqrt((R_sweep + j .* omega .* L) .* (G + j .* omega .* C));
                Asweep = real(gamma_sweep);
                beta = imag(gamma); % phase constant
                wl = 2 * pi / beta; % wavelength
                u_p = w / beta; % propagation velocity
                G_dB = 10 * log10(exp(-2 * alpha * l)); % gain attenuation in dB

                % display
                disp(['alpha = ', num2str(alpha), ' Np/m'])
                disp(['beta = ', num2str(beta), ' rad/m'])
                disp(['Z_0 = ', num2str(z_0), ' Ohms'])
                disp(['G_db = ', num2str(G_dB), ' dB'])
                
                % display the frequency response
                sys=tf([C*l 1], [R*C*l 1]);
                bode(sys)
                
                %display the frequency sweep
                figure(2);
                plot(sweep, Asweep);
                title('Attenuation vs. Frequency');
                xlabel('Frequency (Hz)');
                ylabel('Attenuation (Np/m)');

            else
                disp('You will have to restart the program.');
                return;
            end
        end 

    else
    disp('You will have to restart the program.');
    return;
   
    end
end

% ***LOSSLESS LOSSLESS LOSSLESS LOSSLESS LOSSLESS***
if choice_loss == 2 % system is LOSSLESS
    disp('The system is now lossless!')
    check_lossless = input('Is this correct? [Y or N]: ', "s"); % confirm with user of choice

    if check_lossless == "y" || check_lossless == "Y"
        msg_load_check = "Do you plan to connect a load?"; % ask the user if there is a load 
        opts_load_check = ["Yes" "No"]; % options
        choice_load_check = menu(msg_load_check, opts_load_check); % create interactive menu
        
        if choice_load_check == 1 % there IS a load 
            fprintf('----------------------------------------------------\n\n')
            disp('The system has a load!')
            z_l = input('Please enter a value for the load impedance in Ohms: ');
            check_load = input('Is this correct? [Y or N]: ', "s"); % confirm with user of choice
            fprintf('----------------------------------------------------\n\n')
        
            if check_load == "y" || check_load == "Y"
                w = 2 * pi * f; % angular freq
                j = sqrt(-1); % imaginary number
                G = 2 * pi * sigma_d / log(b / a); % conductance
                C = 2 * pi * e_r * e_0 / log(b / a); % capacitance
                L = u_r * u_0 * log(b / a) / (2 * pi); % inductance
                R = (1 / (2 * pi)) * ((1 / a) +(1 / b)) * sqrt(pi * f * u_0 / sigma_c); % resistance
                beta = w * sqrt(L * C); % phase constant
                gamma = j * beta; % propagation constant
                omega = 2 * pi .* sweep;
                R_sweep = (1 / (2 * pi)) * ((1 / a) +(1 / b)) * sqrt(pi .* sweep * u_0 / sigma_c);
                gamma_sweep = sqrt((R_sweep + j .* omega .* L) .* (G + j .* omega .* C));
                Asweep = real(gamma_sweep);
                z_0 = sqrt(L / C); % characteristic impedance
                wl = 2 * pi / beta; % wavelength
                u_p = w / beta; % propagation velocity

                reflect = (z_l - z_0) / (z_l + z_0); % reflection coefficient
                z_in = z_0 * (z_l + j * z_0 * tan(beta * l)) / (z_0 + j * z_l * tan(beta * l)); % input impedance
                v_in = v_s * (z_in / (z_in + z_s)); % input voltage
                v_0p = v_in / (exp(gamma * l) + reflect * exp(-gamma * l)); % received voltage 
                v_0n = reflect * v_0p; % reflected voltage
                v_l = v_0p + v_0n; % load voltage
                VSWR = (1 + abs(reflect)) / (1 - abs(reflect)); % voltage standing wave ratio

                % display
                disp(['gamma = ', num2str(gamma), ' rad/m'])
                disp(['Z_0 = ', num2str(z_0), ' Ohms'])
                disp(['Reflection coefficient = ', num2str(reflect)])
                disp(['VSWR = ', num2str(VSWR)])
                disp(['Z_in = ', num2str(z_in), ' Ohms'])
                disp(['V_in = ', num2str(v_in), ' V'])
                disp(['V_L = ', num2str(v_l), ' V'])

                % display the frequency response
                sys=tf([C*l 1], [R*C*l 1]);
                bode(sys)
                
                %display the frequency sweep
                figure(2);
                plot(sweep, Asweep);
                title('Attenuation vs. Frequency');
                xlabel('Frequency (Hz)');
                ylabel('Attenuation (Np/m)');
                
            else
                disp('You will have to restart the program.');
                return;
               
            end
        end

        if choice_load_check == 2 % there ISNT a load 
            fprintf('----------------------------------------------------\n\n')
            disp('The system does not have a load!')
            check_load = input('Is this correct? [Y or N]: ', "s"); % confirm with user of choice
            fprintf('----------------------------------------------------\n\n')

            if check_load == "y" || check_load == "Y"
                w = 2 * pi * f; % angular freq
                j = sqrt(-1); % imaginary number
                G = 2 * pi * sigma_d / log(b / a); % conductance
                C = 2 * pi * e_r * e_0 / log(b / a); % capacitance
                L = u_r * u_0 * log(b / a) / (2 * pi); % inductance
                R = (1 / (2 * pi)) * ((1 / a) +(1 / b)) * sqrt(pi * f * u_0 / sigma_c); % resistance
                beta = w * sqrt(L * C); % phase constant
                gamma = j * beta; % propagation constant
                omega = 2 * pi .* sweep;
                R_sweep = (1 / (2 * pi)) * ((1 / a) +(1 / b)) * sqrt(pi .* sweep * u_0 / sigma_c);
                gamma_sweep = sqrt((R_sweep + j .* omega .* L) .* (G + j .* omega .* C));
                Asweep = real(gamma_sweep);
                z_0 = sqrt(L / C); % characteristic impedance
                wl = 2 * pi / beta; % wavelength
                u_p = w / beta; % propagation velocity

                % display
                disp(['gamma = ', num2str(gamma), ' rad/m'])
                disp(['Z_0 = ', num2str(z_0), ' Ohms'])
                % display the frequency response
                sys=tf([C*l 1], [R*C*l 1]);
                bode(sys)
                
                %display the frequency sweep
                figure(2);
                plot(sweep, Asweep);
                title('Attenuation vs. Frequency');
                xlabel('Frequency (Hz)');
                ylabel('Attenuation (Np/m)');

            else
                disp('You will have to restart the program.');
                return;
            end
        end 

    else
        disp('You will have to restart the program.');
        return;
       
    end
end
