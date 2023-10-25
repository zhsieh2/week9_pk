% Define pharmacokinetic parameters for each organ
% V: volume (mL), Q: blood flow rate (mL/min), CL: clearance (mL/min)
% Heart
parameters.heart = struct('V', 300, 'Q', 100, 'CL', 10); 
% Lungs
parameters.lungs = struct('V', 800, 'Q', 500, 'CL', 50);
% Brain
parameters.brain = struct('V', 1500, 'Q', 20, 'CL', 5);
% GI Tract
parameters.GI_tract = struct('V', 2000, 'Q', 100, 'CL', 80);
% Liver
parameters.liver = struct('V', 1200, 'Q', 150, 'CL', 60);
% Pancreas
parameters.pancreas = struct('V', 80, 'Q', 10, 'CL', 4);
% Kidneys
parameters.kidneys = struct('V', 260, 'Q', 250, 'CL', 70);
% Muscle
parameters.muscle = struct('V', 4500, 'Q', 400, 'CL', 200);
% Subcutaneous Fat
parameters.subcutaneous_fat = struct('V', 3000, 'Q', 50, 'CL', 20);
% Visceral Fat
parameters.visceral_fat = struct('V', 2000, 'Q', 30, 'CL', 15);

% Define drug administration parameters for oral route
dose = 50; % mg
bioavailability = 0.9;
% Define time span for the simulation
tspan = [0 24]; % hours
% Define initial conditions
C0 = zeros(1, length(fieldnames(parameters))); % Initial concentrations in each compartment
C0(1) = dose * bioavailability / parameters.heart.V; % Initial concentration in the heart
% Set up and solve the ODE system
[t, C] = ode45(@(t, y) multicompartment_model(t, y, parameters, dose, bioavailability, fieldnames(parameters)), tspan, C0);

% Set a predefined set of distinguishable colors
colors = lines(length(fieldnames(parameters)));

% Plot drug concentration over time for each organ on the same figure
figure;
hold on;
organ_names = fieldnames(parameters);
for i = 1:length(organ_names)
    plot(t, C(:, i), 'DisplayName', organ_names{i}, 'Color', colors(i, :));
end
hold off;
title('Drug Concentration in Various Organs');
xlabel('Time (hours)');
ylabel('Concentration (mg/L)');
legend;
grid on;

% Define the ODE system function
function dCdt = multicompartment_model(t, C, parameters, dose, bioavailability, organ_names)
    n = length(C);
    dCdt = zeros(n, 1);
    for i = 1:n
        if i == 1
            dCdt(i)=(dose*bioavailability/parameters.(organ_names{i}).V-(parameters.(organ_names{i}).CL/parameters.(organ_names{i}).V)*C(i));
        else
            dCdt(i)=(parameters.(organ_names{i}).Q/parameters.(organ_names{i}).V)*C(i-1)-(parameters.(organ_names{i}).CL/parameters.(organ_names{i}).V)*C(i);
        end
    end
end
