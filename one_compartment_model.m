function one_compartment_model()
    % Define parameters and initial conditions
    k_elim = 0.05; % Elimination rate constant
    
    V1 = 100;    % Volume of compartment 1 (mL)
    
    C1_0 = 10;   % Initial concentration in compartment 1 (mg/mL)
    
    initial_conditions = C1_0;
    
    % Time span for simulation (changed to 24 hours)
    tspan = [0, 24]; % Define the time span (0 to 24 hours)
    
    % Solve the system of ODEs using ODE45
    [t, C1] = ode45(@ode_system, tspan, initial_conditions);
    
    % Plot the concentration profile
    figure;
    plot(t, C1, 'black--', 'LineWidth', 2);
    xlabel('Time (hours)');
    ylabel('Concentration (mg/mL)');
    title('One-Compartment Pharmacokinetic Model');
    grid on;
    
    function dCdt = ode_system(t, C1)
        % System of ODEs representing the one-compartment model
        dC1dt = -k_elim * C1;
        dCdt = dC1dt;
    end
end