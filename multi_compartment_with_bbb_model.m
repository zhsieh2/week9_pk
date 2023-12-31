function multi_compartment_with_bbb_model
    % Define parameters and initial conditions
    k12 = 0.1;   % Transfer rate constant from compartment 1 to compartment 2
    k21 = 0.05;  % Transfer rate constant from compartment 2 to compartment 1
    k23 = 0.2;   % Transfer rate constant from compartment 2 to compartment 3
    k32 = 0.1;   % Transfer rate constant from compartment 3 to compartment 2
    k34 = 0.15;  % Transfer rate constant from compartment 3 to compartment 4
    k43 = 0.1;   % Transfer rate constant from compartment 4 to compartment 3
    k_elim = 0.05; % Common elimination rate constant
    k_bbb = 0.005;  % Transfer rate constant from compartment 3 (plasma) to the brain
    k_brain_elim = 0.01;  % Elimination rate constant for the brain
    
    V1 = 100;    % Volume of compartment 1 (well-perfused organs, mL)
    V2 = 200;    % Volume of compartment 2 (adipose, mL)
    V3 = 150;    % Volume of compartment 3 (plasma, mL)
    V4 = 50;     % Volume of compartment 4 (other tissues, excluding brain, mL)
    V_brain = 80;  % Volume of the brain compartment (mL)
    
    C1_0 = 10;   % Initial concentration in compartment 1 (mg/mL)
    C2_0 = 0;    % Initial concentration in compartment 2 (mg/mL)
    C3_0 = 0;    % Initial concentration in compartment 3 (mg/mL)
    C4_0 = 0;    % Initial concentration in compartment 4 (mg/mL)
    C_brain_0 = 0;  % Initial concentration in the brain compartment (mg/mL)
    
    initial_conditions = [C1_0, C2_0, C3_0, C4_0, C_brain_0];
    
    % Time span for simulation (changed to 24 hours)
    tspan = [0, 24]; % Define the time span (0 to 24 hours)
    
    % Solve the system of ODEs using ODE45
    [t, C] = ode45(@ode_system, tspan, initial_conditions);
    
    % Calculate the total concentration in the entire body
    total_concentration = sum(C(:, 1:4), 2);
    
    % Plot the concentration profiles, including the total concentration
    figure;
    plot(t, C(:, 1), 'r', 'LineWidth', 2);
    hold on;
    plot(t, C(:, 2), 'g', 'LineWidth', 2);
    plot(t, C(:, 3), 'b', 'LineWidth', 2);
    plot(t, C(:, 4), 'm', 'LineWidth', 2);
    plot(t, C(:, 5), '--', 'LineWidth', 2);

    plot(t, total_concentration, 'k--', 'LineWidth', 2); % Total concentration
    xlabel('Time (hours)');
    ylabel('Concentration (mg/mL)');
    legend('Compartment 1', 'Compartment 2', 'Compartment 3 (Plasma)', 'Compartment 4 (Tissues)', 'Brain', 'Total Concentration');
    title('Multi-Compartment Pharmacokinetic Model with BBB');
    grid on;
    
    function dCdt = ode_system(t, C)
        % System of ODEs representing the multi-compartment model with BBB
        C1 = C(1);
        C2 = C(2);
        C3 = C(3);
        C4 = C(4);
        C_brain = C(5);
        
        dC1dt = k21 * C2 - k12 * C1 - k_elim * C1;
        dC2dt = k12 * C1 - k21 * C2 - k23 * C2 + k32 * C3;
        dC3dt = k23 * C2 - k32 * C3 - k34 * C3 + k43 * C4 - k_bbb * C3;
        dC4dt = k34 * C3 - k43 * C4 - k_elim * C4;
        dC_brain_dt = k_bbb * C3 - k_brain_elim * C_brain; % Brain elimination
        
        dCdt = [dC1dt; dC2dt; dC3dt; dC4dt; dC_brain_dt];
    end
end
