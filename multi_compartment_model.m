function multi_compartment_model()
    % Define parameters and initial conditions
    k12 = 0.1;   % Transfer rate constant from compartment 1 to compartment 2
    k21 = 0.05;  % Transfer rate constant from compartment 2 to compartment 1
    k23 = 0.2;   % Transfer rate constant from compartment 2 to compartment 3
    k32 = 0.1;   % Transfer rate constant from compartment 3 to compartment 2
    k34 = 0.15;  % Transfer rate constant from compartment 3 to compartment 4
    k43 = 0.1;   % Transfer rate constant from compartment 4 to compartment 3
    k_elim = 0.05; % Elimination rate constant
    
    C1_0 = 100;   % Initial concentration in compartment 1 (mg/mL)
    C2_0 = 0;    % Initial concentration in compartment 2 (mg/mL)
    C3_0 = 0;    % Initial concentration in compartment 3 (mg/mL)
    C4_0 = 0;    % Initial concentration in compartment 4 (mg/mL)
    
    initial_conditions = [C1_0, C2_0, C3_0, C4_0];
    
    % Time span for simulation
    tspan = [0, 24]; % Define the time span (0 to 24 hours)
    
    % Solve the system of ODEs using ODE45
    [t, C] = ode45(@ode_system, tspan, initial_conditions);

    % Calculate the total concentration in the entire body
    total_concentration = sum(C, 2);
    
    % Plot the concentration profiles
    figure;
    plot(t, C(:, 1), 'r', 'LineWidth', 1.5);
    hold on;
    plot(t, C(:, 2), 'g', 'LineWidth', 1.5);
    plot(t, C(:, 3), 'b', 'LineWidth', 1.5);
    plot(t, C(:, 4), 'm', 'LineWidth', 1.5);
    plot(t, total_concentration, 'k--', 'LineWidth', 2); % Total concentration

    xlabel('Time (hours)');
    ylabel('Concentration (mg/mL)');
    xlim([tspan(1), tspan(2)]);
    legend('Compartment 1', 'Compartment 2', 'Compartment 3', 'Compartment 4', 'Total Concentration');
    title('Multi-Compartment Pharmacokinetic Model');
    grid on;

    hold off;

    % Calculate the log of concentration
    log_concentration = log(C);
    log_total_concentration = log(total_concentration);
    
    figure;
    plot(t, log_concentration(:, 1), 'r', 'LineWidth', 2);
    hold on;
    plot(t, log_concentration(:, 2), 'g', 'LineWidth', 2);
    plot(t, log_concentration(:, 3), 'b', 'LineWidth', 2);
    plot(t, log_concentration(:, 4), 'm', 'LineWidth', 2);
    plot(t, log_total_concentration, 'k--', 'LineWidth', 2); % Log total concentration
    xlabel('Time (hours)');
    ylabel('Log Concentration (ln[mg/mL])');
    xlim([tspan(1), tspan(2)]);
    legend('Compartment 1', 'Compartment 2', 'Compartment 3', 'Compartment 4', 'Log Total Concentration');
    grid on;

    hold off;

    % Plot the concentrations vs. the log of time using semilogx
    figure;
    semilogx(t, C(:, 1), 'r', 'LineWidth', 2);
    hold on;
    semilogx(t, C(:, 2), 'g', 'LineWidth', 2);
    semilogx(t, C(:, 3), 'b', 'LineWidth', 2);
    semilogx(t, C(:, 4), 'm', 'LineWidth', 2);
    semilogx(t, total_concentration, 'k--', 'LineWidth', 2); % Total concentration
    xlabel('Log Time (ln[hours])');
    ylabel('Concentration (mg/mL)');
    xlim([tspan(1), tspan(2)]);
    legend('Compartment 1', 'Compartment 2', 'Compartment 3', 'Compartment 4', 'Total Concentration');
    title('Concentration vs. Log Time');
    grid on;
    hold off;
    
    function dCdt = ode_system(t, C)
        % System of ODEs representing the multi-compartment model
        C1 = C(1);
        C2 = C(2);
        C3 = C(3);
        C4 = C(4);
        
        dC1dt = k21 * C2 - k12 * C1 - k_elim * C1;
        dC2dt = k12 * C1 - k21 * C2 - k23 * C2 + k32 * C3;
        dC3dt = k23 * C2 - k32 * C3 - k34 * C3 + k43 * C4;
        dC4dt = k34 * C3 - k43 * C4 - k_elim * C4;
        
        dCdt = [dC1dt; dC2dt; dC3dt; dC4dt];
    end
end
