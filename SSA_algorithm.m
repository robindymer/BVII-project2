
x = zeros(10e7, 9);
t = zeros(10e7, 1);
j = 1;
x(j,:) = [0 0 1 0 1 0 0 0 0]; % initial concentration
t(j) = 0;
Tf = 100; % final time

% Parameter that is subject to change (is tested for 0.2 and 0.08)
delta_R = 0.08;
% Parameters for the system in the following order:
%   alfa_A  alfa'_A  alfa_R  alfa'_R  beta_A  beta_R  teta_A  teta_R ...
%   gamma_A  gamma_R  gamma_C  delta_M_R  delta_M_A  delta_A  delta_R
params = [50, 500, 0.01, 50, 50, 5, 50, 100, 1, 1, 2, 0.5, 10, 1, delta_R];

while t(j,1) < Tf
    w = prop_vilar(x(j,:), params);
    a0 = sum(w);
    % Generate tau
    u = rand(1);
    tau = -log(u)/a0;
    % Generate the reaction that should occur
    v = rand(1);
    p = w / a0;
    F = cumsum(p); % the cumulative probability function
    r = find(F > v, 1);
    
    % Update
    j = j + 1;
    nr = nr_vilar();
    x(j,:) = x(j-1,:) + nr(r,:);
    t(j) = t(j-1) + tau;
end
% Since we (possibly) have a tail of zeros,
% cleanse the matrices
t_clean = t(1:j);
% Activator protein A
x_clean = x(1:j,1);
% Represant protein R
x2_clean = x(1:j,9);
hold on
title("Gillespie's algorithm, dR = 0.08")
ylabel("Concentration")
xlabel("Time")
plot(t_clean, x_clean, 'b-');
plot(t_clean, x2_clean, 'r-');
legend('Activator A', 'Repressor R')
hold off