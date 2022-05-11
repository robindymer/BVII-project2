
x = zeros(10e5, 9);
t = zeros(10e5, 1);
j = 1;
x(j,:) = [0 0 1 0 1 0 0 0 0]; % initial concentration
t(j) = 0;
Tf = 100; % final time

% parameters
alfa_A = 50;
alfaP_A = 500; % P for prime
alfa_R = 0.01; 
alfaP_R = 50; % P for prime
beta_A = 50;
beta_R = 5;
teta_A = 50;
teta_R = 100;
gamma_A = 1;
gamma_R = 1;
gamma_C = 2;
delta_M_R = 0.5;
delta_M_A = 10; 
delta_A = 1;
delta_R = 0.2; % This should be tested for 0.08 aswell

while t(j,1) < Tf
    w = prop_vilar(x(j,:), [alfa_A, alfaP_A, alfa_R, alfaP_R, beta_A, beta_R, teta_A, teta_R, gamma_A, gamma_R, gamma_C, delta_M_R, delta_M_A, delta_A, delta_R]);
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
% Since we (possibly) have a bunch of zeros at the tail,
% cleanse the matrices
t_clean = t(1:j);
% Activator protein A
x_clean = x(1:j,1);
% Represant protein R
x2_clean = x(1:j,9);
hold on
plot(t_clean, x_clean, 'bo');
plot(t_clean, x2_clean, 'ro');
hold off