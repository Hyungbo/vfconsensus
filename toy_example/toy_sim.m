% Number of agents and compute the ring Laplacian matrix
N = 3;
L = ringLaplacian(N);
k = 1;

plotskip = 100;        % Interval for recording data points for plotting
dt       = 1e-3;         % Time step
Steps = 600;
totalSteps = plotskip * Steps;  % Total number of simulation steps 

% Initial states and parameters
x0 = rand(N,1);
theta10 = 0.9 + (1:N)'/10;
theta20 = 1 - (1:N)'/10;

% Define psi function (t is replaced with z)
psi1 = @(x,t) -x;
psi2 = @(x,t) sin(t)*ones(N,1);


% Options
flag_coupling = 1;  % 1 or 0
flag_adapt = 1;     % 1 or 0


% Reset states and parameters for each simulation mode
x  = x0;
theta1 = theta10;
theta2 = theta20;

count = 0;
rec_t = zeros(1,Steps);
rec_x = zeros(N,Steps);
rec_theta1 = zeros(N,Steps);
rec_theta2 = zeros(N,Steps);

% Simulation loop
while count < totalSteps
    % Compute state derivatives
    t = count*dt;
    Psi1 = psi1(x,t);
    Psi2 = psi2(x,t);
    dx = Psi1 .* theta1 + Psi2 .* theta2 - (flag_coupling)*k*L*x;
    dtheta1 = (flag_adapt) * (1/sqrt(k)) * Psi1 .* (-k*L*x);
    dtheta2 = (flag_adapt) * (1/sqrt(k)) * Psi2 .* (-k*L*x);

    % Update states using Euler method
    x = x + dx * dt;
    theta1 = theta1 + dtheta1 * dt;
    theta2 = theta2 + dtheta2 * dt;

    count = count + 1;

    % Record data at specified intervals
    if mod(count, plotskip) == 0
        index = count/plotskip;
        rec_t(index) = t;
        rec_x(:,index) = x;
        rec_theta1(:,index) = theta1;
        rec_theta2(:,index) = theta2;
    end
end

figure(1)
plot(rec_t, rec_x, 'LineWidth', 1.5)
title('x')

figure(2)
subplot(1, 2, 1)    
plot(rec_t, rec_theta1, 'LineWidth', 1.5)
title('theta1')
subplot(1, 2, 2)    
plot(rec_t, rec_theta2, 'LineWidth', 1.5)
title('theta2')



%%
function L = ringLaplacian(n)
% ringLaplacian - Returns the Laplacian matrix of a ring (cycle) graph with n nodes.
%
% Input:
%   n - Number of nodes in the ring graph.
%
% Output:
%   L - n-by-n Laplacian matrix.
%
% Special cases:
%   For n = 1, L is defined as 0.
%   For n = 2, L is [1 -1; -1 1].

    if n == 1
        L = 0;
        return;
    elseif n == 2
        L = [1 -1; -1 1];
        return;
    end

    % Create the adjacency matrix for an n-node ring graph
    A = zeros(n);
    for i = 1:n
        j = mod(i, n) + 1;  % Wrap-around index for the ring structure
        A(i, j) = 1;
        A(j, i) = 1;        % The graph is undirected
    end

    % Degree matrix: each node has degree 2 for n >= 3 in a ring graph
    D = diag(sum(A, 2));

    % Laplacian matrix: L = D - A
    L = D - A;
end