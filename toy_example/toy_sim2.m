N = 3; % number of agents
L = ringLaplacian(N); % compute the ring Laplacian matrix
k = 1; % coupling gain
n = 2; % system order
p = 3; % number of parameter

plotskip = 100;        % Interval for recording data points for plotting
dt       = 1e-3;       % Time step
Steps = 200;
totalSteps = plotskip * Steps;  % Total number of simulation steps 

% Initial states and parameters
x0 = rand(n*N,1);
theta10 = 0.9 + (1:N)'/10;
theta20 = 1 - (1:N)'/10;
theta30 = 1 - (1:N)'/10;

% Define psi function
psi = @(x,t) [ -x(1), sin(t), 1; sin(2*t), -x(2), 1 ];

% Options
flag_coupling = 1;  % 1 or 0
flag_adapt = 1;     % 1 or 0

% Reset states and parameters for each simulation mode
x = x0;
theta = [theta10, theta20, theta30]';
theta = theta(:);

count = 0;
rec_t = zeros(1,Steps);
rec_x = zeros(n*N,Steps);
rec_theta = zeros(p*N,Steps);

% Simulation loop
while count < totalSteps
    % Compute state derivatives
    t = count*dt;
    Psi = zeros(n*N,p*N);
    for i = 1:N
        Psi((i-1)*n+1:i*n, (i-1)*p+1:i*p) = psi(x((i-1)*n+1:i*n),t);
    end
    dx = Psi * theta - (flag_coupling) * k * kron(L,eye(2)) * x;
    dtheta = (flag_adapt) * (1/sqrt(k)) * Psi' * (-k*kron(L,eye(2))*x);

    % Update states using Euler method
    x = x + dx * dt;
    theta = theta + dtheta * dt;

    count = count + 1;

    % Record data at specified intervals
    if mod(count, plotskip) == 0
        index = count/plotskip;
        rec_t(index) = t;
        rec_x(:,index) = x;
        rec_theta(:,index) = theta;
    end
end

figure(1)
subplot(1, 2, 1)    
plot(rec_t, rec_x(1:2:end,:), 'LineWidth', 1.5)
title('x1')
subplot(1, 2, 2)    
plot(rec_t, rec_x(2:2:end,:), 'LineWidth', 1.5)
title('x2')

figure(2)
subplot(1, 3, 1)    
plot(rec_t, rec_theta(1:3:end,:), 'LineWidth', 1.5)
title('theta1')
subplot(1, 3, 2)    
plot(rec_t, rec_theta(2:3:end,:), 'LineWidth', 1.5)
title('theta2')
subplot(1, 3, 3)    
plot(rec_t, rec_theta(3:3:end,:), 'LineWidth', 1.5)
title('theta3')



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