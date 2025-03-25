%% Initialization
clear; close all; clc;

% Simulation parameters
K        = 0.5;         %= 1;  % Coupling strength
plotskip = 2000;        % Interval for recording data points for plotting
dt       = 1e-4;         % Time step
totalSteps = plotskip * 90;  % Total number of simulation steps

% Number of agents and compute the ring Laplacian matrix
N = 3;
L = ringLaplacian(N);

% Initial states and parameters
mu0 = (1:N)';          % Initial mu values
nu0 = (1:N)';          % Initial nu values
z0  = [1; 1; 1];       % Initial z state
y0  = [0; -2; 2];      % Initial y state

% Simulation modes:
% 0 - No couplings, 1 - With couplings, 2 - With couplings and adaptation
modes = 0:2;

figure;
for mode = modes
    % Set coupling and adaptation based on the mode
    k = 0;
    flag_adapt = false;
    switch mode
        case 1
            k = K;
        case 2
            k = K;
            flag_adapt = true;
    end

    % Reset states and parameters for each simulation mode
    z  = z0;
    y  = y0;
    mu = mu0;
    nu = nu0;
    
    count = 0;
    rec_t = [];
    rec_y = [];  % Records of y values

    % Simulation loop
    while count < totalSteps
        % Compute state derivatives
        dz = -z + y;
        dy = (ones(N,1) - mu.*(z.^2 - 1)).*(-z + y) - nu.*z - k*L*y;
        
        if flag_adapt
            dmu = 0.05*( -(z.^2 - ones(N,1)) .* (-z+y) .* (-k*L*y) );
            dnu = 0.05*( -z .* (-k*L*y) );
        end
        
        % Update states using Euler method
        z = z + dz * dt;
        y = y + dy * dt;
        if flag_adapt
            mu = mu + dmu * dt;
            nu = nu + dnu * dt;
        end
        
        count = count + 1;
        
        % Record data at specified intervals
        if mod(count, plotskip) == 0
            rec_t = [rec_t, (count - plotskip/2)*dt];
            rec_y = [rec_y, y];
        end
    end
    
    % Plot results for the current mode
    subplot(3,1, mode+1)
    plot(rec_t, rec_y, 'LineWidth', 1.5)
    xlabel('Time'); ylabel('y');
    switch mode
        case 0
            title('Without Couplings')
        case 1
            title('With Couplings')
        case 2
            title('With Couplings + Adaptation')
    end
end

%% Function: ringLaplacian
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
