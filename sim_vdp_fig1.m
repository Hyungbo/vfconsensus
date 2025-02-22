clear all
clc
close all

global SimStep
SimStep = 0;

%rng('default')

K = .5;
plotskip = 2000;
dt= 0.0001;
N = 3;                  % number of agents
L = ringLaplacian(N);

mu = [1:N]';           % initial parameters
nu = [1:N]';   % initial parameters
%z0 = 10*(rand(N,1)-0.5);     % random initial condition
%y0 = 10*(rand(N,1)-0.5);     % random initial condition
z0 = [1; 1; 1];
y0 = [0; -2; 2];
z = z0;
y = y0;

k = 0;
flag_adapt = 0;
display_str = [ "No interaction", "Coupled", "Adaptation" ];


Cstr = [ "ro", "go", "bo", "mo", "ko" ];
Cchr = [ "r", "g", "b", "m", "k" ];


figure
for SimStep = 0:2
   switch SimStep
        case 1
            k = K;
        case 2
            flag_adapt = 1;
    end

    z = z0;
    y = y0;
    count = 0;
    rec_x = [];
    rec_t = [];

    while count < plotskip * 60
        
        % define derivatives
        dz = -z+y;
        dy = ( ones(N,1) - mu.*(z.^2-ones(N,1)) ).*(-z+y) - nu.*z - k*L*y;
        if flag_adapt == 1
            dmu = 0.05*( (ones(N,1) - z.^2) .* (-z+y) .* (-k*L*y) );
            dnu = 0.05*( -nu .* z .* (-k*L*y) );
        end
    
        % update states 
        z = z + dz*dt;
        y = y + dy*dt;
        if flag_adapt == 1
            mu = mu + dmu*dt;
            nu = nu + dnu*dt;
        end
        count = count + 1;
        
        
        % draw
        if mod(count,plotskip)==0
            
            rec_t = [rec_t, (count-plotskip/2)*dt];
            rec_x = [rec_x, y];
            
        end
    end

    subplot(3,1,SimStep+1)
    plot(rec_t,rec_x)
    switch SimStep
        case 0
            title('without couplings')
        case 1
            title('with couplings')
        case 2
            title('with couplings + adaptation')
    end

end



function L = ringLaplacian(n)
% ringLaplacian - Returns the Laplacian matrix of a ring (cycle) graph with n nodes.
%
% Syntax: L = ringLaplacian(n)
%
% Input:
%   n - Number of nodes in the ring graph.
%
% Output:
%   L - n-by-n Laplacian matrix of the ring graph.
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

    % Initialize the adjacency matrix for an n-node ring graph
    A = zeros(n);
    for i = 1:n
        % Connect node i to node j (with wrap-around)
        j = mod(i, n) + 1;  % Next node (wraps from n to 1)
        A(i, j) = 1;
        A(j, i) = 1;        % Since the graph is undirected
    end

    % Degree matrix: each node has degree 2 in a ring graph (for n>=3)
    D = diag(sum(A, 2));

    % Laplacian matrix: L = D - A
    L = D - A;
end

