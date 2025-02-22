% During the animation, hit the space bar to switch from uncoupled
% -> coupled -> coupled + adaptation

clear all
clc
close all

global SimStep
SimStep = 0;

rng('default')

K = 1;
plotskip = 2000;
dt = 0.0001;
N = 5;                  % number of agents
L = ringLaplacian(N);

mu = [1:N]';           % initial parameters
nu = [1:N]';           % initial parameters
z0 = rand(N,1)-0.5;     % random initial condition
y0 = rand(N,1)-0.5;     % random initial condition
z = z0;
y = y0;

k = 0;
flag_adapt = 0;
display_str = [ "No interaction", "Coupled", "Adaptation" ];

Cstr = [ "ro", "go", "bo", "mo", "ko" ];
Cchr = [ "r", "g", "b", "m", "k" ];

hFig = figure;
set(hFig,'WindowKeyPressFcn',@keyPressCallback);
axis([-3,3,-12,12])
grid on
hold on

% Save Avi file
%v = VideoWriter('animation.avi');
%v.FrameRate = 5;  % 5fps
%open(v);

count = 0;
while true
    
    % define derivatives
    dz = -z + y;
    dy = ( ones(N,1) - mu.*(z.^2 - ones(N,1)) ).*(-z + y) - nu.*z - k*L*y;
    if flag_adapt == 1
        dmu = (ones(N,1) - z.^2) .* (-z + y) .* (-k*L*y);
        dnu = -nu .* z .* (-k*L*y);
    end

    % update states 
    z = z + dz*dt;
    y = y + dy*dt;
    if flag_adapt == 1
        mu = mu + dmu*dt;
        nu = nu + dnu*dt;
    end
    count = count + 1;
    
    switch SimStep
        case 1
            k = K;
        case 2
            flag_adapt = 1;
        case 3
            close;
            break;
    end
    
    % draw and save frame
    if mod(count,plotskip) == 0
        t1 = text(0.8, 10, sprintf('Time: %3.2f / %s', count*dt, display_str(SimStep+1) ));
        for j = 1:N
            p(j) = plot(z(j), y(j), Cstr(j), 'MarkerSize', 10, 'MarkerFaceColor', Cchr(j));
        end
        drawnow;
        
        % Save avi file
        %frame = getframe(hFig);
        %writeVideo(v, frame);
        
        pause(0.1);
        delete(t1);
        for j = 1:N
            delete(p(j))
        end
    end
end

% Save avi file
%close(v);

function keyPressCallback(source, eventdata)
    global SimStep
    % Hit space bar, increase SimStep
    keyPressed = eventdata.Key;
    if strcmpi(keyPressed,'space')
        SimStep = SimStep + 1;
    end
end 

function L = ringLaplacian(n)
    % ringLaplacian - Returns the Laplacian matrix of a ring (cycle) graph with n nodes.
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
        j = mod(i, n) + 1;  % Next node (with wrap-around)
        A(i, j) = 1;
        A(j, i) = 1;        % undirected graph
    end

    % Degree matrix: each node has degree 2 in a ring graph (for n>=3)
    D = diag(sum(A, 2));

    % Laplacian matrix: L = D - A
    L = D - A;
end
