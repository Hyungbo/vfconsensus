%% A simulation of coupled pendulums with timed mode switches and banners.

% Number of agents and compute the ring Laplacian matrix
N = 3;
L = ringLaplacian(N);
k = 1;

plotskip   = 100;      % Interval for recording data points for plotting/animation
dt         = 1e-3;     % Time step
T_end      = 50;       % ----- total simulation time (sec)
totalSteps = round(T_end/dt);

% ----- Mode switch times -----
T1 = 10;               % coupling starts at 10s
T2 = 20;               % adaptation starts at 20s

% ----- Banner settings -----
bannerDur = 1.5;       % seconds to show the banner

% ----- Parameters -----
g = 9.8;
m = 0.2;
b = 0.5;
A = 0.5;
epsilon = 0.1;

% Initial states and parameters
z0 = [0.1; 0; -0.1];   % Initial z state, N by 1
w0 = zeros(N,1);
w0(1:2) = [-0.5; 0.5]; % Initial w state, N by 1
l0 = 0.9 + (1:N)'/10;  % initial pendulum length
th10 = 1./l0.^2;       % N by 1
th20 = 1./l0;          % N by 1

% Define escape torque
escape_torque = @(theta,velocity) (abs(theta) < epsilon) .* sign(velocity) * A;

% Define psi function (t is replaced with z)
psi = @(w,z) [-(b/m)*(w-z) + (1/m)*escape_torque(z,w-z), -g*sin(z)];

% Reset states and parameters for each simulation mode
z   = z0;
w   = w0;
th1 = th10;
th2 = th20;

count    = 0;
rec_t    = [];
rec_z    = [];  
rec_w    = [];
rec_th1  = [];
rec_th2  = [];

flag_coupling = 0;
flag_adapt    = 0;

% =================== Simulation loop ===================
while count < totalSteps
    t = dt * count;

    % Timed mode switching
    if (t >= T1), flag_coupling = 1; end
    if (t >= T2), flag_adapt    = 1; end

    % Compute state derivatives
    dz = w - z;        % theta_dot = angular velocity
    P  = psi(w,z);     % N x 2
    Tm = [th1, th2]';  % 2 x N

    tmp = zeros(N,1);
    for l=1:N
        tmp(l) = P(l,:)*Tm(:,l);
    end

    dw   = tmp + (w - z) - (flag_coupling)*k*L*w;
    dth1 = (flag_adapt) * (sqrt(k)) * P(:,1) .* (-L*w);
    dth2 = (flag_adapt) * (sqrt(k)) * P(:,2) .* (-L*w);

    % Euler update
    z   = z   + dz   * dt;
    w   = w   + dw   * dt;
    th1 = th1 + dth1 * dt;
    th2 = th2 + dth2 * dt;

    % Record data at specified intervals
    if mod(count, plotskip) == 0
        rec_t   = [rec_t, t];
        rec_w   = [rec_w, w];
        rec_z   = [rec_z, z];
        rec_th1 = [rec_th1, th1];
        rec_th2 = [rec_th2, th2];
    end

    count = count + 1;
end

% =================== Quick plots ===================
figure;
subplot(1, 2, 1)
plot(rec_t, rec_z, 'LineWidth', 1.5)
title('z (angle)')
xlabel('t (s)')
subplot(1, 2, 2)
plot(rec_t, rec_th1, 'LineWidth', 1.5)
title('\theta_1 parameters')
xlabel('t (s)')

% =================== Animation ===================
hFig = figure;
axis equal;
axis([-1.2 5.2 -2.5 0.2]);
hold on;
grid on;
title('Pendulum Test');
xlabel('x (m)');
ylabel('y (m)');

l = 1./rec_th1;   % lengths from recorded th1
z = rec_z;        % angles from recorded z

% Initialize graphics objects
pendulum_line = gobjects(N,1);
pendulum_mass = gobjects(N,1);
for i=1:N
    pendulum_line(i) = plot([2*(i-1), 2*(i-1) + l(i,1) * sin(z(i,1))], ...
                            [0, -l(i,1) * cos(z(i,1))], 'k-', 'LineWidth', 2);
    pendulum_mass(i) = plot(2*(i-1) + l(i,1) * sin(z(i,1)), ...
                            -l(i,1) * cos(z(i,1)), ...
                            'ro', 'MarkerSize', 10, 'MarkerFaceColor', 'r');
end
h_time   = text(-1, -2.3, 't = 0.00');
h_banner = [];  % will hold banner text handle

% ----- VideoWriter (MP4) -----
v = VideoWriter('pendulum_animation.mp4','MPEG-4');
v.FrameRate = 30;   % frames per second
v.Quality   = 95;
open(v);

for kf = 1:length(rec_t)
    % Update pendulums
    for i=1:N
        set(pendulum_line(i), 'XData', [2*(i-1), 2*(i-1) + l(i,kf) * sin(z(i,kf))], ...
                              'YData', [0, -l(i,kf) * cos(z(i,kf))]);
        set(pendulum_mass(i), 'XData', 2*(i-1) + l(i,kf) * sin(z(i,kf)), ...
                              'YData', -l(i,kf) * cos(z(i,kf)));
    end
    set(h_time, 'String', sprintf('t = %.2f s', rec_t(kf)) );

    % ----- Show banners around switch times -----
    tnow = rec_t(kf);
    showCoupling   = (tnow >= T1) && (tnow <= T1 + bannerDur);
    showAdaptation = (tnow >= T2) && (tnow <= T2 + bannerDur);

    % Delete previous banner (if any)
    if ~isempty(h_banner) && isvalid(h_banner)
        delete(h_banner);
        h_banner = [];
    end

    if showCoupling || showAdaptation
        % center of current axes
        xl = xlim; yl = ylim;
        xc = mean(xl); yc = mean(yl);
        if showCoupling
            bannerText = 'coupling starts';
        else
            bannerText = 'adaptation starts';
        end
        % Draw a prominent banner
        % Note: If your MATLAB version doesn't support RGBA in BackgroundColor,
        %       change [1 1 1 0.7] to [1 1 1].
        try
            h_banner = text(xc, yc, bannerText, ...
                'HorizontalAlignment','center', ...
                'VerticalAlignment','middle', ...
                'FontSize', 22, 'FontWeight','bold', ...
                'BackgroundColor',[1 1 1 0.7], ...
                'Margin',10);
        catch
            h_banner = text(xc, yc, bannerText, ...
                'HorizontalAlignment','center', ...
                'VerticalAlignment','middle', ...
                'FontSize', 22, 'FontWeight','bold', ...
                'BackgroundColor',[1 1 1], ...
                'Margin',10);
        end
    end

    drawnow;

    % Capture and write one frame per recorded time step
    frame = getframe(hFig);
    writeVideo(v, frame);
end

% Finish video
close(v);

%% --------------- Helper: Laplacian of ring graph ---------------
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
