% Animating a mass spring system
% Author: Dr. David Montgomery, Colorado School of Mines
% Code tested in Matlab R2023b
% 
% Equation of Motion for mass spring system is:
% m*x'' + b*x' + k*x = F0*sin(gamma*t), 
% x(0) = P0, x'(0) = V0,
%
% where m is mass, b is the damping constant, k is the spring constant.
% The forcing function has amplitude F0 and forcing frequency gamma.
%
% Choose from the folowing examples:
% example = 'Freestyle' to edit parameters freely
%         = 'Undamped' for undamped free vibrations
%         = 'Underdamped' for underdamped free vibrations
%         = 'Overdamped' for overdamped free vibrations
%         = 'Critical' for critically damped free vibrations
%         = 'Forced' for underdamped forced vibrations
%         = 'Beats' for undamped forced vibrations w/beats phenomenom
%         = 'Resonance' for undamped forced vibrations w/resonance
% 
% If you want to save a video of the animation set saveVideo = 'true'

clear;
close all;

% -------------------------------------------------------------------------
% --- Edit the parameters between dashed lines ----------------------------
% -------------------------------------------------------------------------
example = 'beats';
saveVideo = 'false'; 

% Physical parameters (only if example = 'freestyle')
m = 1; % mass [kg]
b = 3; % damping constant [kg/s]
k = 1; % spring constant [kg/s^2]

% External forcing function F(t) = F0 * sin(gamma*t)
F0 = 1;
gamma = 0;

% Initial conditions
P0 = 1;
V0 = 0;

% Time parameters
dt = 0.25; % time step for ODE solver [s]
tf = 60; % Final time for simulation [s]

% Length of unstretched spring l, and length of stretching s
l = 3; % [m]
s = 2; % [m]

% Height/length of the mass (for visual purposes only)
lm = 2; % [m]

% -------------------------------------------------------------------------
% --- Don't edit below this line ------------------------------------------
% -------------------------------------------------------------------------

% These are the other possible examples
disp(' ')
disp(['example = ',example])
if strcmpi(example,'undamped')
    [m,b,k,F0,gamma,P0,V0,dt,tf] = undamped();
    titleName = 'Undamped Free Vibrations:  ';
elseif strcmpi(example,'underdamped') 
    [m,b,k,F0,gamma,P0,V0,dt,tf] = underdamped();
    titleName = 'Underdamped Free Vibrations:  ';
elseif strcmpi(example,'overdamped') 
    [m,b,k,F0,gamma,P0,V0,dt,tf] = overdamped();
    titleName = 'Overdamped Free Vibrations:  ';
elseif strcmpi(example,'critical') 
    [m,b,k,F0,gamma,P0,V0,dt,tf] = critical();
    titleName = 'Critically Damped Free Vibrations:  ';
elseif strcmpi(example,'forced') 
    [m,b,k,F0,gamma,P0,V0,dt,tf] = forced();
    titleName = 'Forced Vibrations Underdamped:  ';
elseif strcmpi(example,'beats') 
    [m,b,k,F0,gamma,P0,V0,dt,tf] = beats();
    titleName = 'Forced Beats:  ';
elseif strcmpi(example,'resonance') 
    [m,b,k,F0,gamma,P0,V0,dt,tf] = resonance();
    titleName = 'Forced Resonance:  ';
else
    titleName = '';
end

% Total length of spring after stretching
L = l + s; 

% Define the forcing function
F = @(t) F0*sin(gamma*t);

% Reduce order of DE and solve system of 1st order eqts numerically
% u1' = u2
% u2' = -2*lambda*u2 - omega^2*u1 + F(t)/m
omega = sqrt(k/m);
lambda = b/m/2;
f = @(t,u) [u(2); -2*lambda*u(2) - omega^2*u(1) + F(t)/m];
IC = [P0 V0];
tspan = 0:dt:tf;

[t,x] = ode45(f,tspan,IC);
x = x(:,1); % Get first column as u(:,1) = x(t)

% Function that provides visualization of spring
spr_d = 0;
if b ~= 0 % shift left if there is damping
    spr_d = -lm/4;
end
spring = @(y,spr_b,spr_s) 0.2*sin(spr_b*(y-spr_s))+spr_d;

% Define the video file name and settings
if strcmpi(saveVideo,'true')
    folderName = 'animations';
    videoFileName = [titleName(1:end-3)];
    % Check if the folder exists, if not, create it
    if ~exist(folderName, 'dir')
        mkdir(folderName);
    end
    
    % Define the full path for the video file
    videoFilePath = fullfile(folderName, [videoFileName, '.mp4']);
    
    % Create the VideoWriter object
    videoObj = VideoWriter(videoFilePath, 'MPEG-4');
    videoObj.FrameRate = 10; % Adjust the frame rate as needed
    videoObj.Quality = 75; % Adjust quality (0-100), 
    open(videoObj);

    % Figure Position settings
    % figPosition = [183 295 1253 533];

    % font size, marker size, etc
    fSize = 12;
    mSize = 8;
    lWidth = 2;
else 
    % Figure position if not saving video
    figPosition = [80 83 1356 745];
    fSize = 16;
    mSize = 20;
    lWidth = 4;
end

% Animate the motion
ax = figure();
if ~strcmpi(saveVideo,'true')
    set(ax, 'Position', figPosition);
end
ylim([-L+min(F(t)), max(x)+lm])
for i = 1:length(t)
    clf;

    subplot(1,2,1) % Spring animation on left
    % Display spring
    n = 6; % number of periods per length l
    spr_b = 2*pi*n/(L-F(t(i))+x(i)-lm/2); % frequency
    s = -L+F(t(i)); % phase shift
    spr_y = linspace(-L+F(t(i)),x(i),500);
    spr_x = spring(spr_y,spr_b,s);
    plot(spr_x,spr_y,'-k',LineWidth=2)
    hold on
    
    % Display damper if present
    if b ~= 0 % shift left if there is damping
        dw = 0.2;
        dc = lm/5;
        d_l = dc-dw*1.2;
        d_r = dc+dw*1.2;

        plot([d_l d_l],[-L+F(t(i)) x(i)],'-k',LineWidth=2)
        plot([d_r d_r],[-L+F(t(i)) x(i)],'-k',LineWidth=2)
        plot([dc dc],[x(i)-0.875*l x(i)],'-k',LineWidth=2)
        plot([dc-0.9*dw, dc+0.9*dw],[x(i)-0.875*l x(i)-0.875*l],'-k',LineWidth=2)
    end
    % plot equilibrium line y = 0
    plot([-2*lm,2*lm],[0 0],'k--')

    % plot fixed line
    plot([-2*lm,2*lm],[-L+F(t(i)) -L+F(t(i))],'k-',LineWidth=8)

    % plot upper and lower bounds for force
    if any(F(t))
        plot([-2*lm,2*lm],[-L+max(F(t)) -L+max(F(t))],'k:',LineWidth=1)
        plot([-2*lm,2*lm],[-L+min(F(t)) -L+min(F(t))],'k:',LineWidth=1)
    end

    % plot rectangle for mass
    rectangle('Position',[-lm,x(i)-lm/2,2*lm,lm],'FaceColor','k')
    
    % figure properties
    set(gca,'XTick',[],'YTick',[])
    ymin = -L+min(F(t));
    ymax = max(x)+lm;
    ylim([ymin, ymax])
    sgtitle([titleName,'Time = ',num2str(t(i),'%.2f'),' s'],...
            'FontWeight','bold','FontSize',fSize+4)
    

    subplot(1,2,2) % Solution on right
    plot(t,x,LineWidth=lWidth)
    hold on
    % Plot marker of solution in time
    plot(t(i),x(i),'ks',MarkerSize=mSize,MarkerFaceColor='k')
    % plot equilibrium line y = 0
    plot([t(1),t(end)],[0 0],'k--')

    xlabel('Time [s]','FontSize',fSize)
    ylabel('Displacement from Equilibrium [m]','FontSize',fSize)
    ylim([ymin, ymax])
    set(gca, 'FontSize', fSize);

    if strcmpi(saveVideo,'true')
        % Write the current frame to the video file
        set(gcf, 'Renderer', 'painters'); % or 'painters'
        writeVideo(videoObj, getframe(gcf));
    else
        if i == 1
            pause(2); % hold IC for visualization
        else
            pause(0.1);
        end
    end

    
end

if strcmpi(saveVideo,'true')
    % Close the video file
    close(videoObj);
end

% ------------------------------------------------------------------------
% Functions that set parameters for various scenarios
% ------------------------------------------------------------------------
function [m,b,k,F0,gamma,P0,V0,dt,tf] = undamped()
    m = 1; % mass [kg]
    b = 0; % damping constant [kg/s]
    k = 1; % spring constant [kg/s^2]
    
    % External forcing function F(t) = F0 * sin(gamma*t)
    F0 = 0;
    gamma = 0;
    
    % Initial conditions
    P0 = 2;
    V0 = 0;

    % Time parameters
    dt = 0.25; % time step for ODE solver
    tf = 60; % Final time for simulation
end

function [m,b,k,F0,gamma,P0,V0,dt,tf] = underdamped()
    [m,b,k,F0,gamma,P0,V0,dt,tf] = undamped();

    b = k/5; % damping constant [kg/s]

    % Time parameters
    dt = 0.25; % time step for ODE solver
    tf = 45; % Final time for simulation
end

function [m,b,k,F0,gamma,P0,V0,dt,tf] = overdamped()
    [m,b,k,F0,gamma,P0,V0,dt,tf] = undamped();

    b = 4*k; % damping constant [kg/s]

    % Time parameters
    dt = 0.25; % time step for ODE solver
    tf = 20; % Final time for simulation
end

function [m,b,k,F0,gamma,P0,V0,dt,tf] = critical()
    [m,b,k,F0,gamma,P0,V0,dt,tf] = undamped();

    b = 2*k; % damping constant [kg/s]

    % Time parameters
    dt = 0.25; % time step for ODE solver
    tf = 20; % Final time for simulation
end

function [m,b,k,F0,gamma,P0,V0,dt,tf] = forced()
    [m,b,k,F0,gamma,P0,V0,dt,tf] = underdamped();
    
    % External forcing function F(t) = F0 * sin(gamma*t)
    F0 = 1;
    gamma = 3/4;

    % Initial conditions
    P0 = 0;
    V0 = 0;
end

function [m,b,k,F0,gamma,P0,V0,dt,tf] = beats()
    [m,b,k,F0,gamma,P0,V0,dt,tf] = forced();
    
    % No damping
    b = 0;

    % Natural frequency 
    k = 3;
    omega = sqrt(k/m);

    % External forcing function F(t) = F0 * sin(gamma*t)
    F0 = 1;
    gamma = 0.875*omega;

    % Initial conditions
    P0 = 0;
    V0 = 0;

    % Time parameters
    dt = 0.25;
    tf = 60;
end

function [m,b,k,F0,gamma,P0,V0,dt,tf] = resonance()
    [m,b,k,F0,gamma,P0,V0,dt,tf] = forced();
    
    % No damping
    b = 0;

    % Natural frequency 
    k = 2;
    omega = sqrt(k/m);

    % External forcing function F(t) = F0 * sin(gamma*t)
    F0 = 0.5;
    gamma = omega;

    % Initial conditions
    P0 = 0;
    V0 = 0;

    % Time parameters
    dt = 0.25;
    tf = 60;
end
