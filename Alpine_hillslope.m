%% Alpine hillslopes
% Computational modelling class GEOL5700-004, Project 3
% Kelly Kochanski, 08-Feb-2016


%%-------------------------------------------------------------------------
% This code describes the 1D evolution of an alpine hillslope between two 
%  glacial valleys. 

% Soil is produced by WEATHERING, which is a function of depth to bedrock.

% The HILLSLOPE evolves by frost creep. Soil expands
%  perpendicular to slope each winter, then falls vertically. Expansion 
%  height is limited by soil thickness and freezing depth, after Anderson
%  2002. 

% GLACIAL VALLEYS incise and form steep walls (modelled as vertical) during
%  periods of glacial advance. These walls then weather during interglacial
%  periods.

% AIM: The Wind River Mountains of WY form low, parabolic hillslopes which
%  are completely detached from the glacial valleys, far below. However, in
%  the White Mountains in NH, glacial valleys (which were covered by the
%  same Laurentide ice sheet) are smooth.
%  I will investigate which conditions (weathering rate, frost creep rate,
%  and duration of glacial periods) form disconnected hillslopes as in WY,
%  and which form smooth hillslope/valley systems as in NH.
%%-------------------------------------------------------------------------
% HILLSLOPE VARIABLES AND FUNCTIONS
% At every point on the hillslope, we track:
%   x       - horizontal position, measured from hillcrest, m, vector in i
%   zB(x,t) - elevation of bedrock,                         m, vector in i
%   zS(x,t) - elevation of soil/regolith surface,           m, vector in i
%   h (x,t) - soil/regolith thickness (zS-zB),              m, vector in i
%   S (x,t) - dzS/dx, slope                                 m/m, vector, i
%
%   w(h)    - weathering rate                               m/yr
%   q(S)    - regolith movement rate                        m/yr

% PARAMETERS
rho_r = 3300; % density of rock                             kg/m^3
rho_s = 1500; % density of soil/regolith                    kg/m^3

%% Weathering rate
% This function is a humped distribution, where the weathering rate
% increases with soil depth h while h<h0, then decreases exponentially.
h0    = 1; %                                              m, NON-ZERO  
max_w = 20*10^-6; % maximum weathering rate,                m/yr
w = @(h) 10^-5; %                                          m/yr
% According to Anderson 2002, this formula should be 
% w(h) = h.*exp(-h./h0).*(max_w.*exp(1)./h0);
% But this keeps blowing up. I couldn't figure it out and gave up.

%% Frost creep - effective diffusivity
% Following Anderson 2002, Fig. 9
max_diffusivity = 0.02; % m^2/yr
k_eff = @(h) max_diffusivity.*(1-exp(-h.*2));

%% Glacial incision
% Glacial incision is either "on" during glacial periods,
% or "off" during interglacials.
%incision in valley during glacial,                         m/yr
% This increases away from the edge to make a more U-like valley
glacial_incision = @(x,N) (-x+N*dx)/10^4; 
glacial_duration = 1000; %duration of glacial period,       yr
glacial_starts   = [-11000, -23000, -35000]; %last glacial started...yrs ago

%% Geometry of hillslope
% We will model half the hill (and assume it is symmetric)
% The hillcrest is at x=0; the valley is at x>=x_hill.
x_hill   = 400; % Horizontal range of hillslope             m
x_valley = 100; % Horizontal range of glacial valley        m

% Horizontal step size for numerical grid
dx        = 5; %                                            m
xvector   = 0:dx:(x_hill+x_valley); %                       m, vector
is_valley = xvector > x_hill; %                             true/false
N         = length(xvector);

% Initial conditions 
%(gentle parabola - like a steady state solution - with thin soil and 
% a little noise)
zB0 = - xvector.^2/10^3;%                                  m, vector
zS0 = zB0 + 10 + 1*rand(1,N); %                            m, vector
zB0 = zB0 + 10*rand(1,N);

%% Start running the problem
% Time 
dt       = 10; % time step, initial guess                   yrs
start_time = 5*10^4; % years ago                            yrs
end_time = 5*10^4; % years in future
dt       = dx.^2/4/max_diffusivity; %highest stable step    yrs

% Measure time in years before present
t = - start_time; 
zB = zB0; zS = zS0;
plot(xvector, zB0, 'k'); hold on
plot(xvector, zS0, 'g');
set(gca, 'fontsize', 14)
legend('Bedrock', 'Soil')
ylabel('Elevation (m)')
xlabel('Distance from hillcrest (m)')
title('Initial topography')
disp('Showing initial bedrock and soil topography. Press any key to accept.')
pause()
hold off

% Start the iterations
while t < end_time;
    % Are we in a glacial period?
    is_glaciated = false;
    for i=1:length(glacial_starts)
        if t > glacial_starts(i) && t < glacial_starts(i) + glacial_duration
            is_glaciated = true;
            disp(sprintf('Ice age! t = %d kyr', round(t/1000)))
        end
    end
    
    % Determine soil thickness
    
    h = zS - zB;
    % Weather the bedrock everywhere
    weathered_rock = w(h).*dt;
    % Update rock/surface elevations and soil thickness
    % Do not let soil thickness fall below zero
    zB = zB - weathered_rock;
    zS = max(zB, zS + weathered_rock*(rho_r/rho_s - 1));
    h  = zS - zB;
    
    % Calculate diffusivities
    diffusivities = k_eff(h);
    % Calculate soil fluxes and move soil along hillslope
    zS(1)            = zS(2);
    dzdx(1:N-1)      = diff(zS)/dx;
    % x = L boundary: blocked (assuming non-glacial stream is tiny)
    dzdx(N)          = 0.0;
    soil_flux        = - diffusivities.*dzdx;
    dqdx             = diff(soil_flux)/dx;
    zS(2:end)        = zS(2:end) - dqdx;
    % x = 0 boundary: zero slope
    zS(1)            = zS(2);   
    
    % If glaciated, in valley glacier incises and takes away lots of soil
    zS = zS - is_glaciated.*is_valley.*dt.*h;
    zB = zB - is_glaciated.*is_valley.*glacial_incision(xvector,N).*dt;
    
    % Make sure soil never has negative thickness
    zS               = max(zS, zB);
    
    % Plot things every n_plot time steps
    n_plot = 10;
    if round(mod(t/dt, n_plot)) == 1
        hold on
        set(gca, 'fontsize', 14)
        plot(xvector, zB, 'b')
        plot(xvector, zS, 'g')
        ylabel('Elevation (m)')
        xlabel('Distance from hillcrest (m)')
        title(sprintf('Time: %d kyr', round(t/1000)))
        pause(0.1)
    end
    
    t = t + dt;
        
end

hold on
set(gca, 'fontsize', 14)
plot(xvector, zB, 'k', 'linewidth', 2)
plot(xvector, zS, 'k', 'linewidth',2)
title(sprintf('Time: %d kyr', round(t/1000)))
disp('-------------------')
disp(sprintf('Simulation ran from %d kyr before present to %d kyr after present.',...
    round(start_time/1000), round(end_time/1000)))
disp(' The hillslope experienced three periods of glaciation.')
disp(' The unglaciated hillslope, which was originally rough, is now')
disp('covered with a smooth, parabolic soil surface due to frost creep.')
disp(' The glaciated area is deeply eroded. Soil has filled in some corners,')
disp('but a steep cliff remains.')
disp('--------------------')
disp(' I hypothesized that some mountain ranges have cliff-edged glacial')
disp('valleys, while some have smoother valleys, because areas with high')
disp('soil diffusivity quickly filled their valleys with smooth soil.')
disp(' However, I was unable to keep the upper soil attached to the valley for any')
disp('physically reasonable combination of parameters. This is therefore a')
disp('reasonable model for the Wind River Range, but not for the smoother')
disp('White Mountains.')
disp(' My guess? Frost creep is very slow. In the White Mountains, much of')
disp('the erosion is done by landslides on steep sections. Landslides move')
disp('regolith far more efficiently than frost creep does. Perhaps the Wind')
disp('River Range has more resilient rocks and fewer landslides, or was')
disp('more recently glaciated and has not had time to erode')