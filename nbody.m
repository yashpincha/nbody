% N-body Simulation Parameters
numParticles = 10;         % Number of particles
currentTime = 0;           % Current time of the simulation
simulationEndTime = 10;     % Time at which the simulation ends
timeStep = 0.01;            % Timestep
softeningLength = 0.1;      % Softening length
gravitationalConstant = 1;  % Newton's Gravitational Constant
realTimePlotting = 1;       % Switch on (1) for plotting as the simulation progresses

%{
The realTimePlotting variable is used to control whether the simulation is visually plotted in real-time as it progresses. 
When realTimePlotting is set to 1 (true), the simulation will display the positions of particles and the energy plot dynamically as the 
simulation evolves through each time step. This can be useful for visualizing the behavior of the simulation in real-time.
If realTimePlotting is set to 0 (false), the plots will only be displayed at the end of the simulation.
%}

% Generate Initial Conditions
rng(42);                   % Set the random number generator seed
totalMass = 20;            % Total mass of particles is 20
mass = totalMass * ones(numParticles, 1) / numParticles;
positions = randn(numParticles, 3);  % Randomly selected positions and velocities
velocities = randn(numParticles, 3);

% Convert to Center-of-Mass Frame
velocities = velocities - mean((mass * [1 1 1]) .* velocities) / mean(mass);

% Calculate initial gravitational accelerations
acceleration = getAcceleration(positions, mass, gravitationalConstant, softeningLength);

% Calculate initial energy of the system
[KE, PE] = getEnergy(positions, velocities, mass, gravitationalConstant);

% Number of timesteps
numTimesteps = ceil(simulationEndTime / timeStep);

%{
The ceil function in MATLAB is used to round each element of an array to the nearest integer greater than or equal to that element. 
In other words, it rounds numbers up to the next whole number. If a number is already an integer, ceil leaves it unchanged.
%}

% Save energies and particle orbits for plotting trails
savedPositions = zeros(numParticles, 3, numTimesteps + 1);
savedPositions(:, :, 1) = positions;
savedKineticEnergy = zeros(numTimesteps + 1, 1);
savedKineticEnergy(1) = KE;
savedPotentialEnergy = zeros(numTimesteps + 1, 1);
savedPotentialEnergy(1) = PE;
allTimes = (0:numTimesteps) * timeStep;

% Simulation Main Loop
figureHandle = figure('position', [0 0 600 800]);

for timestep = 1:numTimesteps
    
    % (1/2) Kick
    velocities = velocities + acceleration * timeStep / 2;
    
    % Drift
    positions = positions + velocities * timeStep;
    
    % Update accelerations
    acceleration = getAcceleration(positions, mass, gravitationalConstant, softeningLength);
    
    % (1/2) Kick
    velocities = velocities + acceleration * timeStep / 2;
    
    % Update time
    currentTime = currentTime + timeStep;
    
    % Get energy of the system
    [KE, PE] = getEnergy(positions, velocities, mass, gravitationalConstant);
    
    % Save energies, positions for plotting trail
    savedPositions(:, :, timestep + 1) = positions;
    savedKineticEnergy(timestep + 1) = KE;
    savedPotentialEnergy(timestep + 1) = PE;
    
    % Plot in real time
    if (realTimePlotting) || (timestep == numTimesteps)
        subplot(3, 1, 1:2)
        xx = savedPositions(:, 1, max(timestep - 50, 1):timestep);
        yy = savedPositions(:, 2, max(timestep - 50, 1):timestep);
        plot(xx(:), yy(:), '.', 'color', [0.8, 0.6, 1]);
        hold on
        plot(positions(:, 1), positions(:, 2), 'm.', 'markersize', 14);
        hold off
        axis square
        axis([-2 2 -2 2])
        
        subplot(3, 1, 3)
        plot(allTimes, savedKineticEnergy, 'r.')
        hold on
        plot(allTimes, savedPotentialEnergy, 'b.')
        plot(allTimes, savedKineticEnergy + savedPotentialEnergy, 'k.')
        hold off
        axis([0 simulationEndTime -300 300])
        
        drawnow
    end
end

% Add labels/legend
subplot(3, 1, 3)
xlabel('Time')
ylabel('Energy')
legend('Kinetic Energy (KE)', 'Potential Energy (PE)', 'Total Energy (Etot)', 'Location', 'Northeast');
