function [KE, PE] = getEnergy(positions, velocities, masses, gravitationalConstant)
%   positions:   N x 3 matrix representing the positions of particles [x, y, z]
%   velocities:   N x 3 matrix representing the velocities of particles
%   masses:      N x 1 vector representing the masses of particles
%   gravitationalConstant:   Newton's Gravitational constant
%   kineticEnergy:   The kinetic energy of the system
%   potentialEnergy:   The potential energy of the system

% Kinetic Energy Calculation:
KE = 0.5 * sum(sum(masses .* velocities.^2));

% Potential Energy Calculation:
% Extracting individual components of positions

x = positions(:, 1);
y = positions(:, 2);
z = positions(:, 3);

% Creating matrices for pairwise particle separations: r_j - r_i
dx = x' - x;
dy = y' - y;
dz = z' - z;

% Creating a matrix for particle pairwise distances
r = sqrt(dx.^2 + dy.^2 + dz.^2);

% Summing over the upper triangle to count each interaction only once

PE = gravitationalConstant * sum(sum(triu(-(masses * masses')./r, 1)));
end
