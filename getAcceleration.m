function [acceleration] = getAcceleration (positions, masses, gravitationalConstant, softeningLength)
%   positions:   N x 3 matrix representing the positions of particles [x, y, z]
%   masses:      N x 1 vector representing the masses of particles
%   gravitationalConstant:   Newton's Gravitational constant
%   softeningLength:   Softening length to avoid divisions by zero
%   acceleration:   N x 3 matrix representing the accelerations along each axis

%   Extracting individual components of positions

x = positions(:, 1);
y = positions(:, 2);
z = positions(:, 3);

% Creating matrices for pairwise particle separations: r_j - r_i

dx = x' - x;
dy = y' - y;
dz = z' - z;

% Creating a matrix for 1/r^3 for all pairwise particle separations
inv_r3 = (dx.^2 + dy.^2 + dz.^2 + softeningLength.^2).^(-3/2);

% Fixing diagonal values (representing self-interaction) which result in 1/0=Infinity. Setting them to 0.

inv_r3(isinf(inv_r3)) = 0;

% Calculating acceleration components along each axis
ax = gravitationalConstant * (dx .* inv_r3) * masses;
ay = gravitationalConstant * (dy .* inv_r3) * masses;
az = gravitationalConstant * (dz .* inv_r3) * masses;

% Combining the acceleration components into a single matrix
acceleration = [ax ay az];
end
