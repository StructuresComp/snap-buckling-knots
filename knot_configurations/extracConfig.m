clc; clear all; close all

data = importdata("../datafiles/simDER_c0.50.txt");

nv = 301;

Nsteps = length(data)/nv;

for i = 1000
    config = data((i-1) * nv + 1:i*nv, :);
end

plot3(config(:,1), config(:,2), config(:,3));
axis equal;