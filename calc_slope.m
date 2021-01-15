function [S] = calc_slope(Z,L)
% calc_slope.m calculates local channel slope (S) from vectors of
% elevation (Z) and distance (L)
%
% example
% S = calc_slope(Z,L);
%
% Author: Sean F. Gallen
% Date modified: 02/17/2020

% make a dummy slope vector based on the elevation vector
S = nan(size(Z));

% calculate the slope for all nodes except for the outlet (final) node
S(1:end-1) = abs(Z(1:end-1)-Z(2:end))./abs(L(2:end) - L(1:end-1));

% Calculate slope at the outlet
S(end) = S(end-1);

end

