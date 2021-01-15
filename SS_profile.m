function [S, Lo, Z] = SS_profile(U, K, m, n, A, L)
% SS_profile.m takes inputs of uplift rate (U), erodibility (K), the
% drainage area exponent (m), the slope exponent (n), and vectors of
% drainage area (A) and channel distance (L) both ordered from the channel 
% head to the outlet.
% 
% SS_profile uses Flint's Law to predict the steady-state elevation of a 
% river profile.
%
% There are three outputs:
% S - the steady state slope from Flint's Law
% Lo - profile distance with respect to the outlet and ordered from the
% outlet to the channel head.
% Z - the elevation of the river network ordered from channel head to
% outlet.
%
% Example:
% [S, Lo, Z] = SS_profile(U, K, m, n, A, L);
%
% Author: Sean F. Gallen
% Date modified: 02/17/2020

% calculate the steady-state slope from Flint's Law
ks = (U/K).^(1/n);

S = ks.*A.^(-m/n);

% Integrate to get elevation
% (1) set up variables
Lo = max(L) - L;

% (2) flip vectors
Lo = fliplr(Lo);
So = fliplr(S);

% (3) integrate
Z = zeros(size(Lo));

for i = 2:length(Z)
    dz = ((So(i) + So(i-1))/2)*(Lo(i)-Lo(i-1));
    Z(i) = Z(i-1) + dz;
end

% (4) flip Z vector
Z = fliplr(Z);

