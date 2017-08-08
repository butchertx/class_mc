function [ S1, S2 ] = split_action( actions, alpha, gamma, sx, Nt, Nx )
%SPLIT_ACTION split an array of actions into their contributions
%   action = S1*gamma + S2*alpha
%   S1 = -Nt*(1 - 2*sx)
S1 = -Nx*Nt*(ones(size(sx)) - 2.0*sx);
S2 = 1.0/alpha * (actions - gamma*S1);


end

