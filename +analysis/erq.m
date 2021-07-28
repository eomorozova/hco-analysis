function ERQ = erq(V, Vth)

% Inputs
% V     membrane potential
% Vth   synaptic threshold

% Output
% ERQ   Escape to Release Quotient (ERQ)

ERQ = (mean(V)-Vth)/mean(V);