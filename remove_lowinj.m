function [deltan_rev,tau_rev] = remove_lowinj(deltan,tau,cutoff)

indices = find(deltan<cutoff);
deltan_rev = deltan;
tau_rev = tau;
deltan_rev(indices) = [];
tau_rev(indices) = [];