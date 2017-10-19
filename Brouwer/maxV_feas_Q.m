function [feas_Q, vbus] = maxV_feas_Q(V, V0, dV )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

Vm = abs(V);
vbus = find(abs(Vm - V0) > dV, 1);%violated buses
feas_Q = isempty(vbus);

end

