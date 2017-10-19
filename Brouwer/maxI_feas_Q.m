function [feas_Q, vline] = maxI_feas_Q(Yf, V, Imax )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

Im = abs(Yf*V);
vline = find(Im > Imax, 1);%violated lines
feas_Q = isempty(vline);

end

