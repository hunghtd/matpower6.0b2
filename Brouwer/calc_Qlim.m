function [Qmax, Qmin] = calc_Qlim(mpc, s0, SBASE, genlist)
%UNTITLED8 Summary of this function goes here
%   Detailed explanation goes here
define_constants

Qmax = mpc.gen(:, QMAX)/SBASE;
Qmin = mpc.gen(:, QMIN)/SBASE;

q = imag(s0);
qmax_vgen = find(q(genlist) > Qmax);
qmin_vgen = find(q(genlist) < Qmin);
Qmax(qmax_vgen ) = 1e+7; %if reactive power maximum limits are violated, raise the limits
Qmin(qmin_vgen ) = -1e+7; %if reactive power minimum limits are violated, lower the limits
end

