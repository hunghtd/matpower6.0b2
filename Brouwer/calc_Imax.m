function Imax = calc_Imax(mpc, Yf, Vpfbase, rate)
%UNTITLED7 Summary of this function goes here
%   Detailed explanation goes here
% This function calculates the maximum currents or the thermal limits which
% equals to ``rate'' times the nominal currents
define_constants

I = Yf*Vpfbase;
Imax = abs(I)*rate;

on = find(mpc.branch(:, BR_STATUS) > 0);  %% which lines are on?

if length(on) == length(Imax)
    mpc.branch(on, end+1) = Imax; %% append Imax values to the branch data
else
    disp('Something is wrong with the lines')
end

end

