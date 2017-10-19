function  lammax = find_lamax(violation, cpf, Yf, Imax, V0, dV, Y, genlist, Qmax, Qmin, step_max_I, step_max_V, step_max_Q)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

switch violation
    case 1 %current limits are violated
             lammax = lamax_Imax(cpf, Yf, Imax, step_max_I);
    case 2 %voltage limits are violated
            lammax = lamax_Vmax(cpf, step_max_V, V0, dV);
    otherwise %reactive power limits are violated
           lammax = lamax_Qmax(cpf, step_max_Q, Y, genlist, Qmax, Qmin);
end

end

