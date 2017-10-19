function lammax = lamax_Vmax(cpf, step_max_V, V0, dV)
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here
if step_max_V > 0
    [checknextstep, vbus] = maxV_feas_Q( cpf.cpf.V(:,step_max_V+1), V0, dV );

    if checknextstep
       error('the critial step is wrong. Check bisection iteration')
    else
        Vv_flr = abs(cpf.cpf.V(vbus,step_max_V));
        Vv_ceil = abs(cpf.cpf.V(vbus,step_max_V+1));

        if Vv_flr > Vv_ceil %increase power -> voltage decrease, lower voltage violation
           Vmin = V0(vbus)-dV;
           fac = (Vv_flr - Vmin) / (Vv_flr - Vv_ceil);
        else
             Vmax = V0(vbus)+dV;
             fac = (Vmax - Vv_flr) / (Vv_ceil - Vv_flr);
        end 
    end

    lammax = fac * (cpf.cpf.lam(step_max_V+1) - cpf.cpf.lam(step_max_V)) + cpf.cpf.lam(step_max_V);  
else
    lammax = 0;
end

end

