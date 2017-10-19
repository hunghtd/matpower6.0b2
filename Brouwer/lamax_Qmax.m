function lam_max = lamax_Qmax(cpf, step_max_Q, Y, genlist, Qmax, Qmin)
%UNTITLED6 Summary of this function goes here
%   Detailed explanation goes here
define_constants

if step_max_Q > 0
    [checknextstep, vgen, maxmin] = maxQ_feas_Q( Y, cpf.cpf.V(:,step_max_Q+1), genlist, Qmax , Qmin );
    vbus = genlist(vgen);

    if checknextstep
       error('the critial step is wrong. Check bisection iteration')
    else
        calc_qgen = @(Y,V) imag(V(vbus)*conj(Y(vbus, :)*V));

        Qv_flr = calc_qgen(Y, cpf.cpf.V(:,step_max_Q));
        Qv_ceil = calc_qgen(Y, cpf.cpf.V(:,step_max_Q+1));

        switch maxmin
            case 1 %maximum Qlim is violated, Qv_flr < Qv_ceil
                if Qv_flr < Qv_ceil
                     fac = (Qmax(vgen) - Qv_flr) / (Qv_ceil - Qv_flr);
                end
           case 2 %minimum Qlim is violated
               if Qv_flr > Qv_ceil
                    fac = (Qv_flr - Qmin(vgen)) / (Qv_flr - Qv_ceil);
               end
        end
    end

    lam_max = fac * (cpf.cpf.lam(step_max_Q+1) - cpf.cpf.lam(step_max_Q)) + cpf.cpf.lam(step_max_Q);     

else
    lam_max = 0;
end

end

