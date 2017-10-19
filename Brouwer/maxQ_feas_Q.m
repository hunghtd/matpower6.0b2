function [feas_Q, vgen, maxmin] = maxQ_feas_Q( Y, V, genlist, Qmax , Qmin)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
% this function is to check the feasibility of reactive power generation
% feas_Q = 1: feasible, 0: infeasible
% vgen: vilated generator
% maxmin = 1: maximum Qgen is violated; maxmin = 2, minimum Qgen is violated

q = imag(V.*conj(Y*V));
vgen_max = find(q(genlist) > Qmax, 1);
vgen_min = find(q(genlist) < Qmin, 1);

feas_Q = 1;

if ~isempty(vgen_max)
    feas_Q = 0;
    maxmin = 1; %maximum Qlim is violated
    vgen = vgen_max;
else    
   if ~isempty(vgen_min)
       feas_Q = 0;
       maxmin = 2; %minimum Qlim is violated
       vgen = vgen_min;
   end
end

end

