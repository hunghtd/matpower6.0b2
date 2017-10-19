function resR = CPFfun(testcase, nump)
define_constants;
mpc = loadcase(testcase);
mpc = ext2int(mpc);
SBASE = mpc.baseMVA;

slack = find(mpc.bus(:, BUS_TYPE) == 3);
pq = find(mpc.bus(:, BUS_TYPE) == 1);
% pv = find(mpc.bus(:, BUS_TYPE) == 2);
% nslack = size(slack, 1);
% npq = size(pq, 1);
% npv = size(pv, 1);
% pvpq = [pv; pq];

Y = makeYbus(mpc);
nb = size(Y,1);


% base solution
pfbase = runpf(mpc);
Vpfbase = pfbase.bus(:,VM).*(exp(1i*pfbase.bus(:,VA)*pi/180));
Vbase =  Vpfbase(1:end ~= slack);
s0 = Vpfbase.*conj(Y*Vpfbase);

% target case
targetcase = mpc;
basecase = mpc;

% CPF
PM = 100;
mpopt = mpoption('cpf.stop_at', 'NOSE', 'cpf.step', 0.01);

lambdamax = [];

resR = zeros(nump,2);
phi = linspace(0, pi/2, nump);

for k = 1:nump
  unitvec = [];
  unitvec(1) =  - cos(phi(k));
  unitvec(2) =  - sin(phi(k));
    if nb==5
      targetcase.bus(pq(1), PD) = basecase.bus(pq(1), PD) -  PM*SBASE*unitvec(1)';
      targetcase.bus(pq(1), QD) = basecase.bus(pq(1), QD) -  PM*SBASE*unitvec(2)';
    else
        targetcase.bus(pq(1:2), PD) = basecase.bus(pq(1:2), PD) -  PM*SBASE*unitvec(1:2)';
    end

    %if numel(pq) > 2
      %targetcase.bus(pq(3:end), PD) = 2*basecase.bus(pq(3:end), PD);
    %end

    %targetcase.bus(pq, QD) = 2*basecase.bus(pq, QD);

    %targetcase.gen(:, PG) = 1.001*basecase.gen(:, PG);
    if nb ==118
         [Rcpf, suc, steps, lam] = runRcpf_v2(basecase, targetcase, 1, 10000);
         temp = Rcpf.bus(:, PD)/SBASE;
         resR(k,:) = temp([pq(1:2)]);
    else
            cpf = runcpf(basecase, targetcase, mpopt);
            lammax =  PM*cpf.cpf.max_lam;
            if nb==5 
                resR(k, :) =  [(-real(s0(pq(1)))'- lammax*unitvec(1))  (-imag(s0(pq(1)))' - lammax*unitvec(2)) ];
            else
                resR(k, :) =  -real(s0(pq(1:2)))- lammax*unitvec(1:2)';
            end
    end
   
end
end
