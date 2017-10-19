function resR = constrainedCPF(testcase, busp, nump, rate, dV, stepsize, adaptive, limitQQ, limitIQ, Angle)

define_constants;
mpc = loadcase(testcase);

for k = 1:numel(busp)
    busp(k) = find(mpc.bus(:,1)==busp(k));
end

%renumbering the bus index
mpc = ext2int(mpc);
busp = mpc.order.bus.e2i(busp);%index of new buses after renumbering

%if the chosen buses are not pq buses, reselect those
pq = find((mpc.bus(:, BUS_TYPE) == 1) & (mpc.bus(:, PD) >0));
if sum(ismember(busp, pq)) ~= 2
    busp = pq(1:2);
end

SBASE = mpc.baseMVA;
genlist = mpc.gen(:, GEN_BUS);

[Y, Yf, ~] = makeYbus(mpc);

% base solution
pfbase = runpf(mpc);
Vpfbase = pfbase.bus(:,VM).*(exp(1i*pfbase.bus(:,VA)*pi/180));
s0 = Vpfbase.*conj(Y*Vpfbase);
V0 = abs(Vpfbase);

% calculate branch currents & Imax
Imax = calc_Imax(mpc, Yf, Vpfbase, rate);

% add qmax, qmin
[Qmax, Qmin] = calc_Qlim(mpc, s0, SBASE, genlist);

%% CPF
targetcase = mpc;
basecase = mpc;

mpopt = mpoption('cpf.stop_at', 'NOSE', 'cpf.step', stepsize, 'cpf.adapt_step', adaptive);

resR = zeros(nump,2); %save the loadability limits

%define # of quadrants to plot
if nargin < 10
    phi = linspace(0, pi/2, nump);
elseif nargin == 10
    phi = linspace(0, Angle, nump);
else
    error('requires at most 10 optional inputs')
end

% now run the CPF
for k = 1:nump
        targetcase.bus(busp, PD) = (1+ [cos(phi(k)); sin(phi(k))]) .* basecase.bus(busp, PD);
        targetcase.bus(busp, QD) = (1+ [cos(phi(k)); sin(phi(k))]) .* basecase.bus(busp, QD);

        cpf = runcpf(basecase, targetcase, mpopt);
        
       % find the loadability level satisfying thermal limits, voltage limits, and reactive power limits
       lam_max = calc_lamax(cpf, Yf, Imax, V0, dV, Y, genlist, Qmax, Qmin, limitQQ, limitIQ);
       
       resR(k, :) = (1+lam_max*[cos(phi(k)); sin(phi(k))]) .* basecase.bus(busp, PD)/SBASE; 
end

end
