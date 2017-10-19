function  polytop_plot(resB, resR, basep, busp)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
addpath(genpath('/Users/hunghtd/Dropbox (Personal)/Journal/Feasibility region/matpower6.0b2'));
A = [1 0; 0 1; -1 0; 0 -1];
b = [resB(1, 1); resB(2, 1); -resB(1, 2); -resB(2, 2)];

V = con2vert(A, b);
x = V(:,1);
y = V(:,2);
h = convhull(x,y);

plot(x(h), y(h), 'linewidth', 3)

hold on
plot(resR(:,1), resR(:,2), 'r','linewidth', 2.5)
scatter(basep(1), basep(2), [],'k', '*', 'LineWidth',4)
%axis([basep(1) inf basep(2) inf])
%title(['Solvability boundaries for IEEE ',num2str(nb),' bus system'])
xlabel(['Active power at load bus ', num2str(busp(1)), ' (P_{', num2str(busp(1)), '})'])
ylabel(['Active power at load bus ', num2str(busp(2)), ' (P_{', num2str(busp(2)), '})'])
set(gca,'Fontname','Times New Roman','fontsize',24)
legend( 'Inner approximation', 'Feasibility boundary', 'Base operating point')

end

