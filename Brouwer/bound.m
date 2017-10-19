close all
addpath(genpath('/Users/hunghtd/Dropbox (Personal)/Journal/Feasibility region/matpower6.0b2'));
xm = pi/2;
x = linspace(-xm, xm);
ycos = sin(x);
plot(x,ycos,'r','linewidth', 2.5);
hold on

y0 = zeros(100,1);
plot(x,y0,'k--','linewidth', 2.5)

hold on

A = [1 0; 0 1; -1 0; 0 -1];
b = [xm; sin(xm); xm; sin(xm)];

V = con2vert(A, b);
x = V(:,1);
y = V(:,2);
h = convhull(x,y);
plot(x(h), y(h), 'linewidth', 3)

xticks([-xm 0 xm])
xticklabels({'-\pi/2','0','\pi/2'})
yticks([-sin(xm) 0 sin(xm)])
yticklabels({'-f^-','0', 'f^+'})
set(gca,'Fontname','Times New Roman','fontsize',30)
legend( 'sin(x)')

%%

plot(resR(:,1), resR(:,2), 'r','linewidth', 2.5)
scatter(basep(1), basep(2), [],'k', '*', 'LineWidth',4)
%axis([basep(1) inf basep(2) inf])
%title(['Solvability boundaries for IEEE ',num2str(nb),' bus system'])
xlabel(['Active power at load bus ', num2str(busp(1)), ' (P_{', num2str(busp(1)), '})'])
ylabel(['Active power at load bus ', num2str(busp(2)), ' (P_{', num2str(busp(2)), '})'])
set(gca,'Fontname','Times New Roman','fontsize',24)
legend( 'Inner approximation', 'Feasibility boundary', 'Base operating point')