clear;clc;
close all;

%select data type
% type = 'Loose';
% type = 'Dense';
% type = 'Steel';
type = 'Sub_Steel';

n = 5; % kernel size
load([type '_' num2str(n) '_Dimensional_porosity.mat']);
%% Prepare variables
lenz = length(convex_por_z);
lenx = length(convex_por_x);
leny = length(convex_por_y);
z = [1:lenz];
x = [1:lenx];
y = [1:leny];

%% Plot directional porosity

% x-axis direction
h=figure
hold on

% scatter(x,porosity_x,12,'k','.')
plot(x,convex_por_x,'LineWidth' ,2)

set(gca, 'FontSize',16)
legend('X-direction','Location','best','FontSize',16)
box on

ax = gca;
ax.YLabel.String = 'Porosity';
ax.XLabel.String = 'Position (\it\mu\itm)';
ax.LineWidth=1
set(h,'Units','Inches');
xlim([0 lenx(end)])
ylim([0 1])



% y-axis direction
h=figure
hold on

% scatter(x,porosity_x,12,'k','.')
plot(y,convex_por_y,'LineWidth' ,2)

set(gca, 'FontSize',16)
legend('Y-direction','Location','best','FontSize',16)
box on

ax = gca;
ax.YLabel.String = 'Porosity';
ax.XLabel.String = 'Position (\it\mu\itm)';
ax.LineWidth=2
set(h,'Units','Inches');
xlim([0 leny(end)])
ylim([0 1])


% z-axis direction
h=figure
hold on

% scatter(x,porosity_x,12,'k','.')
plot(z,convex_por_z,'LineWidth' ,2)

set(gca, 'FontSize',16)
legend('Z-direction','Location','best','FontSize',16)
box on

ax = gca;
ax.YLabel.String = 'Porosity';
ax.XLabel.String = 'Position (\it\mu\itm)';
ax.LineWidth=2
set(h,'Units','Inches');
xlim([0 lenz(end)])
ylim([0 1])

