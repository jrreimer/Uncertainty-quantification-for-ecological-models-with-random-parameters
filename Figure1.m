% recreates Figures 1a, b
% requires `NP_det.m` in the working directory
% saves `Figure1a.pdf` and `Figure1b.pdf` in working directory

clear 
close all
global a b c d e rho

%_______________________________Define colors______________________________
algaecolordet = 1/255*[118,176,65]; % color for deterministic P results (green)
nutrientcolordet = 1/255*[255,201,20]; % color for deterministic N results (yellow)
algaecolorRDE = 1/255*[125,91,166]; % RDE algae color (purple); alternatively pink: [252, 100, 113]
nutrientcolorRDE = 1/255*[23,190,187]; % RDE nutrient color (blue)

%_________________________________Parameters_______________________________

a=0; % alpha
b=1; % beta
c=0.8; % gamma
d=0.1; % delta
e=0; % eta
rho = 0.4; % rho

% Initial conditions 
xhat = 0.22; % nutrient initial conditions
yhat = 0.001; % algae initial conditions
tmax = 200; % max time
tspan = 0:tmax; % time span for solver

%__________________________________________________________________________
% Figure 1a - Deterministic solution

[Tdet,Zdet] = ode45(@NP_det,tspan,[xhat;yhat]); %=nut,y=algae

fig=figure;
set(fig,'defaultAxesColorOrder',[0 0 0; 0 0 0]);
hold on

yyaxis left % nutrients on left axis
plot(Tdet,Zdet(:,1),'color',nutrientcolordet,'linewidth',3) 
ylim([0, 0.265])
ylabel('nutrients')
xlabel('time (days)')

yyaxis right % algae on right axis
plot(Tdet,Zdet(:,2),'color',algaecolordet,'linewidth',3)
ylim([0, 0.045]);
ylabel('algae')

set(gca,'fontsize',18) 
legend('N, nutrients','P, algae','Location','northeast')
legend boxoff    

% add labels to plot
NW = [min(xlim) max(ylim)]+[diff(xlim) -diff(ylim)]*0.02;
text(NW(1), NW(2), '(a)', 'VerticalAlignment','top', 'HorizontalAlignment','left','fontsize',18) 
text(Tdet(Zdet(:,2)==max(Zdet(:,2)))-3, max(Zdet(:,2))*1.01, '\bullet','FontSize',18,'color',algaecolordet)
text(Tdet(Zdet(:,2)==max(Zdet(:,2)))-3, max(Zdet(:,2))*1.05, 'P_{max}','FontSize',14)
yyaxis left
text(max(xlim)*.97, Zdet(end,1)*1.04, '\bullet','FontSize',18,'color',nutrientcolordet)
text(max(xlim)*.91, Zdet(end,1)*1.3, 'N*','FontSize',18)
hold off

set(gcf, 'units','inches','Position',  [2, 2, 6, 4])
set(gcf, 'units','inches','PaperSize', [6 4]); 
saveas(gcf,'Figure1a.pdf')

%__________________________________________________________________________
% Figure 1b - Monte Carlo Simulations 

% analytic expression for E[P_max]
EPmax = (-d/(2*rho))*log((b+rho)/(b-rho))*(1-log(d/(c*xhat*sqrt(b^2-rho^2)))) + yhat+c*xhat;

Nsim = 100; % number of MC simulations
rng(123); % set seed for random samples
uvals = -1 + 2*rand(Nsim,1); % generate vector of random values between [-1,1]
bvals = b + rho*uvals; % scale uvals to become values of beta
sample_ts = [0:tmax]; 
MCsims = zeros(length(sample_ts),Nsim); % each column saves one simulation
MCNvals = zeros(length(sample_ts),Nsim);
for j=1:Nsim
    b = bvals(j);
    [Tvals,Zvals] = ode45(@NP_det,sample_ts,[xhat;yhat]); 
    MCsims(:,j) = Zvals(:,2); % algae population is second column of Zvals
    MCNvals(:,j) = Zvals(:,1); % nutrient population is first column of Zvals
end
MCmean = mean(MCsims,2);
MCNmean = mean(MCNvals,2);

fig=figure;
set(fig,'defaultAxesColorOrder',[0 0 0; 0 0 0]);
hold on

yyaxis right % algae on right axis
plot(Tvals,MCsims,'-','color',[algaecolorRDE,0.05],'Linewidth',1)
algaeaxis = ylim;
ylim([0 0.045]);
hm = plot(Tvals,MCmean,'-','color',algaecolorRDE,'Linewidth',3);
xlim([0,200])
xlabel('time (days)')
ylabel('algae')

yyaxis left % nutrients on left axis
plot(Tvals,MCNvals,'-','color',[nutrientcolorRDE,0.05],'linewidth',1)
hm2 = plot(Tvals,MCNmean,'-','color',nutrientcolorRDE,'Linewidth',3);
ylim([0, 0.265])
ylabel('nutrients')

set(gca,'fontsize',18) 
legend([hm2 hm],'E[N]','E[P]','Interpreter', 'none','Location','best')

NW = [min(xlim) max(ylim)]+[diff(xlim) -diff(ylim)]*0.02;
text(NW(1), NW(2), '(b)', 'VerticalAlignment','top', 'HorizontalAlignment','left','fontsize',18)
yyaxis right
text(Tvals(MCmean==max(MCmean))*0.95, max(MCmean)*1.02, '\bullet','FontSize',18,'color',algaecolorRDE)
text(Tvals(MCmean==max(MCmean))*.7, max(MCmean)*0.9, 'max(E[P])','Interpreter','tex','FontSize',14)
text(Tvals(MCmean==max(MCmean))*0.95, EPmax*1.01, '\bullet','FontSize',18,'color',algaecolorRDE);
text(Tvals(MCmean==max(MCmean))*0.95, EPmax*1.08, 'E[P_{max}]','Interpreter','tex','FontSize',14)
yyaxis left
text(max(xlim)*.97, MCNmean(end,1)*1.04, '\bullet','FontSize',18,'color',nutrientcolorRDE)
text(max(xlim)*.85, MCNmean(end,1)*1.3, 'E[N*]','Interpreter','tex','FontSize',14)
% can play with labels using Matlab plot GUI

legend boxoff                   
hold off
set(gcf, 'units','inches','Position',  [2, 2, 6, 4])
set(gcf, 'units','inches','PaperSize', [6 4]); 

saveas(gcf,'Figure1b.pdf')
