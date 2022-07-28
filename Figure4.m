% recreates Figures 4a,b,c
% requires `NP_det.m`, `NP_unc.m`, and `scenario1.mat` in the working directory
% saves `Figure4a.pdf`, `Figure4b.pdf`, and `Figure4c.pdf` in working directory

%__________________________________________________________________________
clear
close all
global A B C a b c d e rho N

%_______________________________Define colors______________________________
algaecolor = 1/255*[125,91,166]; % algae color (purple)
nutrientcolor = 1/255*[23,190,187]; % nutrient color (cyan)
sdshading = 1/255*[94,101,114]; % standard deviation shading color (grey)

%_________________________________Parameters_______________________________

a=0; % alpha
b=1; % beta
c=0.8; % gamma
d=0.1; % delta
e=0; % eta
rho = 0.2; % rho

% Initial conditions 
xhat = 0.22; % nutrient initial conditions
yhat = 0.001; % algae initial conditions
tmax = 160; % max time
sample_ts = 0:tmax; % time span for solver

%__________________________________________________________________________
% % Figure 4a and c - Monte Carlo Simulations & nonintrusive PC

load('scenario1.mat') % load data created by `figure-4.py`

%__________________________________________________________________________
% Figure 4b - Intrusive Polynomial Chaos

%%% Define Legendre polynomials via recurrence relationship %%%
N = 5; % highest polynomial degree 
Np = N+1; % number of polynomials (n in paper notation)
M(Np,Np) = 0;
M(1,1) = 1;
M(2,2) = 1;
for j=3:Np
    M(1,j) = -(j-2)*M(1,j-2)/(j-1);
    for i=2:Np
        M(i,j)=((2*j-3)*M(i-1,j-1)-(j-2)*M(i,j-2))/(j-1);
    end
end

%%% Integrate powers over [-1,1] with the weight function 1/2 %%%
for i=1:2*Np
    n = 2*(i-1);
    Intx(n+1) = 1/(n+1);
end

%%% Construct A aand B matrices, and C vector %%%
for lp=1:Np
    D(lp) = conv(M(:,lp),M(:,lp))'*Intx(1:2*N+1)';
    for ip=1:Np
        for jp=1:Np
            V = conv(conv(M(:,ip),M(:,jp)),M(:,lp));
            A(ip,jp,lp) = V'*Intx(1:3*N+1)'/D(lp);
            V = conv(conv(M(:,2),M(:,ip)),conv(M(:,jp),M(:,lp)));
            B(ip,jp,lp) = V'*Intx(1:4*N+1)'/D(lp);
        end
    end
    V = M(:,lp)'*Intx(1:N+1)'; 
    C(lp) = V/D(lp); % column vector 
end

%%% Stochastic solution (gPC) %%%
xinit(Np) = 0;
xinit(1) = xhat;
yinit(Np) = 0;
yinit(1) = yhat;
zinit = [xinit,yinit];
[T,Z] = ode15s(@NP_unc,sample_ts,zinit');

%%% determine expected values and confidence limits %%%
Tlen = length(T);
for i=1:Tlen
    % for nutrients
    Xbar(i) = Z(i,1:Np)*M'*Intx(1:Np)'; % mean 
    Xvar(i) = Z(i,1:Np).^2*D' - (Xbar(i))^2; % recall D = [<phi_0^2>, <phi_1^2>, <phi_2^2>, ...]  % variance
    Xsd(i) = sqrt(Xvar(i)); % standard deviation
    % for algae
    Ybar(i) = Z(i,Np+1:2*Np)*M'*Intx(1:Np)'; % mean 
    Yvar(i) = Z(i,Np+1:2*Np).^2*D' - (Ybar(i))^2; % recall D = [<phi_0^2>, <phi_1^2>, <phi_2^2>, ...] % variance
    Ysd(i) = sqrt(Yvar(i)); % standard deviation
end

%_______________________ Make Figures _____________________________________
time = 0:0.05:160; % shortened time horizon for plotting

%%% Figure 4a - Monte Carlo %%%
fig=figure;
set(fig,'defaultAxesColorOrder',[0 0 0; 0 0 0]);
hold on

yyaxis left % nutrients on left axis
hm=plot(time,[xhat mc_mean2(1,:)],'color',nutrientcolor,'linewidth',3);
timevec = [time, fliplr(time)];
fillvec = [[xhat mc_mean2(1,:)]+[0 mc_stdev2(1,:)], fliplr([xhat mc_mean2(1,:)]-[0 mc_stdev2(1,:)])];
hf=fill(timevec, fillvec, sdshading, 'LineStyle','none');
set(hf,'facealpha',.2)
ylim([0,0.3])
yyaxis right % algae on right axis
hm2=plot(time,[yhat mc_mean2(2,:)],'color',algaecolor,'linewidth',3);
fillvec = [[yhat mc_mean2(2,:)]+[0, mc_stdev2(2,:)], fliplr([yhat mc_mean2(2,:)]-[0, mc_stdev2(2,:)])];
hf2=fill(timevec, fillvec, sdshading, 'LineStyle','none');
set(hf2,'facealpha',.2)
ylim([0,0.05])
xlim([0,max(time)])
xlabel('time (days)')
yyaxis left
ylabel('nutrients')
yyaxis right
ylabel('algae')
set(gca,'fontsize',18) 
legend([hm,hm2],'nutrient mean \pm 1 stdev','algae mean \pm 1 stdev','Location','northeast')
legend boxoff                  
NW = [min(xlim) max(ylim)]+[diff(xlim) -diff(ylim)]*0.02;
text(NW(1), NW(2), '(a)', 'VerticalAlignment','top', 'HorizontalAlignment','left','fontsize',18)
title('Monte Carlo, 10,000 samples','fontsize',18,'FontWeight','normal')
set(gcf, 'units','inches','Position',  [2, 2, 6, 4])
set(gcf, 'units','inches','PaperSize', [6 4]); 
%saveas(gcf,'Figure4a.pdf')
hold off

%%% Figure 4b - Intrusive PC %%%
fig=figure;
set(fig,'defaultAxesColorOrder',[0 0 0; 0 0 0]);
hold on

yyaxis left % nutrients on left axis
hm=plot(T, Xbar, 'color',nutrientcolor,'linewidth',3);
fillvec = [[Xbar-Xsd], fliplr([Xbar+Xsd])];
hf=fill([T', fliplr(T')], fillvec, sdshading, 'LineStyle','none');
set(hf,'facealpha',.2)
ylim([0,0.3])

yyaxis right % algae on right axis
hm2=plot(T,Ybar,'color',algaecolor,'linewidth',3);
fillvec = [[Ybar-Ysd], fliplr([Ybar+Ysd])];
hf2=fill([T', fliplr(T')], fillvec, sdshading, 'LineStyle','none');
set(hf2,'facealpha',.2)
ylim([0,0.05])
xlim([0,max(time)])

xlabel('time (days)')
yyaxis left
ylabel('nutrients')
yyaxis right
ylabel('algae')
set(gca,'fontsize',18) 
NW = [min(xlim) max(ylim)]+[diff(xlim) -diff(ylim)]*0.02;
text(NW(1), NW(2), '(b)', 'VerticalAlignment','top', 'HorizontalAlignment','left','fontsize',18)
title('intrusive polynomial chaos','fontsize',18,'FontWeight','normal')
set(gcf, 'units','inches','Position',  [2, 2, 6, 4])
set(gcf, 'units','inches','PaperSize', [6 4]); 
%saveas(gcf,'Figure4b.pdf')
hold off

%%% Figure 4c - Nonintrusive PC %%%

fig=figure;
set(fig,'defaultAxesColorOrder',[0 0 0; 0 0 0]);
hold on

yyaxis left % nutrients on left axis
hm=plot(time,[xhat pce_mean(1,:)],'color',nutrientcolor,'linewidth',3);
timevec = [time, fliplr(time)];
fillvec = [[xhat pce_mean(1,:)]+[0 pce_stdev(1,:)], fliplr([xhat pce_mean(1,:)]-[0 pce_stdev(1,:)])];
hf=fill(timevec, fillvec, sdshading, 'LineStyle','none');
set(hf,'facealpha',.2)
ylim([0,0.3])

yyaxis right % algae on right axis
hm2=plot(time,[yhat pce_mean(2,:)],'color',algaecolor,'linewidth',3);
fillvec = [[yhat pce_mean(2,:)]+[0 pce_stdev(2,:)], fliplr([yhat pce_mean(2,:)]-[0, pce_stdev(2,:)])];
hf2=fill(timevec, fillvec, sdshading, 'LineStyle','none');
set(hf2,'facealpha',.2)
ylim([0,0.05])
xlim([0,max(time)])
xlabel('time (days)')

yyaxis left
ylabel('nutrients')
yyaxis right
ylabel('algae')
set(gca,'fontsize',18) 
NW = [min(xlim) max(ylim)]+[diff(xlim) -diff(ylim)]*0.02;
text(NW(1), NW(2), '(c)', 'VerticalAlignment','top', 'HorizontalAlignment','left','fontsize',18)
title('nonintrusive polynomial chaos','fontsize',18,'FontWeight','normal')
set(gcf, 'units','inches','Position',  [2, 2, 6, 4])
set(gcf, 'units','inches','PaperSize', [6 4]); 
%saveas(gcf,'Figure4c.pdf')
hold off
