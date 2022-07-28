% recreates Figures 3a, b
% saves `Figure3a.pdf` and `Figure3b.pdf` in working directory

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
variance = 0.0033:0.001:0.0533; % variance for rho = 0.1 to 0.4

% Initial conditions 
xhatbig = 0.22; % bigger nutrient initial condition
xhatsmall = 0.16; % smaller nutrient initial condition 
yhat = 0.001; % algae initial conditions

%__________________________________________________________________________

% Figure 3a  
xhat = xhatsmall; 
syms Nstarl Nstaru
syms btilde ntilde
Nstareqn = ntilde == (d/(btilde*c))*log(ntilde/xhat)+yhat/c+xhat;
EPtotal = zeros(1,length(variance)); 
ct = 1;

for v = variance
    rho = sqrt(3*v); % relates rho (spread of uniform distribution) to the variance 
    bl = b-rho;
    bu = b+rho; 
    eqn1 = Nstarl==(d/(bu*c))*log(Nstarl/xhat)+yhat/c+xhat;
    eqn2 = Nstaru==(d/(bl*c))*log(Nstaru/xhat)+yhat/c+xhat;
    Nl = min(double(solve(eqn1,Nstarl)));
    Nu = min(double(solve(eqn2,Nstaru)));
    fnuniform = zeros(100,1);
    fnbeta = zeros(100,1); 
    nrange = linspace(Nl,Nu,100);
    for j = 1:100
        n = nrange(j);
        btildeval = double(vpasolve(subs(Nstareqn,ntilde,n),btilde,[bl*0.999 bu*1.001])); % finds solution in interval [bl bu]
        % added in the 0.999 and 1.001 avoid some rounding issues at the ends
        fnuniform(j) = (1/(2*rho))/abs((d*n/(btildeval*d-c*n*btildeval^2))*log(n/xhat)); % pdf for uniform distribution b
        Efnuniform(j) = n*(1/(2*rho))/abs((d*n/(btildeval*d-c*n*btildeval^2))*log(n/xhat)); % integrand of E[N^*] function
        x = (btildeval-b+rho)/(2*rho); % shift value to go into pdf for beta
        betapdf = (1/(2*rho))*x*(1-x)/beta(2,2); % have to scale by 1/2nu so values are high enough (since normal side is 1, but now it's 2nu)
        fnbeta(j) = betapdf/abs((d*n/(btildeval*d-c*n*btildeval^2))*log(n/xhat)); % pdf for beta distribution
    end
    ENstar = trapz(nrange,Efnuniform); 
    EPtotal(ct) = yhat+c*(xhat-ENstar);
    ct=ct+1;
end

% find deterministic P_total
Nstardet = double(vpasolve(subs(Nstareqn,btilde,b),ntilde,[Nl Nu])); % deterministic Nstar; restricted to Nl and Nu for rho=0.4
Ptotaldet = yhat+c*(xhat-Nstardet);

% plot
figure
hold on
line1 = plot(variance, EPtotal,'-', 'linewidth',3,'color',algaecolorRDE,'markerfacecolor',algaecolorRDE);
line2 = yline(Ptotaldet,':','linewidth',3,'color',algaecolordet);
xlabel('variance of B')
ylabel('total algae, P_{total}')
set(gca,'fontsize',18) 
ylim([0.047 0.062])
xlim([variance(1) variance(end)])
%legend('E[P_{total}]','P_{total}','Interpreter','tex','Location','northeast')
lgd = legend([line2 line1],'constant, \beta','random variable, B','Interpreter','tex','Location','northeast');
title(lgd,'growth rate:','FontSize',18)
legend boxoff
NW = [min(xlim) max(ylim)]+[diff(xlim) -diff(ylim)]*0.02;
text(NW(1), NW(2), '(a)', 'VerticalAlignment','top', 'HorizontalAlignment','left','fontsize',18)
hold off

set(gcf, 'units','inches','Position',  [2, 2, 6, 4])
set(gcf, 'units','inches','PaperSize', [6 4]); 
saveas(gcf,'Figure3a.pdf')

%__________________________________________________________________________
% Figure 3b

xhat = xhatbig;
syms Nstarl Nstaru
syms btilde ntilde
Nstareqn = ntilde == (d/(btilde*c))*log(ntilde/xhat)+yhat/c+xhat;

EPtotal = zeros(1,length(variance)); % E[P_total]
ct = 1;

for v = variance
    rho = sqrt(3*v);
    bl = b-rho;
    bu = b+rho; 
    eqn1 = Nstarl==(d/(bu*c))*log(Nstarl/xhat)+yhat/c+xhat;
    eqn2 = Nstaru==(d/(bl*c))*log(Nstaru/xhat)+yhat/c+xhat;
    Nl = min(double(solve(eqn1,Nstarl)));
    Nu = min(double(solve(eqn2,Nstaru)));
    fnuniform = zeros(100,1);
    fnbeta = zeros(100,1); 
    nrange = linspace(Nl,Nu,100);
    for j = 1:100
        n = nrange(j);
        btildeval = double(vpasolve(subs(Nstareqn,ntilde,n),btilde,[bl*0.999 bu*1.001])); % finds solution in interval [bl bu]
        % added in the 0.999 and 1.001 avoid some rounding issues at the ends
        fnuniform(j) = (1/(2*rho))/abs((d*n/(btildeval*d-c*n*btildeval^2))*log(n/xhat)); % pdf for uniform distribution b
        Efnuniform(j) = n*(1/(2*rho))/abs((d*n/(btildeval*d-c*n*btildeval^2))*log(n/xhat)); % integrand of E[N^*] function
        x = (btildeval-b+rho)/(2*rho); % shift value to go into pdf for beta
        betapdf = (1/(2*rho))*x*(1-x)/beta(2,2); % have to scale by 1/2nu so values are high enough (since normal side is 1, but now it's 2nu)
        fnbeta(j) = betapdf/abs((d*n/(btildeval*d-c*n*btildeval^2))*log(n/xhat)); % pdf for beta distribution
    end
    ENstar = trapz(nrange,Efnuniform); 
    EPtotal(ct) = yhat+c*(xhat-ENstar);
    ct=ct+1;
end

% find deterministic P_total
Nstardet = double(vpasolve(subs(Nstareqn,btilde,b),ntilde,[Nl Nu])); % deterministic Nstar; restricted to Nl and Nu for rho=0.4
Ptotaldet = yhat+c*(xhat-Nstardet);

%%% PLOT %%%
figure
hold on
plot(variance, EPtotal,'-', 'linewidth',3,'color',algaecolorRDE,'markerfacecolor',algaecolorRDE)
yline(Ptotaldet,':','linewidth',3,'color',algaecolordet)
xlabel('variance of B')
ylabel('total algae, P_{total}')
set(gca,'fontsize',18) 
ylim([0.115 0.133])
xlim([variance(1) variance(end)])
NW = [min(xlim) max(ylim)]+[diff(xlim) -diff(ylim)]*0.02;
text(NW(1), NW(2), '(b)', 'VerticalAlignment','top', 'HorizontalAlignment','left','fontsize',18)
hold off

set(gcf, 'units','inches','Position',  [2, 2, 6, 4])
set(gcf, 'units','inches','PaperSize', [6 4]); %Set the paper to have width 5 and height 5.
saveas(gcf,'Figure3b.pdf')
