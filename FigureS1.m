% recreates Figures S1a,b,c
% requires `NP_det.m` in the working directory
% saves `FigureS1a.pdf`, `FigureS1b.pdf`, and `FigureS1c.pdf` in working directory

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
bval=b;
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
% Figure S1a 

rho = 0.2; 
bl = b-rho;
bu = b+rho; 
bvec = linspace(bu,bl,100);
fbuniform = ones(100,1)*(1/(2*rho));
x = (bvec-b+rho)/(2*rho); % shift value to go into pdf for beta
fbbeta = (1/(2*rho))*x.*(1-x)./beta(2,2); % have to scale by 1/2nu so values are high enough (since normal side is 1, but now it's 2nu)
figure
h1 = plot(bvec,fbuniform,'k','linewidth',3);
xlim([bl-.1, bu+.1])
hold on
plot([bl bl], [0 fbuniform(1)],'k','linewidth',3)
plot([bu bu], [0 fbuniform(end)],'k','linewidth',3)
h2 = plot(bvec,fbbeta,'k:','linewidth',3);
ylim([0 max(fbbeta)*1.3])
xlabel('b')
set(gca,'fontsize',18) 
NW = [min(xlim) max(ylim)]+[diff(xlim) -diff(ylim)]*0.02;
text(NW(1), NW(2), '(a)', 'VerticalAlignment','top', 'HorizontalAlignment','left','fontsize',18)
legend([h1 h2],'uniform','beta','Location','best')
legend boxoff                   
box off
hold off

set(gcf, 'units','inches','Position',  [2, 2, 6, 4])
set(gcf, 'units','inches','PaperSize', [6 4]); 
saveas(gcf,'FigureS1a.pdf')

%__________________________________________________________________________
% Figure S1b

bl = b-rho;
bu = b+rho; 
bvec = bl:0.05:bu;
Pmax = (d./bvec).*(log(d./(bvec.*c*xhat))-1) + yhat + c*xhat;
Pmaxl = -d/bl+(d/bl)*log(d/(bl*c*xhat)) + yhat + c*xhat; 
Pmaxu = -d/bu+(d/bu)*log(d/(bu*c*xhat)) + yhat + c*xhat; 
syms btilde ptilde
Pmaxeqn = ptilde==-d/btilde+(d/btilde)*log(d/(btilde*c*xhat)) + yhat + c*xhat; 

fpuniform = zeros(100,1);
fpbeta = zeros(100,1); 
prange = linspace(Pmaxl,Pmaxu,100);
for j = 1:100
    p = prange(j);
    btildeval = double(vpasolve(subs(Pmaxeqn,ptilde,p),btilde,[bl*0.99 bu*1.01])); % finds solution in interval [bl bu]
    fpuniform(j) = (1/(2*rho))/abs(-(d/btildeval^2)*log(d/(btildeval*c*xhat))); % pdf for uniform distribution b
    x = (btildeval-b+rho)/(2*rho); % shift value to go into pdf for beta
    betapdf = (1/(2*rho))*x*(1-x)/beta(2,2); % have to scale by 1/2nu so values are high enough (since normal side is 1, but now it's 2nu)
    fpbeta(j) = betapdf/abs(-(d/btildeval^2)*log(d/(btildeval*c*xhat))); % pdf for beta distribution
end
figure
hold on

%simulate and then make a histogram from simulations to check that the pdf looks ok
Nsim = 10000; % number of MC simulations; should be 10,000 for full plot
sample_ts = 0:tmax; 
MCNvals = zeros(length(sample_ts),Nsim);
rhovals = 0.1:0.1:0.4; % variance of beta considered
b = bval; % ensure b is set correctly
[Tdet,Zdet] = ode45(@NP_det,sample_ts,[xhat;yhat]); %=nut,y=algae
figure
hold on
PtotalRDE = zeros(1,length(rhovals));
ct = 1;

uvals = -1 + 2*rand(Nsim,1); % generate vector of random values between [-1,1]
bvals = b + rho*uvals; % convert this to values of b
sample_ts = [0:tmax]; 
MCsims = zeros(length(sample_ts),Nsim); % each column will save one simulation
for j=1:Nsim
    b = bvals(j);
    [Tvals,Zvals] = ode45(@NP_det,sample_ts,[xhat;yhat]); %=nut,y=algae
    MCsims(:,j) = Zvals(:,2); % algae population is second column of Zdet
    MCNvals(:,j) = Zvals(:,1); % nutrient population is first column of Zdet
end
MCmean = mean(MCsims,2);
hm = plot(Tvals,MCmean,'color',[algaecolorRDE,2*rho],'Linewidth',3);
    
% save values for historgrams
maxPmax = max(MCmean); 
Phist = max(MCsims);
Nhist = MCNvals(end,:);
    
% make histogram of Pmax values
histogram(Phist,'Normalization','pdf','FaceColor',algaecolorRDE,'EdgeColor','none','FaceAlpha',0.2)
histogram(Phist,'Normalization','pdf','DisplayStyle','stairs','EdgeColor',algaecolorRDE)

h1 = plot(prange,fpuniform,'k','linewidth',3);
hold on
plot([Pmaxl Pmaxl], [0 fpuniform(1)],'k','linewidth',3)
plot([Pmaxu Pmaxu], [0 fpuniform(end)],'k','linewidth',3)
h2 = plot(prange,fpbeta,'k:','linewidth',3);
plot([Pmaxl Pmaxl], [0 fpbeta(1)],'k:','linewidth',3)
plot([Pmaxu Pmaxu], [0 fpbeta(end)],'k:','linewidth',3)
ylim([0 max(max([fpbeta fpuniform]))*1.3])
xlim([Pmaxl-.002 Pmaxu+.002])
xlabel('P_{max}(B)')
set(gca,'fontsize',18) 
NW = [min(xlim) max(ylim)]+[diff(xlim) -diff(ylim)]*0.02;
text(NW(1), NW(2), '(b)', 'VerticalAlignment','top', 'HorizontalAlignment','left','fontsize',18)
Pmaxdet = (d/b)*(log(d/(b*c*xhat))-1) + yhat + c*xhat; % deterministic Pmax
EPmax = (-d/(2*rho))*log((b+rho)/(b-rho))*(1-log(d/(c*xhat*sqrt(b^2-rho^2)))) + yhat+c*xhat;
plot(Pmaxdet,range(ylim)*0.05, 'ko','MarkerFaceColor','black') % Pmax(beta)
plot(EPmax, range(ylim)*0.02, 'ko') % E[Pmax(B)]
plot(maxPmax, range(ylim)*0.02, 'ks') %max_t (E[P])
hold off

set(gcf, 'units','inches','Position',  [2, 2, 6, 4])
set(gcf, 'units','inches','PaperSize', [6 4]); 
saveas(gcf,'FigureS1b.pdf')


%__________________________________________________________________________
% Figure S1c

syms Nstarl Nstaru
eqn1 = Nstarl==(d/(bu*c))*log(Nstarl/xhat)+yhat/c+xhat;
eqn2 = Nstaru==(d/(bl*c))*log(Nstaru/xhat)+yhat/c+xhat;
Nl = min(double(solve(eqn1,Nstarl)));
Nu = min(double(solve(eqn2,Nstaru)));
syms btilde ntilde
Nstareqn = ntilde == (d/(btilde*c))*log(ntilde/xhat)+yhat/c+xhat;

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
figure
hold on

%histogram from simulations above to check it looks ok
histogram(Nhist,'Normalization','pdf','FaceColor',nutrientcolorRDE,'EdgeColor','none','FaceAlpha',0.2)
histogram(Nhist,'Normalization','pdf','DisplayStyle','stairs','EdgeColor',nutrientcolorRDE)

h1 = plot(nrange,fnuniform,'k','linewidth',3);
xlim([Nl-.01 Nu+.01])
hold on
plot([Nl Nl], [0 fnuniform(1)],'k','linewidth',3)
plot([Nu Nu], [0 fnuniform(end)],'k','linewidth',3)
h2 = plot(nrange,fnbeta,'k:','linewidth',3);
ylim([0 max(fnuniform)*1.3])
xlabel('N^*(B)')
set(gca,'fontsize',18) 
NW = [min(xlim) max(ylim)]+[diff(xlim) -diff(ylim)]*0.02;
text(NW(1), NW(2), '(c)', 'VerticalAlignment','top', 'HorizontalAlignment','left','fontsize',18)
Nstardet = double(vpasolve(subs(Nstareqn,btilde,b),ntilde,[Nl Nu])); % deterministic Nstar
ENstar = trapz(nrange,Efnuniform); 
plot(Nstardet, range(ylim)*0.02, 'ko','MarkerFaceColor','black')
plot(ENstar, range(ylim)*0.02, 'ko') % E[Nstar]
histogram(Nhist,'Normalization','probability','FaceColor',nutrientcolorRDE,'EdgeColor','none','FaceAlpha',0.2)
histogram(Nhist,'Normalization','probability','DisplayStyle','stairs','EdgeColor',nutrientcolorRDE)
hold off

set(gcf, 'units','inches','Position',  [2, 2, 6, 4])
set(gcf, 'units','inches','PaperSize', [6 4]); 
saveas(gcf,'FigureS1c.pdf')

