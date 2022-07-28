% recreates Figures 2a, b, c


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
bval=b; % 
c=0.8; % gamma
d=0.1; % delta
e=0; % eta
% rho = 0.4; % CAN DELETE?
rhovals = 0.1:0.1:0.4; % variance of beta considered

% Initial conditions 
%xhat = 0.22; % nutrient initial conditions
xhatbig = 0.22; % bigger nutrient initial condition
xhatsmall = 0.16; % smaller nutrient initial condition
yhat = 0.001; % algae initial conditions
tmax = 300; % max time
sample_ts = 0:tmax; % time span for solver

% MC simulation paramters
Nsim = 10000; % number of simulations (10000 takes a long time; can reduce this to 100 to play with code)

%__________________________________________________________________________
% Figure 2a, where N0 < (d/bc)*exp(1/2)

[Tdet,Zdet] = ode45(@NP_det,sample_ts,[xhatsmall;yhat]); 

figure
plot(Tdet,Zdet(:,2),'color',algaecolordet,'linewidth',3);
hold on
% PtotalRDE = zeros(1,length(rhovals)); DELETE THIS?
colscaling = [0.15 0.325 0.6 1]; % for making darker shades of purple

ct = 1;
for rho = rhovals
    b = bval; % ensure b is set correctly 
    uvals = -1 + 2*rand(Nsim,1); % generate vector of random values between [-1,1]
    bvals = b + rho*uvals; % convert this to values of beta
    MCsims = zeros(length(sample_ts),Nsim); % each column will store one simulation
    for j=1:Nsim
        b = bvals(j);
        [Tvals,Zvals] = ode45(@NP_det,sample_ts,[xhatsmall;yhat]); 
        MCsims(:,j) = Zvals(:,2); % algae population is second column of Zdet
    end
    MCmean = mean(MCsims,2);
    hm = plot(Tvals,MCmean,'color',[algaecolorRDE,colscaling(ct)],'Linewidth',3); 
    
    % calculate area under the algae curve, then divide it by time algae
    % live (i.e., by d) to get Ptotal
    %PtotalRDE(ct) = trapz(Tvals,MCmean)*d; % DELETE THIS?
    
    ct = ct+1;
end

hm2 = plot(Tdet,Zdet(:,2),'color',algaecolordet,'linewidth',3);

xlim([0,300])
ax = gca; 
ax.YAxis.Exponent = 0;
ylim([0,max(MCmean)*1.2])
xlabel('time (days)')
ylabel('algae')
set(gca,'fontsize',18) 
lgd = legend([hm2 hm],'constant, \beta','random variable, B','Interpreter','tex','Location','best');
title(lgd,'growth rate:','FontSize',18)
NW = [min(xlim) max(ylim)]+[diff(xlim) -diff(ylim)]*0.02;
text(NW(1), NW(2), '(a)', 'VerticalAlignment','top', 'HorizontalAlignment','left','fontsize',18)
legend boxoff                   
hold off
set(gca,'box','off');
set(gcf, 'units','inches','Position',  [2, 2, 6, 4])
set(gcf, 'units','inches','PaperSize', [6 4]); 
saveas(gcf,'Figure2a.pdf')

%__________________________________________________________________________
% Figure 2b, where N0 > (d/bc)*exp(1/2)

b = bval; 
[Tdet,Zdet] = ode45(@NP_det,sample_ts,[xhatbig;yhat]); 
figure
hold on
% PtotalRDE = zeros(1,length(rhovals)); % DELETE THIS?
ct = 1;
for rho = rhovals
    uvals = -1 + 2*rand(Nsim,1); % generate vector of random values between [-1,1]
    bvals = b + rho*uvals; % convert this to values of b
    sample_ts = [0:tmax]; 
    MCsims = zeros(length(sample_ts),Nsim); % each column will save one simulation
    for j=1:Nsim
        b = bvals(j);
        [Tvals,Zvals] = ode45(@NP_det,sample_ts,[xhatbig;yhat]); 
        MCsims(:,j) = Zvals(:,2); % algae population is second column of Zdet
    end
    MCmean = mean(MCsims,2);
    hm = plot(Tvals,MCmean,'color',[algaecolorRDE,colscaling(ct)],'Linewidth',3);

%     if rho==0.2 % TO DELETE? How is this used in a later figure?
%         maxPmax = max(MCmean); % to plot later in Figure 3
%         Phist = max(MCsims);
%     end
%     if rho==0.2
%         Nhist = MCNvals(end,:);
%     end
    
    % calculate area under the algae curve, then divide it by time algae
    % live (i.e., by d) to get Ptotal
    % PtotalRDE(ct) = trapz(Tvals,MCmean)*d; % DELETE THIS?
    
    ct = ct+1;
end
hm2 = plot(Tdet,Zdet(:,2),'color',algaecolordet,'linewidth',3);
xlim([0,300])
ylim([0,max(Zdet(:,2))*1.2])
xlabel('time (days)')
ylabel('algae')
set(gca,'fontsize',18) 
NW = [min(xlim) max(ylim)]+[diff(xlim) -diff(ylim)]*0.02;
text(NW(1), NW(2), '(b)', 'VerticalAlignment','top', 'HorizontalAlignment','left','fontsize',18)
hold off

set(gcf, 'units','inches','Position',  [2, 2, 6, 4])
set(gcf, 'units','inches','PaperSize', [6 4]); 
saveas(gcf,'Figure2b.pdf')

% disp('compared to area under deterministic algae curve:') DELETE THIS
% LINE AND NEXT?
% trapz(Tdet,Zdet(:,2))

%__________________________________________________________________________
% Figure 2c, nutrient results 

b = bval; 

[Tdet,Zdet] = ode45(@NP_det,sample_ts,[xhatbig;yhat]); 

figure
hold on
MCNvals = []; % clear values/size from before
ct =1; 
for rho = rhovals
    b = bval; % make sure b is set correctly 
    uvals = -1 + 2*rand(Nsim,1); % generate vector of random values between [-1,1]
    bvals = b + rho*uvals; % convert this to values of b
    for j=1:Nsim
        b = bvals(j);
        [Tvals,Zvals] = ode45(@NP_det,sample_ts,[xhatbig;yhat]); %=nut,y=algae
        MCNvals(:,j) = Zvals(:,1); % nutrient population is first column of Zdet
    end
    MCNmean = mean(MCNvals,2);
    hm = plot(Tvals,MCNmean,'color',[nutrientcolorRDE,colscaling(ct)],'Linewidth',3);

    % save values to compare with pdf later in Figure 3 - TO DELETE?
%     if rho==0.2
%         Nhist = MCNvals(end,:);
%     end
    ct = ct+1;
end

hm2 = plot(Tdet,Zdet(:,1),'color',nutrientcolordet,'linewidth',3);

xlim([0,200])
ylim([0, xhatbig*1.2])
xlabel('time (days)')
ylabel('nutrients')

set(gca,'fontsize',18) 
lgd = legend([hm2 hm],'constant, \beta','random variable, B','Interpreter','tex','Location','best');
title(lgd,'growth rate:','FontSize',18)
NW = [min(xlim) max(ylim)]+[diff(xlim) -diff(ylim)]*0.02;
text(NW(1), NW(2), '(c)', 'VerticalAlignment','top', 'HorizontalAlignment','left','fontsize',18)
legend boxoff                   
hold off

set(gcf, 'units','inches','Position',  [2, 2, 6, 4])
set(gcf, 'units','inches','PaperSize', [6 4]); 
saveas(gcf,'Figure2c.pdf')
