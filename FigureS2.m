% recreates Figures S2a-f
% requires `NP_det.m` in the working directory
% saves `FigureS2a.pdf` - `FigureS2f.pdf` in working directory

clear 
close all
global a b c d e rho

%_______________________________Define colors______________________________
algaecolordet = 1/255*[118,176,65]; % color for deterministic P results (green)
nutrientcolordet = 1/255*[255,201,20]; % color for deterministic N results (yellow)
algaecolorRDE = 1/255*[125,91,166]; % RDE algae color (purple); alternatively pink: [252, 100, 113]
nutrientcolorRDE = 1/255*[23,190,187]; % RDE nutrient color (blue)
%__________________________________________________________________________

a=0.00075; 
b=1; %beta
bval = b; % will need this to reset b after Monte Carlo sampling it
c=0.8; % gamma
d=0.1; %0.1 delta

% different eta values considered (see Figure caption)
case1e = 0.007; %1.5*a*c/d; 
case2e = 0.00595; %(a*c/d)*(2-(a*c/(4*d^2)))/2.5;
case3e = 0.005; %0.6*a*c/d;

evals=[case1e case2e case3e case1e case2e case3e]; % 3 different eta values, repeated
N0vals = [0.3 0.3 0.3 0.15 0.15 0.15]; % vector of N0 values

yhat = 0.0005; 

tmax = 350; 
sample_ts = 0:tmax;  
Nsim = 100; % number of MC simulations
labels = ['a' 'b' 'c' 'd' 'e' 'f'];

%__________________________________________________________________________
for k = 1:length(N0vals)
    
b = bval; % make sure b is set correctly
e = evals(k);
[Tdet,Zdet] = ode45(@NP_det,sample_ts,[N0vals(k);yhat]); %=nutrients,y=algae
figure
plot(Tdet,Zdet(:,2),'color',algaecolordet,'linewidth',3);
maxval = 0; % will use this to space axes nicely

hold on

for rho = 0.1:0.05:0.4
    b = bval; % ensure b is set correctly 
    uvals = -1 + 2*rand(Nsim,1); % generate vector of random values between [-1,1]
    bvals = b + rho*uvals; % convert this to values of b
    MCsims = zeros(length(sample_ts),Nsim); % each column will save one simulation
    for j=1:Nsim
        b = bvals(j);
        [Tvals,Zvals] = ode45(@NP_det,sample_ts,[N0vals(k);yhat]); %=nut,y=algae
        MCsims(:,j) = Zvals(:,2); % algae population is second column of Zdet
    end
    MCmean = mean(MCsims,2);
    hm = plot(Tvals,MCmean,'color',[algaecolorRDE,2*rho],'Linewidth',3);
    maxval = max(maxval, max(MCmean));
end
hm2 = plot(Tdet,Zdet(:,2),'color',algaecolordet,'linewidth',3);
if k <= 3
    xlim([0,150])
else
    xlim([0,tmax])
end
maxval = max(maxval,max(Zdet(:,2)));
ylim([0, maxval*1.1])
xlabel('time (days)')
ylabel('algae')
set(gca,'fontsize',18) 
if k==1
    legend([hm2 hm],'P','E_B[P]','Interpreter','tex','Location','best')
    legend boxoff                  
end
NW = [min(xlim) max(ylim)]+[diff(xlim) -diff(ylim)]*0.02;
text(NW(1), NW(2), strcat('(', labels(k), ')'), 'VerticalAlignment','top', 'HorizontalAlignment','left','fontsize',18)
hold off

set(gcf, 'units','inches','Position',  [2, 2, 6, 4])
set(gcf, 'units','inches','PaperSize', [6 4]); 
saveas(gcf,['AnalyticPlot_Open' num2str(k) '.pdf'])
end
