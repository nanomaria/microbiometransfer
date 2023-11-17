clear all
close all

% Standard parameters
alphani = 0.0;
Kn = 90;
ri = -0.5;
i0 = 10;

ymax = 6e-3*ones(100,1);


% Vary i0

%i03 = linspace(1,35,10);
i0 = linspace(1,40,100);
%i02 = linspace(1,35,100);
lambdan = (alphani*Kn - ri)./(i0*Kn.*log(i0));
%lambdan2 = (alphani*Kn - ri)./(Kn.*log(i02));
%lambdan3 = (alphani*Kn - ri)./(Kn.*log(i03));


figure(1)
subplot(2,2,1)
area(i0,ymax,'EdgeColor','b','FaceColor','b','FaceAlpha',.3,'EdgeAlpha',.3)
hold on
area(i0,lambdan,'EdgeColor','w','FaceColor','w','FaceAlpha',1, 'EdgeAlpha',1)
hold on
area(i0,lambdan,'Edgecolor','r','FaceColor','r','FaceAlpha',.3,'EdgeAlpha',.3)
hold on
plot(i0,lambdan,'k','linewidth',1.5)
ylim([0 6e-4])
xlim([2,40])
xlabel('Propagule pressure (i_0)')
ylabel('Microbiome transfer rate (\lambda_n)')
set(gca,'fontsize',12)
title('(a)')

% Vary carrying capacity

i0 = 10;

Kn = linspace(80,160,100);
lambdan = (alphani.*Kn - ri)./(i0*Kn.*log(i0));
ymax = 6e-3*ones(size(Kn,2),1);
figure(1)
subplot(2,2,2)
area(Kn,ymax,'EdgeColor','b','FaceColor','b','FaceAlpha',.3,'EdgeAlpha',.3)
hold on
area(Kn,lambdan,'EdgeColor','w','FaceColor','w','FaceAlpha',1, 'EdgeAlpha',1)
hold on
area(Kn,lambdan,'Edgecolor','r','FaceColor','r','FaceAlpha',.3,'EdgeAlpha',.3)
hold on
plot(Kn,lambdan,'k','linewidth',1.5)
ylim([0 6e-4])
xlabel('Carrying capacity (K_n)')
ylabel('Microbiome transfer rate (\lambda_n)')
set(gca,'fontsize',12)
title('(b)')

% Vary competition

Kn = 90;

alphani = linspace(0.02,0.03,100);
lambdan = (alphani.*Kn - ri)./(i0*Kn*log(i0));

figure(1)
subplot(2,2,3)
area(alphani,ymax,'EdgeColor','b','FaceColor','b','FaceAlpha',.3,'EdgeAlpha',.3)
hold on
area(alphani,lambdan,'EdgeColor','w','FaceColor','w','FaceAlpha',1, 'EdgeAlpha',1)
hold on
area(alphani,lambdan,'Edgecolor','r','FaceColor','r','FaceAlpha',.3,'EdgeAlpha',.3)
hold on
plot(alphani,lambdan,'k','linewidth',1.5)
ylim([0 6e-4])
xlabel('Competitive effect (\alpha_{ni})')
ylabel('Microbiome transfer rate (\lambda_n)')
set(gca,'fontsize',12)
title('(c)')

% Vary ri

alphani = 0.02;
ri = linspace(1,1.5,100);
lambdan = (alphani.*Kn - ri)./(i0*Kn*log(i0));

figure(1)
subplot(2,2,4)
area(ri,ymax,'EdgeColor','b','FaceColor','b','FaceAlpha',.3,'EdgeAlpha',.3)
hold on
area(ri,lambdan,'EdgeColor','w','FaceColor','w','FaceAlpha',1, 'EdgeAlpha',1)
hold on
area(ri,lambdan,'Edgecolor','r','FaceColor','r','FaceAlpha',.3,'EdgeAlpha',.3)
hold on
plot(ri,lambdan,'k','linewidth',1.5)
ylim([0 6e-4])
xlabel('Growth rate (r_i)')
ylabel('Microbiome transfer rate (\lambda_n)')
set(gca,'fontsize',12)
title('(d)')
