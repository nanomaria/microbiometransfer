%%%%%%%%%%%%% Code prepared by Maria M. Martignoni
%%%%%%%%%%%%% Date: Jan 29, 2023
%%%%%%%%%%%%% For questions, contact Maria.MartignoniMseya@mail.huji.ac.il
%%%%%%%%%%%%% Associated manuscript: Microbiome transfer from native to invasive species may increase invasion risk and shorten invasion lag


%%% Code produces figure 1 of the manuscript



close all
clear all


global rm  Km  alphanm alphamn eps epsi rn Kn Ki ri alphani alphain


rm = 2.5;
rn = 1.5;
ri = 1.5;
Ki = 60;
Kn = 90;
Km = 90;
alphain = 0.02;
alphani = 0.02;
alphanm = alphani;
alphamn = alphain;




% Extinction
eps = 0;
epsi = 0;
n0=Kn;
i0 = 5;
im0 = 0;
Tfin = 20;
options = odeset('RelTol',1e-4,'AbsTol',1e-6);
[T,Yext] = ode45(@eq_nondim_K, 0:.1:Tfin, [n0; i0;im0], options);


%microbiome exchange
eps = 1e-6;
options = odeset('RelTol',1e-4,'AbsTol',1e-6);
[T,Ymic] = ode45(@eq_nondim_K, 0:.1:Tfin, [n0; i0;im0], options);

figure(2)
subplot(1,2,1)
plot(T,Yext(:,2)+Yext(:,3),'-','Color',[1 0 0],'linewidth',2)
hold on
plot(T,Ymic(:,3)+Ymic(:,2),'-','Color','#0072BD','linewidth',2)
hold on
%plot(T,Ymic(:,3),':','Color','#0072BD','linewidth',1)
hold on
%legend({'I_0','I_m'}, 'Location','west', 'FontSize',12)
xlabel('Time [arbitrary unit]')
ylabel('Introduced population I')
title('(a)')
%axis([0 Tfin 0 110])
legend({'no microbiome transfer', 'microbiome transfer'},'Location','northwest','FontSize',12)
set(gca,'fontsize',14)
axis([0 Tfin 0 40])








%%% Stochastic simulations

%close all
clear all

% model parameters
rn = 1.5;
ri = 1.5;
rm = 2.5;

Kn = 90;
Ki = 60;
Km = 90;

alphain = 0.02;
alphani = 0.02;
alphamn = 0.02;
alphanm = 0.02;



%eps = 1e-6;

dt = 0.1;
T = 20/dt;

n0 = Kn;
im = 0;

Ntrials = 1; %Number of stochastic trials. Decrease for faster results.

N = zeros(T,Ntrials);
I0 = zeros(T,Ntrials);
Im = zeros(T,Ntrials);


i0_vec = linspace(1,40,400); %% put smaller timesteps for faster results
eps_vec = linspace(0,6e-4,600); %% put smaller timesteps for faster results
wtrials = size(i0_vec,2);
ztrials = size(eps_vec,2);


Pt = zeros(ztrials,wtrials);

for z = 1:ztrials
eps = eps_vec(z);

for w = 1:wtrials

i0 = i0_vec(w);

P_est = zeros(1,Ntrials);


for j = 1:Ntrials

for i = 1:T
      

 N(1,j) = n0;
I0(1,j) = i0;
Im(1,j) = 0;

  
epsNI = poissrnd(eps*N(i,j)*I0(i,j)*dt);
delta_N = ((rn*N(i,j) - rn*N(i,j)^2/Kn - alphain*I0(i,j)*N(i,j)   - alphamn*Im(i,j)*N(i,j)));
delta_I0 = ((ri*I0(i,j) - ri*I0(i,j)*(Im(i,j)+I0(i,j))/Ki  - alphani*N(i,j)*I0(i,j)));
delta_Im = ((rm*Im(i,j) - rm*Im(i,j)*(Im(i,j)+I0(i,j))/Km - alphanm*N(i,j)*Im(i,j)));



N(i+1,j) = N(i,j) + delta_N*dt;
I0(i+1,j) = I0(i,j) + delta_I0*dt - epsNI*dt;
Im(i+1,j) = Im(i,j) + delta_Im*dt + epsNI*dt;



end


if Im(end,j) > 50
P_est(j) = 1;
else 
P_est(j) = 0;
end


end

Pt(z,w) = sum(P_est)/Ntrials;


%figure(3)
%plot(0:T,N,'Color',[0.5,0.5,0.5],'linewidth',1);
%hold on
%h1 =plot(0:T,I0+Im,'--','Color',[0.5,0.5,0.5],'linewidth',1);
%hold on
% legend({'stochastic trajectories','deterministic'},'Location','west', 'FontSize',12)





end

end


% Pt : Matrix with probability of establishment

figure(2)
subplot(1,2,2)
image(i0_vec,eps_vec,Pt,'CDataMapping','scaled')
cMap = interp1([0;1],[1 0 0;0 0 1],linspace(0,1,256));
colorbar
clim([0 1])
ylabel(colorbar,'Probability of establishment','FontSize',14,'Rotation',90);
%c = hot(10);
colormap(cMap);
hColourbar.Label.Position(1) = 10;
ylabel('Microbiome transfer rate (\lambda_n)')
xlabel('Propagule pressure (i_0)')
title('(b)')
set(gca,'fontsize',14)
set(gca,'YDir','normal')  


hold on
i0_cont = linspace(1.01,43,99);
lambdan = (-ri+alphani*Kn)./(Kn.*i0_cont.*log(i0_cont));
plot(i0_cont,lambdan,'k-','Linewidth',1.5)


