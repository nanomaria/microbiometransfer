close all
clear all

% model parameters
rn = 1.5;
Kn = 100;
alphain = 0.02;
alphani = 0.01;
alphamn = 0.02;
alphanm = 0.01;
ri = 1.5;
rm = 1.5;
Ki = 20;
Km = 90; 




dt = 0.01;
T = 50/dt;

n0 = Kn*ri*(rn - alphain*Ki)/(rn*ri-alphain*alphani*Kn*Ki);
i0 = Ki*rn*(ri-alphani*Kn)/(rn*ri-alphain*alphani*Kn*Ki);
im0 = 0;


Ntrials = 500;
N = zeros(T,Ntrials);
I0 = zeros(T,Ntrials);
Im = zeros(T,Ntrials);


%i0_vec = 1:10:Ki+1;
eps = 1e-3;%,2e-6,3e-6,4e-6,5e-6,6e-6];



%% Stochastic simulations

for j = 1:Ntrials

for i = 1:T
      

 N(1,j) = n0;
I0(1,j) = i0;
Im(1,j) = im0;

  
epsNI = poissrnd(eps*N(i,j)*I0(i,j)*dt);

delta_N = ((rn*N(i,j) - rn*N(i,j)^2/Kn - alphain*I0(i,j)*N(i,j)   - alphamn*Im(i,j)*N(i,j)));
delta_I0 = ((ri*I0(i,j) - ri*I0(i,j)*(Im(i,j)+I0(i,j))/Ki  - alphani*N(i,j)*I0(i,j)));
delta_Im = ((rm*Im(i,j) - rm*Im(i,j)*(Im(i,j)+I0(i,j))/Km - alphanm*N(i,j)*Im(i,j)));



N(i+1,j) = N(i,j) + delta_N*dt;
I0(i+1,j) = I0(i,j) + delta_I0*dt - epsNI*dt ;
Im(i+1,j) = Im(i,j) + delta_Im*dt + epsNI*dt ;



end



end



detN = zeros(T,1);
detI0 = zeros(T,1);
detIm = zeros(T,1);


for j = 1:1

for i = 1:T
      

 detN(1,j) = n0;
detI0(1,j) = i0;
detIm(1,j) = im0;



epsNI = (eps*detN(i,j)*detI0(i,j)*dt);
delta_N = ((rn*detN(i,j) - rn*detN(i,j)^2/Kn - alphain*detI0(i,j)*N(i,j)   - alphamn*detIm(i,j)*detN(i,j)));
delta_I0 = ((ri*detI0(i,j) - ri*detI0(i,j)*(detIm(i,j)+detI0(i,j))/Ki  - alphani*detN(i,j)*detI0(i,j)));
delta_Im = ((rm*detIm(i,j) - rm*detIm(i,j)*(detIm(i,j)+detI0(i,j))/Km - alphanm*detN(i,j)*detIm(i,j)));



detN(i+1,j) = detN(i,j) + delta_N*dt;
detI0(i+1,j) = detI0(i,j) + delta_I0*dt - epsNI*dt;
detIm(i+1,j) = detIm(i,j) + delta_Im*dt + epsNI*dt;



end



end


%see figure with all the trials




figure(2)
h1 = plot(0:dt:T*dt,N,'-','Color',[0.8,0.8,0.8],'linewidth',1);
hold on
h2 = plot(0:dt:T*dt,I0+Im,'-','Color',[0.8,0.8,0.8],'linewidth',1);
hold on
h3 = plot(0:dt:T*dt,mean(N,2),':','Color','#EDB120','linewidth',1.5);
hold on
h4 = plot(0:dt:T*dt,mean(I0,2) + mean(Im,2),':','Color','#0072BD','linewidth',1.5);
hold on
h5 = plot(0:dt:T*dt,detN,'-','Color','#EDB120', 'linewidth',2);
hold on
h6 = plot(0:dt:T*dt,detI0+detIm,'-','Color','#0072BD', 'linewidth',2);
hold on
%legend('\lambda_m = 5e-2','\lambda_m = 5e-3','\lambda_m = 0','','','','Location','northwest')
lgd = legend([h5,h3,h6,h4],{'deterministic trajectory (N)','mean of stochastic trajectories (N)','deterministic trajectory (I)','mean of stochastic trajectories (I)'});%,'Location','west', 'FontSize',12}) )
set(lgd,'Location','east', 'FontSize',12);
%plot(0:T,detIm,':','Color','#0072BD','linewidth',2
xlabel('Time [arbitrary unit]')
ylabel('Population')
ylim([0 100])
set(gca,'fontsize',12)




