%%%%%%%%%%%%% Code prepared by Maria M. Martignoni
%%%%%%%%%%%%% Date: Jan 29, 2023
%%%%%%%%%%%%% For questions, contact Maria.MartignoniMseya@mail.huji.ac.il
%%%%%%%%%%%%% Associated manuscript: Microbiome transfer from native to invasive species may increase invasion risk and shorten invasion lag

%Code supports findings shown in Fig. 3 of the manuscript and shows spatial
%dynamics of microbiome exchange

clear all
close all

% SPACE AND TIME
x1 = -100;
x2 = 100;
dx = 1;
dt = 0.001;
x = x1:dx:x2;


 
% PARAMETERS FOR POPULATION GROWTH
% model parameters
rn = 1.5;
Kn = 80;
alphain = 0.2;
alphani = 0.01;
ri = 1.5;
Ki = 20;
Km = 80;
alphanm = alphani;
alphamn = alphain;
rm = ri;


eps = 1e-5;


Tfin = 120;
t = 1:dt:Tfin;
% INITAL DENSITY
    N0  = Kn*((x>=-80).*(x<=-50))+Kn*((x>=0).*(x<=70));
    I00 = 1*((x>=-80).*(x<=-70))+1*((x>=0).*(x<=10));  
    I0m = 0*  ((x>=-100).*(x<=100)); 

%% Initial conditions
N(:,1)  = N0;  % x positions at time zero
I0(:,1) = I00;
Im(:,1) = I0m;



Dm = 0.1*((x>=-70).*(x<=-50))+0.1*((x>=0).*(x<=70));
Di = 0.1*((x>=-70).*(x<=-50))+0.1*((x>=0).*(x<=70));
   

for  j = 1:length(t)
 %% Boundary conditions (no-flux)
 Im(1,j) = Im(2,j);
 Im(length(x),j) = Im(length(x)-1,j);
 N(1,j) = N(2,j);
 N(length(x),j) = N(length(x)-1,j);
 I0(1,j) = I0(2,j);
 I0(length(x),j) = I0(length(x)-1,j);
 

for i = 2:(length(x)-1)

epsNI = poissrnd(eps*N(i,j)*I0(i,j)*dt);

 N(i,j+1)= N(i,j) +  (rn* N(i,j)*(1-      N(i,j)     /Kn) - alphain*I0(i,j)*N(i,j) -alphamn*Im(i,j)*N(i,j))*dt;
I0(i,j+1)= I0(i,j) + (Di(i)*((I0(i+1,j) -2*I0(i,j) + I0(i-1,j))/(dx)^2)+ri*I0(i,j)*(1-(I0(i,j)+Im(i,j))/Ki) - alphani*I0(i,j)*N(i,j) - epsNI)*dt;
Im(i,j+1)= Im(i,j) + (Dm(i)*((Im(i+1,j) -2*Im(i,j) + Im(i-1,j))/(dx)^2)+rm*Im(i,j)*(1-(Im(i,j)+I0(i,j))/Km) - alphanm*Im(i,j)*N(i,j) + epsNI)*dt;     

% Dm*((Im(i+1,j) -2*Im(i,j) + Im(i-1,j))/(dx)^2) 


end
end


%left_color = [0 0 0];
%right_color = [0 0 0];
%set(fig4,'defaultAxesColorOrder',[left_color; right_color]);


figure(3)
for k = 1:(length(t)*dt)
plot(x,N(:,round(k/dt)),'-', 'Color','#EDB120','linewidth', 1.5)
hold on 
plot(x,I0(:,round(k/dt)),'--', 'Color','#0072BD','linewidth', 1.5)
hold on 
plot(x,Im(:,round(k/dt)),':', 'Color','#0072BD','linewidth', 2)
xlabel('Distance [arbitrary unit])')
ylabel('Population')
legend('N','I_0','I_m')
hold off
drawnow
end

% Just plot (without drawing)

figure(1)
subplot(2,2,1)
plot(x,N(:,1),'-', 'Color','#EDB120','linewidth', 1.5)
hold on 
plot(x,I0(:,1),'--', 'Color','#0072BD','linewidth', 1.5)
hold on 
plot(x,Im(:,1),':', 'Color','#0072BD','linewidth', 2)
xlabel('Distance [arbitrary unit])')
ylim([0,100])
ylabel('Population')
legend('N','I_0','I_m')
set(gca,'fontsize',14)

figure(1)
subplot(2,2,2)
plot(x,N(:,round(end/4)),'-', 'Color','#EDB120','linewidth', 1.5)
hold on 
plot(x,I0(:,round(end/4)),'--', 'Color','#0072BD','linewidth', 1.5)
hold on 
plot(x,Im(:,round(end/4)),':', 'Color','#0072BD','linewidth', 2)
xlabel('Distance [arbitrary unit])')
ylabel('Population')
ylim([0,100])
%legend('N','I_0','I_m')
set(gca,'fontsize',14)

figure(1)
subplot(2,2,3)
plot(x,N(:,round(end/3*2)),'-', 'Color','#EDB120','linewidth', 1.5)
hold on 
plot(x,I0(:,round(end/3*2)),'--', 'Color','#0072BD','linewidth', 1.5)
hold on 
plot(x,Im(:,round(end/3*2)),':', 'Color','#0072BD','linewidth', 2)
xlabel('Distance [arbitrary unit])')
ylabel('Population')
ylim([0,100])
set(gca,'fontsize',14)
%legend('N','I_0','I_m')

figure(1)
subplot(2,2,4)
plot(x,N(:,end),'-', 'Color','#EDB120','linewidth', 1.5)
hold on 
plot(x,I0(:,end),'--', 'Color','#0072BD','linewidth', 1.5)
hold on 
plot(x,Im(:,end),':', 'Color','#0072BD','linewidth', 2)
xlabel('Distance [arbitrary unit])')
ylabel('Population')
ylim([0,100])
%legend('N','I_0','I_m')
set(gca,'fontsize',14)



