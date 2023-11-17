%%%%%%%%%%%%% Code prepared by Maria M. Martignoni
%%%%%%%%%%%%% Date: Jan 29, 2023
%%%%%%%%%%%%% for questions, contact Maria.MartignoniMseya@mail.huji.ac.il
%%%%%%%%%%%%% Associated manuscript:
%%%%%%%%%%%%% Microbiome transfer from natives to introduced species facilitates invasion after lag time

%%% Code produces figure 1 of the manuscript

close all
clear all


global rn Kn alphain alphani  ri Ki rm Km alphanm alphamn eps epsi


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


% microbiome exchange rate
eps1 = 1e-7; %between species
epsi1 = 0;  %between invaders

Tfin = 50;

n0 = 20;
i00 = 2;
im0 = 0;


%%% First figure (Fig. 1a), invasion lag


% Find ss for initial conditions
eps = 0;
epsi = 0;

options = odeset('RelTol',1e-4,'AbsTol',1e-6);
[T,Yss] = ode45(@eq_nondim_K, 0:1:Tfin, [n0; i00;0], options);

nss = Yss(end,1);
i0ss = Yss(end,2);

eps = eps1;
epsi = epsi1;

n0 = nss;
i00 = i0ss;
im0 = 0;

options = odeset('RelTol',1e-4,'AbsTol',1e-6);
[T,Y1] = ode45(@eq_nondim_K, 0:1:Tfin, [nss; i0ss;im0], options);
Ym = Y1(:,2)+Y1(:,3);


% Find inflection point
a = diff(diff(Ym));
b = zeros(1,size(a,1));
for  w = 1:(size(a,1)-2)
b(w) = (a(w)>0 && a(w+1)<=0);
end
b= b(10:40);
u = find(b==1);
t_lag = u(1)+10;
Ym_a = Ym(1:t_lag);



figure(1)
subplot(2,2,[1,2])
area(0:t_lag-1,Ym_a,'EdgeColor','y','FaceColor','y','FaceAlpha',.3,'EdgeAlpha',.3)
hold on
h1 = plot(T,Y1(:,1),'Color','#EDB120', 'linewidth',2);
hold on
h2 = plot(T,Y1(:,2),'--','Color','#0072BD','linewidth',1);
hold on
h3 = plot(T,Y1(:,3),':','Color','#0072BD','linewidth',1);
hold on
h4 = plot(T,Ym,'Color','#0072BD','linewidth',2);
hold on
legend([h1,h2,h3,h4],{'N','I_0','I_m','I (total)'}, 'Location','west', 'FontSize',10)
xlabel('Time [arbitrary unit]')
ylabel('Population')
title('(a)')
axis([0 Tfin 0 110])
set(gca,'fontsize',12)



%%% Figures b and c: impact of horizontal and vertical microbiome transfer



dt = 0.01;
T = 50/dt;

n0 = Kn*ri*(rn - alphain*Ki)/(rn*ri-alphain*alphani*Kn*Ki);
i0 = Ki*rn*(ri-alphani*Kn)/(rn*ri-alphain*alphani*Kn*Ki);
im0 = 0;


Ntrials = 1;
N = zeros(T,Ntrials);
I0 = zeros(T,Ntrials);
Im = zeros(T,Ntrials);


eps_vec = [1e-7,1e-5];
ri_vec = 1.4;
epsi_vec = [0.02,0];


for g = 1:size(eps_vec,2)
    eps = eps_vec(g);
   
for z = 1:size(epsi_vec,2)

epsi = epsi_vec(z);




for j = 1:Ntrials

for i = 1:T
      

 N(1,j) = n0;
I0(1,j) = i0;
Im(1,j) = im0;

  
epsNI = poissrnd(eps*N(i,j)*I0(i,j));
epsII = poissrnd(epsi*Im(i,j)*I0(i,j));

delta_N = ((rn*N(i,j) - rn*N(i,j)^2/Kn - alphain*I0(i,j)*N(i,j)   - alphamn*Im(i,j)*N(i,j)));
delta_I0 = ((ri*I0(i,j) - ri*I0(i,j)*(Im(i,j)+I0(i,j))/Ki  - alphani*N(i,j)*I0(i,j)));
delta_Im = ((rm*Im(i,j) - rm*Im(i,j)*(Im(i,j)+I0(i,j))/Km - alphanm*N(i,j)*Im(i,j)));



N(i+1,j) = N(i,j) + delta_N*dt;
I0(i+1,j) = I0(i,j) + delta_I0*dt - epsNI*dt - epsII*dt;
Im(i+1,j) = Im(i,j) + delta_Im*dt + epsNI*dt + epsII*dt;



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



epsNI = (eps*detN(i,j)*detI0(i,j));
epsII = (epsi*detI0(i,j)*detIm(i,j));
delta_N = ((rn*detN(i,j) - rn*detN(i,j)^2/Kn - alphain*detI0(i,j)*N(i,j)   - alphamn*detIm(i,j)*detN(i,j)));
delta_I0 = ((ri*detI0(i,j) - ri*detI0(i,j)*(detIm(i,j)+detI0(i,j))/Ki  - alphani*detN(i,j)*detI0(i,j)));
delta_Im = ((rm*detIm(i,j) - rm*detIm(i,j)*(detIm(i,j)+detI0(i,j))/Km - alphanm*detN(i,j)*detIm(i,j)));



detN(i+1,j) = detN(i,j) + delta_N*dt;
detI0(i+1,j) = detI0(i,j) + delta_I0*dt - epsNI*dt -epsII*dt;
detIm(i+1,j) = detIm(i,j) + delta_Im*dt + epsNI*dt + epsII*dt;



end



end

if z == 2
det_1 = detI0+detIm;
% Find inflection point
a = diff(diff(det_1));
b = zeros(1,size(a,1));
for  w = 1:(size(a,1)-2)
b(w) = (a(w)>0 && a(w+1)<=0);
end
u = find(b==1);
t_lag_1 = u(1);
det_1 = det_1(1:t_lag_1+1);
end

% (N, I0+Im are to see the stochastic trajectories)

figure(1)
subplot(2,2,3)
if z == 2 && g == 1
area(0:dt:t_lag_1*dt,det_1,'EdgeColor','y','FaceColor','y','FaceAlpha',.3,'EdgeAlpha',.3)
hold on
end
%plot(0:T,N,'Color',[0.5,0.5,0.5],'linewidth',1);
%hold on
%h1 =plot(0:T,I0+Im,'--','Color',[0.5,0.5,0.5],'linewidth',z);
%hold on
%plot(0:T,mean(N,2),'Color','r','linewidth',1);
%hold on
%plot(0:T,mean(I0,2) + mean(Im,2),':','Color','r','linewidth',z)
%h2 = plot(0:T,detN,'Color','#EDB120', 'linewidth',2);
%hold on
%area(0:t_lag_2,det_2,'EdgeColor','y','FaceColor','y','FaceAlpha',.3,'EdgeAlpha',.3)
if g == 1
f1 = plot(0:dt:T*dt,detI0+detIm,'-','Color','#0072BD', 'linewidth',z);
elseif g == 2
f2 = plot(0:dt:T*dt,detI0+detIm,'-','Color',"#4DBEEE", 'linewidth',z);
legend([f1,f2],'Large \lambda_n', 'Small \lambda_n','Fontsize',10)%Small \lambda_n','','Large \lambda_n', '', 'location', 'southeast','Fontsize',10) %,['Increased horizontal transfer' newline 'among introduced individuals'],['Decreased growth rate' newline '(vertical transmission)'] 
end
hold on
title('(b)')
%legend([h1,h2],['Large microbiome' newline ' transfer rate (\lambda_n)'], ['Small microbiome' newline 'transfer rate (\lambda_n)'],'Fontsize',10)%Small \lambda_n','','Large \lambda_n', '', 'location', 'southeast','Fontsize',10) %,['Increased horizontal transfer' newline 'among introduced individuals'],['Decreased growth rate' newline '(vertical transmission)'] 
      %,'decreased growth rate (vertical transmission)','increased interspecific horizontal transmission','das','das')
xlabel('Time [arbitrary unit]')
ylabel('Introduced population (I_0 + I_m)')
set(gca,'fontsize',12)

end
end



% Plot now for smaller growth rate r

for g = 1:size(eps_vec,2)
    eps = eps_vec(g);
 
epsi =0 ;
for w = 1:size(ri_vec,2)

ri = ri_vec(w);
rm = ri;

n0 = Kn*ri*(rn - alphain*Ki)/(rn*ri-alphain*alphani*Kn*Ki);
i0 = Ki*rn*(ri-alphani*Kn)/(rn*ri-alphain*alphani*Kn*Ki);
im0 = 0;

for j = 1:Ntrials

for i = 1:T
      

 N(1,j) = n0;
I0(1,j) = i0;
Im(1,j) = im0;

  
epsNI = poissrnd(eps*N(i,j)*I0(i,j));
epsII = poissrnd(epsi*Im(i,j)*I0(i,j));

delta_N = ((rn*N(i,j) - rn*N(i,j)^2/Kn - alphain*I0(i,j)*N(i,j)   - alphamn*Im(i,j)*N(i,j)));
delta_I0 = ((ri*I0(i,j) - ri*I0(i,j)*(Im(i,j)+I0(i,j))/Ki  - alphani*N(i,j)*I0(i,j)));
delta_Im = ((rm*Im(i,j) - rm*Im(i,j)*(Im(i,j)+I0(i,j))/Km - alphanm*N(i,j)*Im(i,j)));



N(i+1,j) = N(i,j) + delta_N*dt;
I0(i+1,j) = I0(i,j) + delta_I0*dt - epsNI*dt - epsII*dt;
Im(i+1,j) = Im(i,j) + delta_Im*dt + epsNI*dt + epsII*dt;



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




epsNI = (eps*detN(i,j)*detI0(i,j));
epsII = (epsi*detI0(i,j)*detIm(i,j));
delta_N = ((rn*detN(i,j) - rn*detN(i,j)^2/Kn - alphain*detI0(i,j)*N(i,j)   - alphamn*detIm(i,j)*detN(i,j)));
delta_I0 = ((ri*detI0(i,j) - ri*detI0(i,j)*(detIm(i,j)+detI0(i,j))/Ki  - alphani*detN(i,j)*detI0(i,j)));
delta_Im = ((rm*detIm(i,j) - rm*detIm(i,j)*(detIm(i,j)+detI0(i,j))/Km - alphanm*detN(i,j)*detIm(i,j)));



detN(i+1,j) = detN(i,j) + delta_N*dt;
detI0(i+1,j) = detI0(i,j) + delta_I0*dt - epsNI*dt -epsII*dt;
detIm(i+1,j) = detIm(i,j) + delta_Im*dt + epsNI*dt + epsII*dt;



end



end


%see figure with all the trials

%figure(1)
%subplot(2,3,[4,5])
%if g == 1
%plot(0:dt:T*dt,detI0+detIm,':','Color','#0072BD', 'linewidth',1.5)
%elseif g == 2
%plot(0:dt:T*dt,detI0+detIm,':','Color','#77AC30', 'linewidth',1.5)
%end
%legend('',['Early microbiome transfer' newline 'from natives'], '',['Late microbiome transfer' newline 'from natives'],'Fontsize',10)%Small \lambda_n','','Large \lambda_n', '', 'location', 'southeast','Fontsize',10) %,['Increased horizontal transfer' newline 'among introduced individuals'],['Decreased growth rate' newline '(vertical transmission)'] 
      %,'decreased growth rate (vertical transmission)','increased interspecific horizontal transmission','das','das')
%legend({'',['Early microbiome' newline 'transfer'],'',['Late microbiome' newline 'transfer'] ... %,['Increased horizontal transfer' newline 'among introduced individuals'],['Decreased growth rate' newline '(vertical transmission)'] 
%     }, 'location', 'southeast','Fontsize',10) %,'decreased growth rate (vertical transmission)','increased interspecific horizontal transmission','das','das')
%xlabel('Time [arbitrary unit]')
%ylabel('Introduced population I')
%set(gca,'fontsize',12)

end
end


%%% Figure % reduction in invasion lag time


Tfin = 250;



ri_v = [1.1,1.4,1.7,2.0];
epsi_v = linspace(0,0.04,8);
eps =  1e-7;

Msteps = size(epsi_v,2);
Rsteps = size(ri_v,2);
T_eps = zeros(size(ri_v,2),size(epsi_v,2));



for r = 1:Rsteps


ri = ri_v(r);
rm = ri;

n0 = Kn*ri*(rn - alphain*Ki)/(rn*ri-alphain*alphani*Kn*Ki);
i0 = Ki*rn*(ri-alphani*Kn)/(rn*ri-alphain*alphani*Kn*Ki);
im0 = 0;




for i = 1:size(epsi_v,2)

epsi = epsi_v(i);

[T,Y1] = ode45(@eq_epsi, 0:.1:Tfin, [n0; i0;im0]);

Yi = Y1(:,2)+Y1(:,3);
Yir = round(Yi);

%figure(2)
%plot(T,Y1(:,1),'Color',[0.5,0.5,0.5],'linewidth',1);
%hold on
%h1 =plot(T,Yi,'--','Color',[0.5,0.5,0.5],'linewidth',1);
%legend({'trajectories'},'Location','west', 'FontSize',12)


% Find inflection points
%a = diff(diff(Yi));
%b = zeros(1,size(a,1)-2);
%for  w = 1:(size(a,1)-2)
%b(w) = (a(w)>0 && a(w+1)<=0);
%end
%u = find(b==1);
%index_infl = u(1);
%figure(4)
%plot(T,Yi,'--','Color',[0.5,0.5,0.5],'linewidth',1);

u = find(Yir==(round(0.8*Km):round(0.8*Km+2)));
index_infl = u(1);

%if r==1
%    b=b(200:700);
%    u = find(b==1);
%index_infl = u(end);
%end

%if r==2
%    b=b(25:80);
%    u = find(b==1);
%index_infl = u(1) + 25;
%end

%if r==3
%    b=b(300:500);
%    u = find(b==1);
%index_infl = u(1) + 300;
%end

%if r==4
%    b=b(200:500);
%    u = find(b==1);
%index_infl = u(1) + 200;
%end

%if r==5
%    b=b(200:400);
%    u = find(b==1);
%index_infl = u(1) + 200;
%end



T_eps(r,i) = index_infl;


%plot the trajectories
%figure(2)
%plot(T,Y1(:,1),'Color',[0.5,0.5,0.5],'linewidth',1);
%hold on
%h1 =plot(T,Yi,'--','Color',[0.5,0.5,0.5],'linewidth',1/i);
%hold on
%legend({'trajectories'},'Location','west', 'FontSize',12)
%hold on


end
end



T_eps = T_eps*0.1; %time intervals are 0.1
T_eps_p = (T_eps - T_eps(:,1) ) ./ T_eps(:,1)*100*(-1);
T_eps_time = -T_eps +T_eps(:,1);

s= summer(size(T_eps,1)+1);
figure(1);
subplot(2,2,4)
for z = 1:size(T_eps_time,1)
plot(epsi_v,(T_eps_p(z,:)),'-','Color',s(z,:),'linewidth',2);
hold on
end
lgd=legend('r_m=1.1','r_m=1.5','r_m=1.9','r_m = 2.3', 'location', 'northwest','Fontsize',10);
title('(c)')
xlabel({'Microbiome transfer rate among' ;'introduced individuals (\lambda_m)'})
ylabel({'% reduction in invasion lag time'})
set(gca,'fontsize',12)




