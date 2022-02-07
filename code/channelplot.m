%% channel plot
clear all;
close all;

Nu = 30; % users
T = 100; % time length
taps_max = 20;
taps_min = 10;
NT = 1; % 1 antenna
Nloop = 200;
M = 16;
K = log2(M);
Taps = randi([taps_min,taps_max],Nu,T);

data = randi([0 1],Nu,T*K);
s = qammod(data',M,'InputType','bit','UnitAveragePower',true);
s = s';

mode = 1; % 0:time-space 1:space

ep = [0.1, 0.05, 0.01];
Nep = length(ep);

SNRdB = [0:1:20];
Ber = zeros(length(SNRdB),Nep);
Ber_dpc = zeros(length(SNRdB),1);
Ber_ideal = zeros(length(SNRdB),1);

Hst = zeros(Nu,T, Nu, T);
sq2 = square(2);

for u1 = 1:Nu
   for t1 = 1:T
       for u2 = 1:Nu
          for t2 = 1:T
             taps = Taps(u2,t1);
             if t2 > t1 || t2 < t1 - taps
                continue
             end
             vi = sqrt(t2/T)*rand();
             ei = sqrt(t2/T)*(rand()+rand()*1i)/sq2;
             Hst(u1,t1,u2,t2) = ei+vi*(randn()+randn()*1i)/sq2;
          end
       end 
   end 
end 

%% plot Hst
plt_Hst = squeeze(Hst(1,100,:,:));

sub_Hst1 = squeeze(Hst(1,1,:,:));
sub_Hst2 = squeeze(Hst(1,10,:,:));
sub_Hst3 = squeeze(Hst(1,50,:,:));
sub_Hst4 = squeeze(Hst(1,100,:,:));

figure
subplot(2,2,1)
mesh(((abs(sub_Hst1))));
xlabel('t^\prime');
xlabh = get(gca,'XLabel');
ylabel('u^\prime');
ylabh = get(gca,'YLabel');
title('(a)','fontsize',14,'fontweight','b');
set(gca,'fontsize',14)
view(-15,15);
subplot(2,2,2)
mesh(((abs(sub_Hst2))));
xlabel('t^\prime');
xlabh = get(gca,'XLabel');
ylabel('u^\prime');
ylabh = get(gca,'YLabel');
title('(b)','fontsize',14,'fontweight','b');
set(gca,'fontsize',14)
view(-15,15);
subplot(2,2,3)
mesh(((abs(sub_Hst3))));
xlabel('t^\prime');
xlabh = get(gca,'XLabel');
ylabel('u^\prime');
ylabh = get(gca,'YLabel');
title('(c)','fontsize',14,'fontweight','b');
set(gca,'fontsize',14)
view(-15,15);
subplot(2,2,4)
mesh(((abs(sub_Hst4))));
xlabel('t^\prime');
xlabh = get(gca,'XLabel');
ylabel('u^\prime');
ylabh = get(gca,'YLabel');
title('(d)','fontsize',14,'fontweight','b');
% title('Time-varying response of user1','fontsize',14,'fontweight','b');
set(gca,'fontsize',14)
view(-15,15);
[h, wd, ht] = tightfig();
% print -opengl -dpdf -r600 Hst_1_10_50_100.pdf

figure
plt_Ht = squeeze(Hst(1,:,1,:));
mesh((abs(plt_Ht')));
xlabel('t');
xlabh = get(gca,'XLabel');
ylabel('t^\prime');
ylabh = get(gca,'YLabel');
set(gca,'fontsize',30)
[h, wd, ht] = tightfig();
% print -opengl -dpdf -r600 Ht_1.pdf

figure
plt_Hs = squeeze(Hst(:,1,:,1));
mesh((abs(plt_Hs)));
xlabel('u^\prime');
xlabh = get(gca,'XLabel');
ylabel('u');
ylabh = get(gca,'YLabel');
set(gca,'fontsize',30)

[h, wd, ht] = tightfig();
% print -opengl -dpdf -r600 Hs_1.pdf

%% Boxplot
stat = zeros(Nu, T);
for i = 1:Nu
   for j = 1:T
      stat(i,j) = abs(Hst(i,j,i,j));
   end
end

figure
boxplot(stat);
hold on
xlabel('t');
ylabel('Mean channel gain');
xtics=[1:10:T];
xticlab={'0' '10' '20' '30' '40' '50' '60' '70' '80' '90' '100'};
set(gca,'XTick',xtics,'XTickLabel',xticlab,'fontsize',30)

[h, wd, ht] = tightfig();
% print -opengl -dpdf -r600 boxplot.pdf


%% eigenfuntions
h_r = reshape(real(Hst), Nu*T, Nu*T);
h_i = reshape(imag(Hst), Nu*T, Nu*T);

H = [h_r, -h_i;h_i, h_r];

[U1,Sig,V1] = svd(H);

Ur = U1(1:Nu*T,:);
Ui = U1(Nu*T+1:end,:)*1i;
U = Ur+Ui;

Vr = V1(1:Nu*T,:);
Vi = V1(Nu*T+1:end,:)*1i;
V = Vr+Vi;
%%
phi_1 = reshape(U(:,1),Nu,T);
phi_2 = reshape(U(:,3),Nu,T);
phi_3 = reshape(U(:,3),Nu,T);

psi_1 = reshape(V(:,1),Nu,T);
psi_2 = reshape(V(:,3),Nu,T);
psi_3 = reshape(V(:,3),Nu,T);
figure;
subplot(2,2,1)
mesh(abs(phi_1));
grid on;
xlabel('u^\prime');
ylabel('t^\prime');
view(-15,50);
set(gca,'fontsize',20)
title('\phi_1 (u^\prime,t^\prime)');
subplot(2,2,2)
mesh(abs(phi_2));
grid on;
xlabel('u^\prime');
ylabel('t^\prime');
view(-15,50);
title('\phi_2 (u^\prime,t^\prime)');
set(gca,'fontsize',20)
subplot(2,2,3)
mesh(abs(psi_1));
grid on;
xlabel('u');
ylabel('t');
view(-15,50);
set(gca,'fontsize',20)
title('\psi_1 (u,t)');
subplot(2,2,4)
mesh(abs(psi_2));
grid on;
xlabel('u');
ylabel('t');
title('\psi_2 (u,t)');
view(-15,50);
set(gca,'fontsize',20)
% print -opengl -dpdf -r600 Eigenfunctions.pdf

%%
figure;
subplot(2,2,1)
mesh(real(phi_1));
grid on;
xlabel('u^\prime');
ylabel('t^\prime');
view(-15,50);
set(gca,'fontsize',20)
title('R_\phi_1 (u^\prime,t^\prime)');
subplot(2,2,2)
mesh(imag(phi_1));
grid on;
xlabel('u^\prime');
ylabel('t^\prime');
view(-15,50);
title('I_\phi_1 (u^\prime,t^\prime)');
set(gca,'fontsize',20)
subplot(2,2,3)
mesh(real(phi_2));
grid on;
xlabel('u');
ylabel('t');
view(-15,50);
set(gca,'fontsize',20)
title('R_\phi_2 (u,t)');
subplot(2,2,4)
mesh(imag(phi_2));
grid on;
xlabel('u');
ylabel('t');
title('I_\phi_2 (u,t)');
view(-15,50);
set(gca,'fontsize',20)

% print -opengl -dpdf -r600 Eigenfunctions2.pdf