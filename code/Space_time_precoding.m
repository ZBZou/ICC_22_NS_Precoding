%% Spatial-time precoding
clear all;
close all;

Nu = 10; % users
T = 100; % time length
taps_max = 20;
taps_min = 10;
NT = 1; 
% data = zeros(Nu, T-taps_max); %data

Taps = randi([taps_min,taps_max],Nu,T);

data = -3 + 2*randi([0,3],Nu,T); %data serie

SNRdB = 10;

sigma2 = NT*0.5*10^(-SNRdB/10); sigma = sqrt(sigma2);

%% NS channel generate
Hst = zeros(Nu,T, Nu, T);
stat = zeros(Nu*T,T);
for u1 = 1:Nu
   for t1 = 1:T
       for u2 = 1:Nu
          for t2 = 1:T
             taps = Taps(u2,t1);
             if t2 > t1 || t2 < t1 - taps
                continue
             end
              vi = sqrt(t2/T)*rand();
              ei = sqrt(t2/T)*rand();
              Hst(u1,t1,u2,t2) = ei + vi*randn();
          end
       end 
       vj = sqrt(t1/T)*randn();
       ej = sqrt(t1/T)*randn();
   end 
end 

% Hst = Hst + 0.01*randn(Nu,T, Nu, T);


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
print -opengl -dpdf -r600 Hst_1_10_50_100.pdf

figure
plt_Ht = squeeze(Hst(1,:,1,:));
mesh((abs(plt_Ht')));
xlabel('t');
xlabh = get(gca,'XLabel');
ylabel('t^\prime');
ylabh = get(gca,'YLabel');
set(gca,'fontsize',30)
[h, wd, ht] = tightfig();
print -opengl -dpdf -r600 Ht_1.pdf

figure
plt_Hs = squeeze(Hst(:,1,:,1));
mesh((abs(plt_Hs)));
xlabel('u^\prime');
xlabh = get(gca,'XLabel');
ylabel('u');
ylabh = get(gca,'YLabel');
set(gca,'fontsize',30)

[h, wd, ht] = tightfig();
print -opengl -dpdf -r600 Hs_1.pdf

%% Boxplot
stat = zeros(Nu, T);
for i = 1:Nu
   for j = 1:T
      stat(i,j) = Hst(i,j,i,j);
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
print -opengl -dpdf -r600 boxplot.pdf

%% GMT

K = reshape(Hst, T*Nu, T*Nu);
[U,temp_sig,V] = svd(K);

DS = diag(temp_sig);
Sig = [];
N = 0;
for i = 1:length(DS)
   if DS(i) < 0.01 
      N = i-1; % N most contributed eigenfunctions
      break;
   end
   N;
   Sig(i) = DS(i);
end 

%% construct X(t)
s = reshape(data,Nu*T,1);
xn = zeros(Nu*T,1);
X = zeros(Nu*T,1);

for i = 1:N
   xn = dot(s, U(:,i))/Sig(i);
   X = X + xn*V(:,i);
end


%%  receive Rx

Tx = reshape(X, Nu, T);
Rx = zeros(Nu, T);
for i = 1:Nu
   for j = 1:T
      Rx(i,j) = sum(sum(squeeze(Hst(i,j,:,:)) .* Tx)) + sigma*randn();
   end
end 

%% plot

plt_X = Tx;
plt_r = Rx;
plt_s = data;

figure;
subplot(2,1,1)
[x1,y1] = meshgrid(1:size(plt_X,1), 1:size(plt_X,2));  % T is your table
plot3(x1,y1,plt_X);
grid on;
xlabel('u');
ylabel('t');
view(-15,50);
title('Transmit signal X(u,t)');
set(gca,'fontsize',20)

subplot(2,1,2)
[x2,y2] = meshgrid(1:size(plt_s,1), 1:size(plt_s,2));  % T is your table
plot3(x2,y2,plt_s,'--','color',[0 0.4470 0.7410]);
hold on;
plot3(x2,y2,plt_r,'.-','color',[0.8500 0.3250 0.0980]);
grid on;
xlabel('u');
ylabel('t');
view(-55,82);
title('Data signal s(u,t) and receive signal r(u,t)');
set(gca,'fontsize',20)

figure
for i = 1:Nu
   subplot(Nu,1,i)
   plot([1:T], plt_s(i,:),'--');
   hold on;
   plot([1:T], plt_r(i,:),'.-');
   xlabel("u_"+"{"+i+"}");
   axis([1 T -5 5])
end

% estimate/division
ne1 = zeros(Nu,1);
ne2 = zeros(Nu,1);

for i = 1:Nu
   hat_r1 = round(plt_r(i,:));
   hat_r2 = round(plt_r(i,taps_max+1:end));
   hat_r1(find(hat_r1 < -3)) = -3;
   hat_r1(find(hat_r1 > 3)) = 3;
   hat_r2(find(hat_r2 < -3)) = -3;
   hat_r2(find(hat_r2 > 3)) = 3;
   ne1(i) = length(find(squeeze(plt_s(i,1:end)) - hat_r1(1:end)));
   ne2(i) = length(find(squeeze(plt_s(i,taps_max+1:end)) - hat_r2(1:end)));
end

BER1 = sum(ne1)/(Nu*T)
BER2 = sum(ne2)/(Nu*(T - taps_max))


%%
phi_1 = reshape(U(:,1),Nu,T);
phi_2 = reshape(U(:,2),Nu,T);
phi_3 = reshape(U(:,3),Nu,T);

psi_1 = reshape(V(:,1),Nu,T);
psi_2 = reshape(V(:,2),Nu,T);
psi_3 = reshape(V(:,3),Nu,T);
figure;
subplot(2,2,1)
mesh(phi_1);
grid on;
xlabel('u^\prime');
ylabel('t^\prime');
view(-15,50);
set(gca,'fontsize',20)
title('\phi_1 (u^\prime,t^\prime)');
subplot(2,2,2)
mesh(phi_2);
grid on;
xlabel('u^\prime');
ylabel('t^\prime');
view(-15,50);
title('\phi_2 (u^\prime,t^\prime)');
set(gca,'fontsize',20)
subplot(2,2,3)
mesh(psi_1);
grid on;
xlabel('u');
ylabel('t');
view(-15,50);
set(gca,'fontsize',20)
title('\psi_1 (u,t)');
subplot(2,2,4)
mesh(psi_2);
grid on;
xlabel('u');
ylabel('t');
title('\psi_2 (u,t)');
view(-15,50);
set(gca,'fontsize',20)

print -opengl -dpdf -r600 Eigenfunctions.pdf