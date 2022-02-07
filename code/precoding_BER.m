%% BER
clear all;
close all;

Nu = 10; % users
T = 100; % time length
taps_max = 20;
taps_min = 10;
NT = 1; % 1 antenna
Nloop = 200;

mode = 0; % 0:time-space 1:space

ep = [0.001]; % epsilon

Nep = length(ep);

SNRdB = [0:1:20];
Ber = zeros(length(SNRdB),Nep);
Ber_dpc = zeros(length(SNRdB),1);
Ber_ideal = zeros(length(SNRdB),1);

Ms = [2, 4, 16, 64]; %modulation

Ber2 = zeros(length(SNRdB),Nep, length(Ms));
Ber_dpc2 = zeros(length(SNRdB),1, length(Ms));
Ber_ideal2 = zeros(length(SNRdB),1, length(Ms));

for l = 1:length(Ms)
   M = Ms(l);
   K = log2(M);
   Taps = randi([taps_min,taps_max],Nu,T);
   data = randi([0 1],Nu,T*K);
   s = qammod(data',M,'InputType','bit','UnitAveragePower',true);
   s = s';
   for n = 1:length(SNRdB)
      sigma2 = NT*0.5*10^(-SNRdB(n)/10); sigma = sqrt(sigma2);
      sq2 = square(2);
      for m = 1:Nloop
         % NS channel generate
         l
         n
         m
         Hst = zeros(Nu,T, Nu, T);
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

         % Eigenprecoding
         if mode == 0 %space-time

               ne_dpc = DPC(Nu, T, s, data, Hst, sigma, M, taps_max,mode);

               Ber_dpc(n) = Ber_dpc(n) + ne_dpc/(Nu*log2(M)*T);
               % EP
               h_r = reshape(real(Hst), Nu*T, Nu*T);
               h_i = reshape(imag(Hst), Nu*T, Nu*T);

               H = [h_r, -h_i;h_i, h_r];

               [U,Sig,V] = svd(H);
               Sig = diag(Sig);

               N = zeros(Nep,1);
               for i = 1:Nep
                  for j = 1:length(Sig)                
                     if Sig(j) >= ep(i)
                        N(i) = j; % N most contributed eigenfunctions
                     end
                  end 
               end

               % construct X(t)
               sr = reshape(real(s),[],1);
               si = reshape(imag(s),[],1);

               S = [sr;si];

               X = zeros(2*Nu*T,Nep);

               for i = 1:Nep
                  for j = 1:N(i)
                     xn = dot(S, U(:,j))/Sig(j);
                     X(:,i) = X(:,i) + xn*V(:,j);
                  end
               end
               Xr = X(1:Nu*T,:);
               Xi = X(Nu*T+1:end,:)*1i;
               X = Xr+Xi;
               %receive Rx

               Rx = zeros(Nu, T, Nep);
               for i = 1:Nu
                  for j = 1:T
                     for q = 1:Nep
                        X_temp = squeeze(X(:,q));
                        Tx = reshape(X_temp, Nu,T);
                        Rx(i,j,q) = sum(sum(squeeze(Hst(i,j,:,:)) .* Tx)) + sigma*(randn()+randn()*1i)/sq2;
                     end
                  end
               end 
         elseif mode == 1 %space
               % DPC
               ne_dpc = DPC(Nu, T, s, data, Hst, sigma, M, 0, mode);
              
               Ber_dpc(n) = Ber_dpc(n) + ne_dpc/(Nu*log2(M)*T);
               % EP

               Rx = zeros(Nu, T, Nep);
               for k = 1:T
                  Hu = squeeze(Hst(:,k,:,k));
                  h_r = reshape(real(Hu), Nu, Nu);
                  h_i = reshape(imag(Hu), Nu, Nu);
                  H = [h_r, -h_i;h_i, h_r];
                  [U,temp_sig,V] = svd(H);
                  Sig = diag(temp_sig);
                  N = zeros(Nep,1);
                  for i = 1:Nep
                     for j = 1:length(Sig)                
                        if Sig(j) >= ep(i)
                           N(i) = j; % N most contributed eigenfunctions
                        end
                     end 
                  end
                  % construct X(t)
                  sr = reshape(real(s(:,k)),[],1);
                  si = reshape(imag(s(:,k)),[],1);

                  S = [sr;si];

                  X = zeros(2*Nu,Nep);

                  for i = 1:Nep
                     for j = 1:N(i)
                        xn = dot(S, U(:,j))/Sig(j);
                        X(:,i) = X(:,i) + xn*V(:,j);
                     end
                  end
                  Xr = X(1:Nu,:);
                  Xi = X(Nu+1:end,:)*1i;
                  X = Xr+Xi;

                  % receive signal
                  for i = 1:Nep
                      Rx(:,k,i) = Hu*X(:,i) + sigma*(randn()+randn()*1i)/sq2;
                  end
               end
            end
         % ideal
         Rx_ideal = s + sigma*(randn(Nu,T)+randn(Nu,T)*1i)/sq2;

         % estimate/division
         hat_s = zeros(Nu, K*T, Nep);
         for i = 1:Nep
            temp_Rx = squeeze(Rx(:,:,i));
            temp_s = qamdemod(temp_Rx',M,'OutputType','bit','UnitAveragePower',true);
            hat_s(:,:,i) = temp_s';
         end

         %ber
         if mode == 1
            cp = 0;
         elseif mode == 0
            cp = taps_max;
         end

         for i = 1:Nep
             ne = myber(data(:,K*cp+1:end),Rx(:,cp+1:end,i),M);
             Ber(n,i) =  Ber(n,i) + ne/(Nu*log2(M)*(T-cp));
         end

         ne_ideal = myber(data(:,K*cp+1:end),Rx_ideal(:,cp+1:end),M);
         Ber_ideal(n) =  Ber_ideal(n) + ne_ideal/(Nu*log2(M)*(T-cp));

      end
      Ber(n,:) = Ber(n,:)/Nloop;
      Ber_dpc(n,:) = Ber_dpc(n,:)/Nloop;
      Ber_ideal(n,:) = Ber_ideal(n,:)/Nloop;
   end
   Ber2(:,:,l) = Ber;
   Ber_dpc2(:,:,l) = Ber_dpc;
   Ber_ideal2(:,:,l) = Ber_ideal;
end

%% save
if mode == 0
   save('BER_space_Time.mat','Ber2');
else
   save('BER_space.mat','Ber2');
end

%%

% semilogy(SNRdB,Ber_dpc2(:,1,3),'-o',SNRdB,Ber2(:,1,3),'-o',SNRdB,Ber2(:,2,3),'-o',...
%    SNRdB,Ber2(:,3,3),'-o',SNRdB,Ber_ideal2(:,1,3),'-o','linewidth',2), grid on

% legend('DPC with equalization','HOGMT-Precoding: \epsilon {=} 10^{-2}','HOGMT-Precoding: \epsilon {=} 10^{-3}',...
%    'HOGMT-Precoding: \epsilon {=} 10^{-4}','Ideal(AWGN)','Location','Southwest','fontsize',16)

%%
figure;
semilogy(SNRdB,Ber2(:,1,1),'-o',SNRdB,Ber2(:,1,2),'-o',SNRdB,Ber2(:,1,3),'-o',...
   SNRdB,Ber2(:,1,4),'-o','linewidth',2), grid on


legend('HOGMT-Precoding:BPSK','HOGMT-Precoding:QPSK',...
   'HOGMT-Precoding:16QAM','HOGMT-Precoding:64QAM','Location','Southwest','fontsize',16)
xlabel ('SNR (dB)')
ylabel('BER (dB)')
set(gca, 'fontsize', 20)
[h, wd, ht] = tightfig();
if mode == 0
   print -opengl -dpdf -r600 Ber_space_time_temp.pdf
elseif mode == 1
   print -opengl -dpdf -r600 Ber_space_temp.pdf   
end
%% plot H_s
% Hs = squeeze(Hst(:,1,:,1));
% mesh((abs(Hs)));
% xlabel('u^\prime');
% xlabh = get(gca,'XLabel');
% ylabel('u');
% ylabh = get(gca,'YLabel');
% set(gca,'fontsize',30)
% 
% [h, wd, ht] = tightfig();
% print -opengl -dpdf -r600 Hs.pdf
