function ne_dpc = DPC (Nu, T, s, data, Hst, sigma, M, taps_max,mode)
% Nu = 4;
% T = 100;
% M = 4;
% mode = 1;
% data = ones(Nu, T);
% SNRdB = [0:1:30];
% Ber_dpc = zeros(length(SNRdB),1);
% taps_max = 20;
% taps_min = 10;
% NT = 1; % 1 antenna
% Nloop = 100;
% % data = zeros(Nu, T-taps_max); %data
% 
% Taps = randi([taps_min,taps_max],Nu,T);
% 
% for n = 1:length(SNRdB)
%    sigma2 = NT*0.5*10^(-SNRdB(n)/10); sigma = sqrt(sigma2);
%    for m = 1:Nloop
%       Hst = zeros(Nu,T, Nu, T);
%       for u1 = 1:Nu
%          for t1 = 1:T
%              for u2 = 1:Nu
%                 for t2 = 1:T
%                    taps = Taps(u2,t1);
%                    if t2 > t1 || t2 < t1 - taps
%                       continue
%                    end
%                     vi = sqrt(t2/T)*rand();
%                     ei = sqrt(t2/T)*rand();
%                     Hst(u1,t1,u2,t2) = randn();
%                 end
%              end 
%          end 
%       end 
K = log2(M);
L = zeros(Nu,Nu,T);
Q = zeros(Nu,Nu,T);
for i = 1:T
   H = squeeze(Hst(:,i,:,i));
   [Q_temp,R_temp] = qr(H');
   L(:,:,i)=R_temp'; Q(:,:,i)=Q_temp';
end

xp = s;
Tx_signal = s;
for i = 1:T
   for m=2:Nu % Eqs.(13.39)(13.41)
       xp(m,i) = xp(m,i) - L(m,1:m-1,i)/L(m,m,i)*xp(1:m-1,i);
   end
   Tx_signal(:,i) = Q(:,:,i)'*xp(:,i);
end

Rx_signal = zeros(Nu,T);
r = zeros(Nu,T);

if mode == 0 % time-space
   for i = 1:Nu
      for j = 1:T
         Rx_signal(i,j) = sum(sum(squeeze(Hst(i,j,:,:)) .* Tx_signal)) + sigma*(randn()+randn()*1i);
      end
   end 
elseif mode == 1 % space
   for i = 1:T
       H = squeeze(Hst(:,i,:,i));
       Rx_signal(:,i) = H*Tx_signal(:,i) + sigma*(randn()+randn()*1i);
   end
end

%% equalization
if mode == 0 
   for i = 1:Nu
      for j = 1:T
         H_tau = squeeze(Hst(i,j,i,:));
         H_tau(j) = 0; % make gains at time j be zeros, then others would be delays
         Delays = dot(H_tau,Rx_signal(i,:)); % delay effects for one user
         Rx_signal(i,j) = Rx_signal(i,j) - Delays;
      end
   end 
end

%% estimation
for i = 1:T
   inv(diag(diag(L(:,:,i))));
   r(:,i) = inv(diag(diag(L(:,:,i))))*Rx_signal(:,i);
end

if mode == 1
   cp = 0;
elseif mode == 0
   cp = taps_max;
end
      
 
ne_dpc = myber(data(:,K*cp+1:end),r(:,cp+1:end),M);
%       Ber_dpc(n) = Ber_dpc(n) + ne_dpc/(Nu*log2(M)*T);
%    end
%    Ber_dpc(n,:) = Ber_dpc(n,:)/Nloop;
% end
% 
% %%
% semilogy(SNRdB,Ber_dpc);