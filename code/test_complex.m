%% test complex space channel
clear all;
close all;

Nu = 10;
T = 100;
Hst = zeros(Nu,T, Nu, T);
M = 16;
K = log2(M);
sq2 = square(2);

for u1 = 1:Nu
   for t1 = 1:T
       for u2 = 1:Nu
          for t2 = 1:T
              vi = sqrt(t2/T)*rand();
              ei = sqrt(t2/T)*(rand()+rand()*1i)/sq2;
              Hst(u1,t1,u2,t2) = ei+vi*(randn()+randn()*1i)/sq2;
          end
       end 
   end 
end 

plt_H = squeeze(Hst(1,:,1,:));
mesh(abs(plt_H));

data = randi([0 1],Nu,T*K);


s = qammod(data',M,'InputType','bit','UnitAveragePower',true);
s = s';

sr = reshape(real(s),[],1);
si = reshape(imag(s),[],1);

S = [sr;si];

h_r = reshape(real(Hst), Nu*T, Nu*T);
h_i = reshape(imag(Hst), Nu*T, Nu*T);

H = [h_r, -h_i;h_i, h_r];

[U,Sig,V] = svd(H);
Sig = diag(Sig);

N = length(Sig);

X = zeros(2*Nu*T,1);

for i = 1:N
   xn = dot(S, U(:,i))/Sig(i);
   X = X + xn*V(:,i);
end

Xr = X(1:Nu*T);
Xi = X(Nu*T+1:end)*1i;

Tx = reshape(Xr+Xi,Nu,T); 

Rx = zeros(Nu, T);
for i = 1:Nu
   for j = 1:T
      H_temp = squeeze(Hst(i,j,:,:));
      Rx(i,j) = sum(sum(H_temp .* Tx)); % no noise
   end
end 

hat_s = qamdemod(Rx',M,'OutputType','bit','UnitAveragePower',true);
hat_s = hat_s';

ber = length(find(round(data - hat_s)))/(Nu*T*K)




