function ne = myber(data,r,M)
hat_s = qamdemod(r',M,'OutputType','bit','UnitAveragePower',true);
hat_s = hat_s';

data_bits = reshape(data,1,[]);
r_bits = reshape(hat_s,1,[]);
ne = length(find(data_bits - r_bits));

