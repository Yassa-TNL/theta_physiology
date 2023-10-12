A = [0 1 1 9 9 9 0 0 0 0]; % clinical
B = [1 1 9 9 9]; % research

flpA = flip(A);
conv_temp = conv(flpA,B, 'full');
figure;plot(conv_temp)
match_lag = find(conv_temp==max(conv_temp));

clinical_region_conv_match = flip(flpA(match_lag-length(B)+1:match_lag));