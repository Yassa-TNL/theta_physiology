% win= round(.135*fs);
% for elec = 1:size(cond1,3)
%     for gamma_freq = find(freq>40)
%        % cond1(gamma_freq,:,elec) = conv(cond1(gamma_freq,:,elec), ones(1,win)/win,'same');
%         cond2(gamma_freq,:,elec) = conv(cond2(gamma_freq,:,elec), ones(1,win)/win,'same');
%         cond3(gamma_freq,:,elec) = conv(cond3(gamma_freq,:,elec), ones(1,win)/win,'same');
%        % cond4(gamma_freq,:,elec) = conv(cond4(gamma_freq,:,elec), ones(1,win)/win,'same');
%     end
% end