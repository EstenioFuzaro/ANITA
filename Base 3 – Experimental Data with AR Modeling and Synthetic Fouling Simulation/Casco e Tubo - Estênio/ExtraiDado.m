clc; close all; clear t;
t(1) = 0;
for i = 2:length(Temp_s_casco)
    t(i,:) = i-1;
end
Temp_s_casco = Temp_s_casco + 6;
 save CascoTubo_E6 Temp_s_casco t