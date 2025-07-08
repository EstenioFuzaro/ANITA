clc; close all; clear t;
t(1) = 0;
for i = 2:length(Temp_s_fria)
    t(i,:) = i-1;
end
 save PlacaPlana_E6 Temp_s_fria t