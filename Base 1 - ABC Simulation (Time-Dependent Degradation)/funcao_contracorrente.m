function dTdt = funcao_contracorrente(t,T,a1,b1,a2,b2,T_eq,T_ef,odefun)
% Sistema de EDO para trocador de calor considerando decaimento do U
dT = zeros(3,1);
% =================== Parâmetros do Modelo  ====================

Tt0=T_eq;
Tc0=T_ef;

dT(1) = T(3)*a1*(T(1) - T(2)) + b1*(Tt0 - T(1));  %Temperatura do fluido quente

dT(2) = T(3)*a2*(T(1) - T(2)) + b2*(Tc0 - T(2));  %Temperatura do fluido frio

dT(3) = odefun;                                        %Valores de U
 % =======.=========================================================
 dTdt=[dT(1);dT(2);dT(3)];