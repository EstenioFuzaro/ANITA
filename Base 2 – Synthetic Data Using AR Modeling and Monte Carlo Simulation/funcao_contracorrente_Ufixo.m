function dTdt = funcao_contracorrente_Ufixo(t,T,a1,b1,a2,b2,T_eq,T_ef,U)
% Sistema de EDO para trocador de calor considerando decaimento do U
dT = zeros(3,1);
% =================== Parâmetros do Modelo  ====================

Tt0=T_eq;
Tc0=T_ef;

dT(1) = U*a1*(T(1) - T(2)) + b1*(Tt0 - T(1));  %Temperatura do fluido quente

dT(2) = U*a2*(T(1) - T(2)) + b2*(Tc0 - T(2));  %Temperatura do fluido frio

%dT(3) = odefun;                                        %Valores de U
 % =======.=========================================================
 dTdt=[dT(1);dT(2)];
end