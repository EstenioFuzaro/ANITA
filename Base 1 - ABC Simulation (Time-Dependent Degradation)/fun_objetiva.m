function [J] = fun_objetiva(Parameters)
%
%  fun_objetiva
%    This function evaluates the error.
%
%  USAGE: [J] = myLikelihood(Parameters)
%__________________________________________________________________________
%  OUTPUTS
% 
%    J : Função objetiva - erro.
%__________________________________________________________________________
%  INPUTS
%    
%    Parameters : Computational model static configuration parameters:
%     - Parameters.varExp : experimental variance 
%     - Parameters.f : experimental force
%     - Parameters.Uexp : experimental displacement
%     - Parameters.DIM : dimension of the problem
%__________________________________________________________________________


%%

% tmp = [];
idx = randperm(length(Parameters.Texp_tubo));


% for ii = 1:length(Parameters.Texp_tubo)
%  tmp(ii) = (norm(Parameters.Texp_tubo(ii) - Parameters.Ts_tubo)./norm(Parameters.Texp_tubo(ii))) + ...
%      +(norm(Parameters.Texp_casco(ii) - Parameters.Ts_casco)./norm(Parameters.Texp_casco(ii)));
% end
% J = sum(tmp)/length(Parameters.Texp_tubo);
ii = idx(1);
J = (norm(Parameters.Texp_tubo(ii) - Parameters.Ts_tubo)./norm(Parameters.Texp_tubo(ii))) + ...
   +(norm(Parameters.Texp_casco(ii) - Parameters.Ts_casco)./norm(Parameters.Texp_casco(ii)));

J = J./2;
%varExp= (Parameters.varExp_tubo+Parameters.varExp_casco)/k;

end