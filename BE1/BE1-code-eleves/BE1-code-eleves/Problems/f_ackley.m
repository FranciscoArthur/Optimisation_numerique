%*******************************************************************************
%   Codage de la fonction dite de "ackley"                                 *
%*****************************************************************************0*
%                                                    %**************************
%                                                    % PARAMETRES EN ENTREE    *
%                                                    %**************************
%
%            x                    % un vecteur de R^2                          %
%
%                                                    %**************************
%                                                    % PARAMETRES EN SORTIE    *
%                                                    %**************************
%
%            fdex                 % f_ackley(x)                              %
%
% Responsable: Y. Diouane (youssef.diouane@isae.fr) -- 2017/2018
% (C) Institut Supérieur de l'Aéronautique et de l'Espace (ISAE-Supaéro)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



function fdex = f_ackley(x)        

%*****************************
%     LES DIFFERENTS "OBJETS"*
%*****************************
                                                     %**************************
                                                     % GLOBALES EN MISE A JOUR *
                                                     %**************************
          
          global f_count  ;                  % nombre     d'evaluations de 
                                             % fdex=f_ackley(x) , sans pitie

                                                     %**************************
                                                     % FONCTIONS MATLAB        *
                                                     %**************************

%         return                             % utiliser help <nom_fonction>
%                      %                     % pour obtenir l'aide en ligne



%*********************
%      CODE          *
%*********************

d = length(x);

c = 2*pi;
b = 0.2;
a = 20;

sum1 = 0;
sum2 = 0;
for ii = 1:d
    xi = x(ii);
    sum1 = sum1 + xi^2;
    sum2 = sum2 + cos(c*xi);
end

term1 = -a * exp(-b*sqrt(sum1/d));
term2 = -exp(sum2/d);

fdex = term1 + term2 + a + exp(1);
f_count          = f_count   + 1                                               ;

return                                                                         ;
