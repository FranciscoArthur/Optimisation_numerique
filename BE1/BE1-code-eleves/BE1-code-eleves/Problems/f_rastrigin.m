%*******************************************************************************
%   Codage de la fonction dite de "rastrigin"                                 *
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
%            fdex                 % f_rastrigin(x)                              %
%
% Responsable: Y. Diouane (youssef.diouane@isae.fr) -- 2017/2018
% (C) Institut Supérieur de l'Aéronautique et de l'Espace (ISAE-Supaéro)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



function fdex = f_rastrigin(x)        

%*****************************
%     LES DIFFERENTS "OBJETS"*
%*****************************
                                                     %**************************
                                                     % GLOBALES EN MISE A JOUR *
                                                     %**************************
          
          global f_count  ;                  % nombre     d'evaluations de 
                                             % fdex=f_rastrigin(x) , sans pitie

                                                     %**************************
                                                     % FONCTIONS MATLAB        *
                                                     %**************************

%         return                             % utiliser help <nom_fonction>
%                      %                     % pour obtenir l'aide en ligne



%*********************
%      CODE          *
%*********************
 
d = length(x);
sum = 0;
for ii = 1:d
	xi = x(ii);
	sum = sum + (xi^2 - 10*cos(2*pi*xi));
end

fdex = 10*d + sum;

f_count          = f_count   + 1                                               ;

return                                                                         ;
