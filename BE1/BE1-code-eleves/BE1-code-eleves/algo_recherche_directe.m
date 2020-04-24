%*******************************************************************************
%Cette fonction algo_recherche_directe implemente l'algorithme de recherche
%directe directionnelle
%                                                                              *
%                                                                              *
%                                                                              *
%*******************************************************************************
%                                                    %**************************
%                                                    % PARAMETRES EN ENTREE    *
%                                                    %**************************
%
%            une_f                % une fonction dont on cherche un minimum    %
%
%            un_x0                % un incontournable point initial            %
%            un_alpha_0           % un pas initial
%            D                    % un ensemble de directions (si possible PSS)
%            un_gamma             % facteur d'expansion 
%            un_beta              % facteur de contraction
%            un_nit_max           % nombre maximum d'iterations autorisees     %
%            un_f_count_max       % nombre maximum d'evaluations de une_f auto %
%                                 % risees                                     %
%            une_tol_x            % seuil de stationnarite des x_k             %
%            une_tol_f            % seuil de stationnarite des f_k             %
%
%                                                    %**************************
%                                                    % PARAMETRES EN SORTIE    *
%                                                    %**************************
%
%            x_opt                % la solution proposee par l'algorithme      %
%            f_opt                % une_f (x_opt, varargin{:})                 %
%            fin                  % la cause de l'arret de l'algorithme        %
%            nit                  % le nombre iterations de l'algorithme       %
%            f_count              % le nombre d'evaluations de une_f           %
%
% Responsable: Y. Diouane (youssef.diouane@isae.fr) -- 2017/2018
% (C) Institut Supérieur de l'Aéronautique et de l'Espace (ISAE-Supaéro)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [                                                             ...
    x_opt  ,                                                           ...
    f_opt  ,                                                           ...
    fin    ,                                                           ...
    nit                                                                ...
    ]                                                                  ...
    = algo_recherche_directe(                                          ...
    une_f          ,                                ...
    un_x0          ,                                ...
    un_alpha_0     ,                                ...
    D              ,                                ...
    un_gamma       ,                                ...
    un_beta        ,                                ...
    un_c           ,                                ...
    un_nit_max     ,                                ...
    un_f_count_max ,                                ...
    une_tol_x      ,                                ...
    une_tol_f                                       ...
    )

%**********************
%    CODE             *
%**********************

global f_count  ;                  % nombre     d'evaluations de la
% fonction a minimiser, sans pitie

f_count    = 0                                                             ;


% tempx sera desormais    %
% vecteur colonne         %
tempx      = un_x0'                                                        ;
alpha      = un_alpha_0                                                    ;
fdex       = feval(une_f,tempx)                                            ;
n          = length(tempx   )                                              ;
r          = size(D,2)                                                     ;
k          =   0                                                           ;
fin        =   0                                                           ;

f_test = @(xk,alphak,d) feval(une_f,xk+alphak*d)-(feval(une_f,xk));

while(fin==0 && k < un_nit_max)
    isSuccess = false;
    
    for search = 1:1:r
        if (f_test(tempx,alpha,D(:,search)) < -(10^-2)/2*(D(1,search)^2+D(2,search)^2)*alpha^2)
            isSuccess = true;
            d = D(:,search);
        end
    end
    
    if (isSuccess == true) % iteration successful
        tempx_next = tempx+alpha*d;
        alpha = un_gamma*alpha;
    else
        tempx_next = tempx;
        alpha = un_beta*alpha;
    end
   
    fdex       = feval(une_f,tempx);
    fdex_next  = feval(une_f,tempx_next);
    
    if(alpha <= 10^-6) 
        fin =  1                                                           ;
    end
    
    if(isSuccess == true) % erreur
        if (abs((fdex_next - fdex)/(tempx_next - tempx)) <= 10^-6)
            fin =  2                                                           ;
        end
    end
    
    if(f_count >= un_f_count_max)
        fin =  4                                                           ;
    end
    
    if(k==un_nit_max)
        fin =  3                                                           ;
    end
    
    tempx = tempx_next;
    k = k +1 ;
end
x_opt    =   tempx                                                         ;
f_opt    =    fdex                                                         ;
nit      =   k                                                             ;
