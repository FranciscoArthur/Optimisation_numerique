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
while(fin==0 && k < un_nit_max)
    %
    % A completer 
    %
    if(f_count >= un_f_count_max)
        fin =  4                                                           ;
    end
    if(k==un_nit_max)
        fin =  3                                                           ;
    end
    k = k +1 ;
end
x_opt    =   tempx                                                         ;
f_opt    =    fdex                                                         ;
nit      =   k                                                             ;
