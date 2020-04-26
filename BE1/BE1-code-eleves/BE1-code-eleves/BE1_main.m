function BE4_main(algo,prob)
addpath('./Problems');
%**************************************************************************
%  Cette fonction constitue le noyau du programme principal relatif a
%  l'implementation de ce BE
% des :
% 1-  algorithme de recherche directe directionnelle
% 2-  une strategie d'evolution
% 3-  une strategie d'evolution globalement convergente
%
%
%    algo: 1, 2, 3 ou 4 pour choisir l'algorithm.
%       algo=1 algorithme de recherche directe directionnelle
%       algo=2 une strategie d'evolution
%       algo=3 une strategie d'evolution globalement convergente
%    prob: 1, 2 ou 3.
%       prob=1 choix de la fonction f1 "rosenbrock"
%       prob=2 choix de la fonction f2 "rastrigin"
%       prob=3 choix de la fonction f3 "ackley"
%
% Responsable: Y. Diouane (youssef.diouane@isae.fr) -- 2017/2018
% (C) Institut Supérieur de l'Aéronautique et de l'Espace (ISAE-Supaéro)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%**************************
% GLOBALES EN MISE A JOUR *
%**************************

global f_count  ;                  % nombre     d'evaluations de la
% fonction a minimiser, sans pitie
global g_count  ;                  % nombre     d'evaluations du
% gradient g_f        , sans pitie
global h_count  ;                  % nombre     d'evaluations de la
% Hessienne h_f       , sans pitie

%***********************************************
%      Pour l'affichage (à ne pas modifier)     *
%***********************************************

ligne_tiret = ['---------------------------------------------------' ]   ;
ligne_tiret = ['|' ligne_tiret  ligne_tiret  '|'];

lentete     = ['| DEPART  |          METHODE         |   PROBLEME     |  FIN  |   F_COUNT   |   NITER    |    F_OPT    |'   ];

les_formats=['|%s| %s |%s |  %3.0f  |   %6.0f    |    %5.0f   |%11.6g  | '                      ];

nom_algo   ={ '    Recherche Directe   ','  Strategie Evolution   ' ,...
    ' Stra. Evol. Convergente'};
nom_prob   ={'   Rosenbrock  ', '   Rastrigin   ' ,'    Ackley     '};
nom_point = {' (1,1.2) ',' (2,2.3) ','(-2.2,1) ',' (2.33,5)'};
%&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&%
%      ETAPE DE RESOLUTION DU PROBLEME D'OPTIMISATION                     %
%&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&%
% en opti non-lineaire il faut un
% point de recherche initial ! x0
x(1,:)=[1 1.2];
x(2,:)=[2 2.3];
x(3,:)=[-2.2 1];
x(4,:)=[2.33 5];
%
% Choix de probleme de minimisation (f1, f2 ou f3)
%
if(prob==1)
    f_fun=@f_rosenbrock;
elseif(prob==2)
    f_fun=@f_rastrigin;
else
    f_fun=@f_ackley;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     Algorithme de recherche directe directionnelle
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if(algo==1)
    %
    % Choix de la base positivement génératrice
    %
    D=[eye(2) -eye(2)]; % la base maximale
    %D=[eye(2) -ones(2,1)]; % la base minimale
    %Q=rand(2); D=[Q -Q]; % la base aleatoire
    %
    for i=1:4 % pour chaque point de part
        top_on(i)           =cputime;
        [x_opt(:,i),f_opt(i),conv(i),ite(i)]= ...
            algo_recherche_directe( ...
            f_fun  ,...
            x(i,:)          ,...
            1,...
            D,...
            2,...
            0.5,...
            0.5,...
            100000,...
            100000,...
            10^-6,...
            10^-6 ...
            );
        % On garde les sorties pour l'afichage ..
        temps(i)                        =cputime - top_on(i)              ;
        x_optimale(:,i,algo)            =x_opt(:,i)                       ;
        t_cpu_time(algo ,i)             = temps(i)                        ;
        t_fin     (algo,i)              =conv(i)                          ;
        t_f_count (algo,i)              =f_count                          ;
        t_nit     (algo,i)              =ite(i)                           ;
        t_f_opt   (algo,i)              =f_opt(i)                         ;
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     Une strategie d'evolution
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if(algo==2)
    for i=1:4 % pour chaque point de part
        top_on(i) =cputime;
        [x_opt(:,i),f_opt(i),conv(i),ite(i)]= ...
            algo_strategie_evolution( ...
            f_fun  ,...
            x(i,:)          ,...
            1,...
            1000,...
            1000,...
            10^-6 ...
            );
        % On garde les sorties pour l'afichage ..
        temps(i)                        =cputime - top_on(i)              ;
        x_optimale(:,i,algo)            =x_opt(:,i)                       ;
        t_cpu_time(algo ,i)             = temps(i)                        ;
        t_fin     (algo,i)              =conv(i)                          ;
        t_f_count (algo,i)              =f_count                          ;
        t_nit     (algo,i)              =ite(i)                           ;
        t_f_opt   (algo,i)              =f_opt(i)                         ;
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     Une strategie d'evolution globalement convergente
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if(algo==3)
    for i=1:4 % pour chaque point de part
        top_on(i) =cputime;
        [x_opt(:,i),f_opt(i),conv(i),ite(i)]= ...
            algo_strategie_evolution_convergente( ...
            f_fun  ,...
            x(i,:)          ,...
            1,...
            0.95, ...
            1000,...
            3000,...
            10^-6 ...
            );
        % On garde les sorties pour l'afichage ..
        temps(i)                        =cputime - top_on(i)              ;
        x_optimale(:,i,algo)            =x_opt(:,i)                       ;
        t_cpu_time(algo ,i)             = temps(i)                        ;
        t_fin     (algo,i)              =conv(i)                          ;
        t_f_count (algo,i)              =f_count                          ;
        t_nit     (algo,i)              =ite(i)                           ;
        t_f_opt   (algo,i)              =f_opt(i)                         ;
    end
end
%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++%
% ON AFFICHE TOUS LES RESULTATS  DE TOUS LES TESTS SOUS FORME D'UN TABLEAU%
%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++%
disp (sprintf     ('\n'                                                 ));
disp (ligne_tiret                                                        );
disp (lentete                                                            );
disp (ligne_tiret                                                        );

for i=1:4 % pour chaque point de départ
    x_optimale(:,i,algo);
    disp (sprintf     ( les_formats  ,                             ...
        nom_point    {i}    ,                                      ...
        nom_algo   {algo}  ,                                       ...
        nom_prob   {prob}   ,                                      ...
        t_fin      (algo,i) ,                                      ...
        t_f_count  (algo,i) ,                                      ...
        t_nit      (algo,i) ,                                      ...
        t_f_opt    (algo,i)  )            )                                  ;
end
disp (ligne_tiret                                                        );


clear           global      %on detruit les variables globales en sortant
