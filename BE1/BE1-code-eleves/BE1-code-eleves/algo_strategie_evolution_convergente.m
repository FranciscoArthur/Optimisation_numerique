%*******************************************************************************
%Cette fonction algo_strategie_evolution_convergente implemente une strategie
% d'evolution globalement convergente
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
%            un_beta              % facteur de contraction
%            un_nit_max           % nombre maximum d'iterations autorisees     %
%            un_f_count_max       % nombre maximum d'evaluations de une_f auto %
%                                 % risees                                     %
%            une_tol_x            % seuil de stationnarite des x_k             %
%
%                                                    %**************************
%                                                    % PARAMETRES EN SORTIE    *
%                                                    %**************************
%
%            x_opt                % la solution proposee par l'algorithme      %
%            f_opt                % une_f (x_opt, varargin{:})                 %
%            fin                  % la cause de l'arret de l'algorithme        %
%            nit                  % le nombre iterations                       %
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
    = algo_strategie_evolution_convergente(                            ...
    une_f          ,                                ...
    un_x0          ,                                ...
    un_sigma_0     ,                                ...
    un_beta        ,                                ...
    un_nit_max     ,                                ...
    un_f_count_max ,                                ...
    une_tol_x                                       ...
    )

%**********************
%    CODE             *
%**********************

global f_count  ;                  % nombre     d'evaluations de la
% fonction a minimiser, sans pitie

f_count    = 0                                                             ;

% --------------------  Initialization --------------------------------
% User defined input parameters (need to be edited)
xmean      = un_x0'                                                        ;
fmean      = feval(une_f,xmean)                                            ;
sigma_es   = un_sigma_0                                                    ;
N          = length(xmean   )                                              ;
% Strategy parameter setting: Selection
lambda = 4+floor(3*log(N));  % population size, offspring number
mu = lambda/2;               % number of parents/points for recombination
weights = log(mu+1/2)-log(1:mu)'; % muXone array for weighted recombination
mu = floor(mu);
weights = weights/sum(weights);     % normalize recombination weights array
mueff=sum(weights)^2/sum(weights.^2); % variance-effectiveness of sum w_i x_i

% Strategy parameter setting: Adaptation
cc = (4 + mueff/N) / (N+4 + 2*mueff/N); % time constant for cumulation for C
cs = (mueff+2) / (N+mueff+5);  % t-const for cumulation for sigma control
c1 = 2 / ((N+1.3)^2+mueff);    % learning rate for rank-one update of C
cmu = min(1-c1, 2 * (mueff-2+1/mueff) / ((N+2)^2+mueff));  % and for rank-mu update
damps = 1 + 2*max(0, sqrt((mueff-1)/(N+1))-1) + cs; % damping for sigma
% usually close to 1
% Initialize dynamic (internal) strategy parameters and constants
pc = zeros(N,1); ps = zeros(N,1);   % evolution paths for C and sigma
B = eye(N,N);                       % B defines the coordinate system
D = ones(N,1);                      % diagonal D defines the scaling
C = B * diag(D.^2) * B';            % covariance matrix C
invsqrtC = B * diag(D.^-1) * B';    % C^-1/2
eigeneval = 0;                      % track update of B and D
chiN=N^0.5*(1-1/(4*N)+1/(21*N^2));  % expectation of
%   ||N(0,I)|| == norm(randn(N,1))

k          =   0                                                           ;
fin        =   0                                                           ;
% -------------------- Generation Loop --------------------------------

while(fin==0 && k < un_nit_max)
    %
    % A Completer
    %
    if(f_count >= un_f_count_max)
        fin =  4                                                           ;
    end
    if(k==un_nit_max)
        fin =  3                                                           ;
    end
    k = k +1 ;
end
x_opt    =   xmean                                                         ;
f_opt    =    feval(une_f,xmean)                                           ;
nit      =   k                                                             ;