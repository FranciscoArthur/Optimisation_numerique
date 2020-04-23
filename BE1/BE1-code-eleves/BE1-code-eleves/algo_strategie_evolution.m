%*******************************************************************************
%Cette fonction algo_strategie_evolution_convergente implemente une strategie
% d'evolution (CMA-ES -- N. Hansen, 2001)
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
%            un_nit_max           % nombre maximum d'iterations autorisees     %
%            un_f_count_max       % nombre maximum d'evaluations de une_f auto %
%                                 % risees                                     %
%            une_tol_x            % seuil de stationnarite des x_k             %
%
%                                                    %**************************
%                                                    % PARAMETRES EN SORTIE    *
%                                                    %**************************
%
%            x_opt                % la solution proposee par l'algoroithme     %
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
    = algo_strategie_evolution(                                        ...
    une_f          ,                                ... 
    un_x0          ,                                ...
    un_sigma_0     ,                                ...
    un_nit_max     ,                                ...
    un_f_count_max ,                                ...
    une_tol_x                                       ...
    )

%**********************
%    CODE             *
%**********************

global f_count  ;                  % nombre     d'evaluations de la
% fonction a minimiser, sans pitie

f_count    = 0                                                                 ;

% -----------------------  Initialization --------------------------------
% User defined input parameters
xmean      = un_x0'                                                            ;
sigma_es      = un_sigma_0                                                     ;
N          = length(xmean   )                                                  ;
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

k          =   0                                                               ;
fin        =   0                                                               ;
% -------------------- Generation Loop --------------------------------

while(fin==0 && k < un_nit_max)
    % Generate and evaluate lambda offspring
    for k=1:lambda,
        arx(:,k) = xmean + sigma_es * B * (D .* randn(N,1)); % m + sigma^ES * Normal(0,C)
        arfitness(k) = feval(une_f, arx(:,k)); % objective function call
    end
    
    % Sort by fitness and compute weighted mean into xmean
    [arfitness, arindex] = sort(arfitness);  % minimization
    xold = xmean;
    xmean = arx(:,arindex(1:mu)) * weights;  % recombination, new mean value
    
    % Cumulation: Update evolution paths
    ps = (1-cs) * ps ...
        + sqrt(cs*(2-cs)*mueff) * invsqrtC * (xmean-xold) / sigma_es;
    hsig = sum(ps.^2)/(1-(1-cs)^(2*f_count/lambda))/N < 2 + 4/(N+1);
    pc = (1-cc) * pc ...
        + hsig * sqrt(cc*(2-cc)*mueff) * (xmean-xold) / sigma_es;
    
    % Adapt covariance matrix C
    artmp = (1/sigma_es) * (arx(:,arindex(1:mu)) - repmat(xold,1,mu));  % mu difference vectors
    C = (1-c1-cmu) * C ...                   % regard old matrix
        + c1 * (pc * pc' ...                % plus rank one update
        + (1-hsig) * cc*(2-cc) * C) ... % minor correction if hsig==0
        + cmu * artmp * diag(weights) * artmp'; % plus rank mu update
    
    % Adapt step size sigma
    sigma_es = sigma_es * exp((cs/damps)*(norm(ps)/chiN - 1));
    
    % Update B and D from C
    if f_count - eigeneval > lambda/(c1+cmu)/N/10  % to achieve O(N^2)
        eigeneval = f_count;
        C = triu(C) + triu(C,1)'; % enforce symmetry
        [B,D] = eig(C);           % eigen decomposition, B==normalized eigenvectors
        D = sqrt(diag(D));        % D contains standard deviations now
        invsqrtC = B * diag(D.^-1) * B';
    end
    
    % Check stopping criterion
    if(f_count >= un_f_count_max)
        fin =  4                                                           ;
    end
    if(k==un_nit_max)
        fin =  3                                                           ;
    end
    if(sigma_es < une_tol_x  )
        fin =  2                                                           ;
    end
    k = k +1 ;
end
x_opt    =   xmean                                                         ;
f_opt    =    feval(une_f,xmean)                                           ;
nit      =   k                                                             ;