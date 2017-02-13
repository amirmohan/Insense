function [z] = Insense(D,m,param)
% This linear program solves:
%
%     min Average(mu^2(z))
%     s.t.
%          sum(z) = m
%          0 =< z < 1
% Given a dictionary D, choose the best m rows
% so that, the resulting sub-dictionary has a low pairwise coherence.
%
% OUTPUT:
% z = vector whose nonzero elements are indices of the optimal selection
%
% INPUTS
% D: Dictionary
% m: Row-sparsity inducing parameter
% param.tol
% param.maxit
% param.gamma: step size


tol = param.tol;
maxit= param.maxit;


[d,n] = size(D);

% usefull matrices for gradient calculation
AA= zeros(n^2,d);
BB= zeros(n,d);
for i = 1:d
    tmp = D(i,:)'*D(i,:);
    BB(:,i) = diag(tmp);
    %tmp(logical(eye(n)))=0;
    AA(:,i)= tmp(:);
end


%----initialization
z = param.init;
OBJ_prev = AVG_mu_B(z,D);

iter=1;
rel_OBJ = tol+10;


%----Projected Descend MAIN loop------
while (iter <= maxit)&&( rel_OBJ >= tol)
    PGRAD = GRAD(z,D,AA); % projected Grad along the Simplex face
    
    % set stepsize (with back tracking or manual)
    if (param.backtrack == 1)
        gamma = BackTrack(z,m,PGRAD,D);
    else
        gamma = param.gamma;
    end
    
    % Gradient descend
    z = z - gamma * PGRAD;
    
    % SBS projection
    z = proxm_fast(z,m); 
    OBJ = AVG_mu_B(z,D);
    rel_OBJ = abs(OBJ-OBJ_prev)/abs(OBJ);
    OBJ_prev = OBJ;
    iter = iter+1;
    
end

%---- Average Coherence 
function OBJ = AVG_mu_B(z,D)
EPS = 1e-10;
OBJ = D'*diag(z)*D;
OBJ = (OBJ.^2+EPS)./( diag(OBJ)*(diag(OBJ))' +EPS^2);
OBJ(logical(eye( size(OBJ,1) ))) = 0; % set diag zero
OBJ = sum(sum(OBJ));

%---Gradient
function [g] = GRAD(z,D,AA)
EPS = 1e-10;
G = D'*diag(z)*D;
g = G;
g = 2*g./( diag(g)*(diag(g))'+ EPS^2);
tmp = (ones(size(G,1),1)*diag(G)') .* (G.^2 +EPS)./( diag(G)*(diag(G))' +EPS^2).^2 ;
tmp(logical(eye( size(tmp,1) ))) = 0; % set diag zero
tmp = sum(tmp,2);
g(logical(eye( numel(tmp) ))) = -2*tmp;

g = AA'*g(:);

function gamma = BackTrack(z,m,Grad,D)
gamma = 1e2;
Energy0 = AVG_mu_B(z,D);
Diff = AVG_mu_B(proxm_fast(z - gamma*Grad,m),D) - Energy0;
iter2 = 1;
while (Diff > -gamma/50*norm(Grad)^2) && (iter2 < 30) 
     
    gamma = gamma/2;
    Energy = AVG_mu_B(proxm_fast(z - gamma*Grad,m),D);
    Diff = Energy - Energy0;
    iter2 = iter2+1;
end

