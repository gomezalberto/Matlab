function [ eVect, eVal, Z, exemplar, stdDev, U, X_norm ] = sPCA( X, L, n_obs, varargin )
%sPCA Implementation of supervised and unsupervised PCA
%   Inputs: X       training data (not normalised) pxn;
%           L       kernel matrix of target variable;
%           n_obs   number of observations
%           varargin:   'normalise' 0 don't normalise input matrix (default)
%                                   1 normalise input matrix
%                       'type'  'normal' (default), 'dual' or 'kernel'
%
% Author: Devis Peressutti, Biomedical Engineering, KCL, 2013

n_optargin = size(varargin,2);

normalise = false;
type = 'normal';

if (~n_optargin)
    s = 'Normal supervised PCA performed on not-normalised input data';
    disp(s);
else
    if n_optargin==2
        if strcmp(varargin{1},'normalise')
            normalise = varargin{2};
        elseif strcmp(varargin{1},'type')
            type = varargin{2};
        end
    elseif n_optargin==4
        if strcmp(varargin{1},'normalise')
            normalise = varargin{2};
            type = varargin{4};
        else
            normalise = varargin{4};
            type = varargin{2};
        end
    else
        error('Wrong number of arguments')
    end
end

% The input matrix is reshaped into a DxN matrix
X = reshape(X,numel(X)/n_obs,n_obs);

% Definition of H (centring matrix)
H = eye(n_obs) - (1/n_obs)*ones(n_obs);

% Mean vector returned
exemplar = mean(X,2);

% % X is centred
X_centred = X*H;

% X is normalised
if normalise
    % disp('Input matrix normalised')
    x_std = std(X_centred,0,2);
    X_norm = X_centred./repmat(x_std,1,n_obs);
else
    % disp('Input matrix not normalised')
    X_norm = X_centred;
    x_std = ones(size(X_centred,1),1);
end

stdDev = x_std;

% Supervised PCA
switch type
    case 'normal'
        % disp('Performing normal supervised PCA')
        
        % Computation of Q
        Q = X_norm*L*transpose(X_norm);
        
        % Computation of eigenvectors and eigenvalues
        [V, E] = eig(full(Q));
        
        % extract eigenvectors and sort by eigenvalues
        eigenVectors = V(:,floor((find(E>0.000001)-1) / size(V,1))+1); % non-zero eigenvalues
        eigenValues = E(find(E>0.000001));
        vv=[eigenValues transpose(eigenVectors)];
        vvs = sortrows(vv,-1);
        eigenVectors = transpose(vvs(:,2:end));
        eigenValues = vvs(:,1);
        perc_expl = 100*eigenValues/sum(eigenValues);
        
        % Get the eigenvalues containing XX% of the total variance
        var = 99;
        i = 1;
        s = sum(perc_expl(1:i));
        while ( s<var )
            i = i+1;
            s = sum(perc_expl(1:i));
        end
        d = i;
        
        % Encode training data
        Z = transpose(transpose(eigenVectors(:,1:d))*X_norm);
        
        eVal = eigenValues;
        eVect = eigenVectors(:,1:d);
        
    case 'dual'
        
        % disp('Performing dual supervised PCA')
        
        % Decomposition of L (LDL decomposition used for semi-positive matrices)
        [l_m,d_m] = ldl(L);
        Del = abs(l_m*sqrt(d_m));
        Psi = X_norm*transpose(Del); 
        
        % Computation of eigenvectors and eigenvalues
        [V, E] = eig(full(transpose(Psi)*Psi));
        
        % extract eigenvectors and sort by eigenvalues
        eigenVectors = V(:,floor((find(E>0.000001)-1) / size(V,1))+1); % non-zero eigenvalues
        eigenValues = E(find(E>0.000001));
        vv=[eigenValues transpose(eigenVectors)];
        vvs = sortrows(vv,-1);
        eigenVectors = transpose(vvs(:,2:end));
        eigenValues = vvs(:,1);
        perc_expl = 100*eigenValues/sum(eigenValues);
        
        % Get the eigenvalues containing XX% of the total variance
        var = 99;
        i = 1;
        s = sum(perc_expl(1:i));
        while ( s<var )
            i = i+1;
            s = sum(perc_expl(1:i));
        end
        d = i;
        
        eVal = eigenValues(1:d);
        eVect = eigenVectors(:,1:d);
        
        U = Psi*eVect*inv(diag(sqrt(eVal)));
        
        % Encode training data
        Z = transpose(transpose(U)*X_norm);
        
    case 'kernel'
        
        % disp('Performing kernel supervised PCA')
        % disp('The input matrix is the kernel matrix of input data')
        
        % Kernel is double centred
        X_norm = H*X_norm;
        
        % Computation of Q
        Q = X_norm*L*X_norm;
        
        % Computation of eigenvectors and eigenvalues
        [V, E] = eig(full(Q),X_norm);
        
        % extract eigenvectors and sort by eigenvalues
        eigenVectors = V(:,floor((find(E>0.000001)-1) / size(V,1))+1); % non-zero eigenvalues
        eigenValues = E(find(E>0.000001));
        vv=[eigenValues transpose(eigenVectors)];
        vvs = sortrows(vv,-1);
        eigenVectors = transpose(vvs(:,2:end));
        eigenValues = abs(vvs(:,1));
        perc_expl = 100*eigenValues/sum(eigenValues);
        
        % Get the eigenvalues containing XX% of the total variance
        var = 99;
        i = 1;
        s = sum(perc_expl(1:i));
        while ( s<var )
            i = i+1;
            s = sum(perc_expl(1:i));
        end
        d = i;
        
        % Encode training data
        Z = transpose(transpose(eigenVectors(:,1:d))*X_norm);
        
        eVal = eigenValues(1:d);
        eVect = eigenVectors(:,1:d);
        U = eVect;
        
    otherwise
        error('Type not implemented')
end

end

