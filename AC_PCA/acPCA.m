function [ obj ] = acPCA( X, Y, lambda, nPC, kernel, bandwidth, eval, opts)
%% 
% acPCA -- simultaneous dimension reduction and adjustment for confounding 
% variation.

% input 
% X: a the n by p data matrix, where n is the number of samples, p is the 
% number of variables
% Missing values in X should be labeled as NaN. If a whole sample in X is
% missing, it should be removed.
% Y: the n by q confounder matrix, where n is the number of samples, q is 
% the number of confounding factors.  
% Missing values in Y should be labeled as NaN.
% lambda: tuning parameter, non-negative
% nPC: number of principal components to compute
% kernel: the kernel to use, should be either 'linear', 'gaussian'.
% bandwidth: bandwidth for gaussian kernel. Provide any number for 'linear'
% kernel, it won't affect the result.
% eval: 1 or 0. 1: evaluate the significance of the PCs.
% opts: optional, some other options:
% opts.centerX: center the columns in X. Default is 1(True).
% opts.centerY: center the columns in Y. Default is 1(True).
% opts.scaleX: scale the columns in X to unit standard deviation. Default is 0(False).
% opts.scaleY: scale the columns in Y to unit standard deviation. Default is 0(False).
% opts.numPerm the number of permutations to evaluate the significance of the PCs. Default is 100. 
% opts.alpha the significance level. Default is 0.05. 
% If the kth PC is not significant, we don't consider the PCs after it.
% If the eigenvalue and variance explained by the PCs give inconsistent results, 
% we choose the maximum number of significant PCs.

% output
% obj:
% obj.v: the principal components, p by nPC matrix
% obj.Xv: the projected data, i.e. X times v
% obj.eigenX: eigenvalues of the PCs in X
% obj.varX: variance of the PCs in X
% obj.varX_perc: percentage of the variance explained by the PCs in X
% obj.eigenXperm: eigenvalues of the PCs in permutation. If eval=0, nan is
% returned
% obj.varXperm: variance of the PCs in permutation. If eval=0, nan is
% returned.
% obj.varXperm_perc: percentage of the variance explained by the PCs in
% permutation. If eval=0, nan is returned.
% obj.sigPC: number of significant PCs. If eval=0, nan is returned.
% ...: input parameters for the function
%%

if nargin < 8
    opts = [];  
    opts.centerX = 1; 
    opts.centerY = 1; 
    opts.scaleX = 1;
    opts.scaleY = 0; 
    opts.numPerm = 100;
    opts.alpha = 0.05;
end
% check whether a whole row in X is missing
Xmis = sum(~isnan(X), 2);
if (sum(Xmis==0))
    error('some samples in X is missing');
end
[nX, p] = size(X);
[nY, ~] = size(Y);
% check whether the number of samples in X and Y match
if (nX~=nY)
    error('The numbers of samples in X and Y do not match')
end
% check whether lambda is non-negative
if (lambda<0)
    error('lambda should be non-negative')
end
% center the X matrix
if (opts.centerX)
    X = X-repmat(mean(X,'omitnan'),nX,1);
end
% center the Y matrix
if (opts.centerY)
    Y = Y-repmat(mean(Y,'omitnan'),nY,1);
end
% scale the X matrix
if (opts.scaleX)
    Xsd = std(X,'omitnan');
    Xsd(Xsd==0) = 1;
    X = X./repmat(Xsd,nX,1);
end
% scale the Y matrix
if (opts.scaleY)
    Ysd = std(Y,'omitnan');
    Ysd(Ysd==0) = 1;
    Y = Y./repmat(Ysd,nY,1);
end
% input the missing values in X and Y with the mean
X(isnan(X)) = mean(mean(X, 'omitnan'));
Y(isnan(Y)) = mean(mean(Y, 'omitnan'));
K = calkernel(Y, kernel, bandwidth);
calAv = @(v) (X' - lambda*X'*K)*(X*v);
optseigs.isreal = 1; 
optseigs.issym = 1;
[V,D] = eigs(calAv,p,nPC,'LA', optseigs);
% the eigenvalues
eigenX = diag(D)';
% the projection
Xv = X*V;
% variance in Xv
varX = var(Xv, 0, 1);
% total variance in X
totvar = sum(var(X, 0, 1));
% percentage
varX_perc = varX/totvar;
% evaluate the significance of PCs
if (eval)
    eigenXperm = zeros(opts.numPerm, nPC);
    varXperm = zeros(opts.numPerm, nPC); 
    varXperm_perc = zeros(opts.numPerm, nPC); 
    for i = 1:opts.numPerm
        Xperm = X(randperm(nX),:);
        calAvperm = @(v) (Xperm' - lambda*Xperm'*K)*(Xperm*v);
        [Vperm,Dperm] = eigs(calAvperm,p,nPC,'LA', optseigs);
        % eigenvalue
        eigenXperm(i,:) = diag(Dperm)';
        % the projection
        Xpermv = Xperm*Vperm;
        % variance in Xpermv
        varXperm(i,:) = var(Xpermv, 0, 1)';
        % percentage
        varXperm_perc(i,:) = varXperm(i,:)/totvar;
    end
    labs1 = find(eigenX < quantile(eigenXperm, 1-opts.alpha));
    if (numel(labs1)==0)
      sigPC1 = 0;
    else
      sigPC1 = min(labs1)-1;
    end
    labs2 = find(varX < quantile(varXperm, 1-opts.alpha));
    if (numel(labs2)==0)
      sigPC2 = 0;
    else
      sigPC2 = min(labs2)-1;
    end
    sigPC = max(sigPC1, sigPC2);
else
    eigenXperm = nan;
    varXperm = nan; 
    varXperm_perc = nan; 
    sigPC = nan;
end
% output
obj.Xv = Xv;
obj.v = V;
obj.eigenX = eigenX;
obj.varX = varX;
obj.varX_perc = varX_perc;
obj.eigenXperm = eigenXperm;
obj.varXperm = varXperm; 
obj.varXperm_perc = varXperm_perc; 
obj.sigPC = sigPC;
obj.lambda = lambda;
obj.kernel = kernel;
if (strcmp(kernel,'gaussian'))
    obj.bandwidth = bandwidth;
end