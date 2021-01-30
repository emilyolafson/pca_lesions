function [ obj ] = acPCAtuneLambda(X, Y, lambdas, nPC, kernel, bandwidth, anov, opts)
%% 
% acPCAcv -- performs cross-validation for acPCA to select the lambda
% parameter

% input 
% X: a the n by p data matrix, where n is the number of samples, p is the 
% number of variables
% Missing values in X should be labeled as NaN. If a whole sample in X is
% missing, it should be removed.
% Y: the n by q confounder matrix, where n is the number of samples, q is 
% the number of confounding factors.  
% Missing values in Y should be labeled as NaN.
% lambdas: a vector of tuning parameters, non-negative. If there is no 0,
% 0 will be added to lambdas.
% nPC: number of principal components to compute
% kernel: the kernel to use, should be either 'linear', 'gaussian'.
% bandwidth: bandwidth for gaussian kernel. Provide any number for 'linear'
% kernel, it won't affect the result.
% anov: 1 or 0. 1: the penalty term has the between groups sum of squares interpretation.
% opts: optional, some other options:
% opts.centerX: center the columns in X. Default is 1(True).
% opts.centerY: center the columns in Y. Default is 1(True).
% opts.scaleX: scale the columns in X to unit standard deviation. Default is 0(False).
% opts.scaleY: scale the columns in Y to unit standard deviation. Default is 0(False).
% opts.perc: the best lambda is defined to be the smallest lambda with 
% R(lambda)<=perc (if anov=1), or R(lambda)<=perc*R(lambda=0) (if anov=F) 
% in the nPC principal components. 
% opts.quiet: 1 or 0. Output the progress of the program. Default is 0.
% Default is 0.

% output
% obj:
% obj.ratio: a vector with the ratios, (v'X'KXv)/(v'X'Xv). Same length as lambdas.
% obj.best_lambda: the best lambda after cross-validation
% ...: input parameters for the function
%%

if nargin < 8
    opts = [];  
    opts.centerX = 1; 
    opts.centerY = 1; 
    opts.scaleX = 1;
    opts.scaleY = 1; 
    opts.perc = 0.001;
    opts.quiet = 0;
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
% check whether lambdas is non-negative
if (sum(lambdas<0))
    error('lambdas should be non-negative')
end

lambdas = sort(lambdas);
% add 0 to lambdas
if (sum(lambdas==0)==0)
    lambdas = [0 lambdas];
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

optseigs.isreal = 1; 
optseigs.issym = 1;
ratio = zeros(length(lambdas), nPC);
for lab = 1:length(lambdas)
    if (~opts.quiet)
        if (rem(lab,round(length(lambdas)/10))==0)
            disp([num2str(round(lab/length(lambdas)*100)) '% completed']);
        end
    end
    lambda = lambdas(lab);
    calAv = @(v) (X' - lambda*X'*K)*(X*v);
    [V,~] = eigs(calAv,p,nPC,'LA', optseigs);
    Xv = X*V;
    ratio(lab, :) = diag(Xv'*K*Xv)./diag(Xv'*Xv);
end

if (anov)
    thres = opts.perc;
else
    thres = max(ratio(1,:))*opts.perc;
end

tmp = find(sum(ratio<=thres, 2) == nPC);
if (numel(tmp)==0)
    warning('lambda is not large enough, increase max(lambdas) or increase perc');
else
    best_lambda = lambdas(min(tmp));
end

obj.ratio = ratio;
obj.best_lambda = best_lambda;
obj.lambdas = lambdas;
obj.kernel = kernel;
if (strcmp(kernel,'gaussian'))
    obj.bandwidth = bandwidth;
end