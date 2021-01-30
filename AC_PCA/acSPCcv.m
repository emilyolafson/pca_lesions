function [ obj ] = acSPCcv( X, Y, X4Y, c1, c2s, v_ini, v_substract, kernel, bandwidth, opts)
%% 
% acSPC -- Perform cross-validation to tune the sparsity parameter c2 in 
% function acSPC

% input 
% X: a the n by p data matrix, where n is the number of samples, p is the 
% number of variables
% Missing values in X should be labeled as NaN. If a whole sample in X is
% missing, it should be removed.
% Y: the n by q confounder matrix, where n is the number of samples, q is 
% the number of confounding factors.  
% Missing values in Y should be labeled as NaN.
% X4Y: 'd' or a matrix of the same dimension as X, this matrix is used to 
% calculate the empirical Hilbert Schmidt criterion. The option 'd' means 
% default and uses X.
% c1: 'd' or a non-negative scalar. Tuning parameter for v'X'KXv. The option 
% 'd' means default and uses v_ini'X4Y'KX4Yv_ini. 
% c2s: a vector of non-negative tuning parameters controlling sparsity.
% v_ini: the initial v. Recommended to be the estimate of the non-sparse 
% version, for a "warm" start.
% v_substract: 'd' or a p by k matrix, the principal components 
% to be subtracted, where k is the number of PCs to be substracted. The 
% option 'd' means default and does not substract anything, i.e. calculate 
% the first principal component.
% kernel: the kernel to use, should be either 'linear', 'gaussian'.
% bandwidth: bandwidth for gaussian kernel. Provide any number for 'linear'
% kernel, it won't affect the result.
% opts: optional, some other options:
% opts.centerX: center the columns in X. Default is 1(True).
% opts.centerY: center the columns in Y. Default is 1(True).
% opts.scaleX: scale the columns in X to unit standard deviation. Default is 0(False).
% opts.scaleY: scale the columns in Y to unit standard deviation. Default is 0(False).
% opts.fold: the fold number for cross-validation. Default is 10.
% opts.plot: 1 or 0. 1: generates the c2 vs. mean squared error(MSE) plot. 
% Default is 1.
% opts.quiet: 1 or 0. 0: output the progress of the program. Default is 0. 
% opts...: other parameters for convergence.

% output
% obj:
% obj.mse: a vector of MSEs for c2s. Mean is taken for the sum of squared 
% errors in each fold.
% obj.c2s: the input c2s.
% obj.best_c2: c2 with the smallest MSE}
% obj.mse_sd: a vector of standard deviations for the MSEs. 
% Standard deviation across the folds is calculated.} 
%%

if nargin < 10
	opts.centerX = 1;
    opts.centerY = 1; 
    opts.scaleX = 0;
	opts.scaleY = 0;
    opts.fold = 10;
	opts.plot = 1; 
	opts.quiet = 0; 
    opts.maxiter = 25;
    opts.delta = 10^-4;
    opts.fold = 10;
    opts.plot = 1;
    opts.quiet = 0;
end

% check whether a whole row in X is missing
Xmis = sum(~isnan(X), 2);
if (sum(Xmis==0))
    error('some samples in X is missing');
end
[nX, ~] = size(X);
[nY, ~] = size(Y);

% check whether the number of samples in X and Y match
if (nX~=nY)
    error('The numbers of samples in X and Y do not match')
end
if (~strcmp(X4Y,'d'))
    % check whether the number of samples in X4Y and Y match
    [nX4Y, ~] = size(X4Y);
    if (nX4Y~=nY)
        error('The numbers of samples in X4Y and Y do not match')
    end
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
if (strcmp(c1,'d'))
    if (strcmp(X4Y,'d'))
        c1 = v_ini'*X'*K*X*v_ini;
    else
        X4Y(isnan(X4Y)) = mean(mean(X4Y, 'omitnan'));
        c1 = v_ini'*X4Y'*K*X4Y*v_ini;
    end
end

% check whether c1 is non-negative
if (c1<0)
    error('c1 should be non-negative')
end

% check whether c2s is non-negative
if (sum(c2s<0))
    error('c2s should be non-negative')
end

lab = unidrnd(opts.fold, size(X));
difsqr = nan(opts.fold, length(c2s));

opts1 = [];  
opts1.centerX = opts.centerX;
opts1.centerY = opts.centerY; 
opts1.scaleX = opts.scaleX;
opts1.scaleY = opts.scaleY;
opts1.maxiter = opts.maxiter;
opts1.delta = opts.delta;
opts1.filter = 0;    
for f = 1:opts.fold
    if (~opts.quiet)
        disp(['Running fold ' num2str(f)]);
    end
    X_cv = X;
    X_cv(lab==f) = mean(mean(X_cv(lab~=f)));
    for i = 1:length(c2s)
        c2 = c2s(i);
        if (~opts.quiet)
            disp([num2str(round(i/length(c2s)*100)) '% completed']);
        end
        if (c2==0)
            difsqr(f, i) = sum(X(lab==f).^2);  
        else
            if (strcmp(X4Y,'d'))
                tmp = acSPC( X_cv, Y, X, c1, c2, v_ini, v_substract, kernel, bandwidth, opts1);
            else
                tmp = acSPC( X_cv, Y, X4Y, c1, c2, v_ini, v_substract, kernel, bandwidth, opts1);
            end
            v = tmp.v;
            u = tmp.u;
            tmp = u'*X_cv*v*u*v';
            difsqr(f, i) = sum((tmp(lab==f) - X(lab==f)).^2);
        end
    end
end

mse = mean(difsqr, 1);
mse_sd = std(difsqr);
[~,I] = min(mse);
best_c2 = c2s(I);

% the c2s vs. mse plot
if (opts.plot) 
    hax=axes; 
    plot(c2s, mse, 'bo--');
    hold on 
    line([best_c2 best_c2],get(hax,'YLim'),'Color',[1 0 0])
    xlabel('c2s') 
    ylabel('Mean Squared Error') 
end


obj.mse = mse;
obj.mse_sd = mse_sd;
obj.c2s = c2s;
obj.best_c2 = best_c2;