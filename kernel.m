function [k] = kernel(X1, X2, kernelChoice, sigma)
if nargin < 4, sigma = 1; end; %Standard deviation
if nargin < 3, kernelChoice = 1; sigma = 1; end; %Standard deviation
%Computes kernel values of vectors in respective COLUMNS in X1 and X2.
%Yields row vector of kernel values.


if kernelChoice == 1
    k = sum(X1.*X2); %Linear kernel = simply the dot product
elseif kernelChoice == 2
    normSquare = sum((X1-X2).*(X1-X2),1);
    k = exp(-1*normSquare/(2*sigma^2)); %Gaussian kernel
elseif kernelChoice == 3
    d = sigma; %The 4th argument is interpreted as the polynomial parameter
    k = (sum(X1.*X2) + 1).^d; %Polynomial kernel
else
    error('kernelChoice must be either 1=Linear or 2=Gaussian or 3=Polynomial');
end %if kernelChoice == 1