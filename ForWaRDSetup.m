% This part sets up the experiment according to the specifications laid out
% in the ExptSetup.m file

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input Signal
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if InputSignal == 1
  load camera; X = camera;
  clear camera
elseif InputSignal == 2
  load lenna256; X = lenna;
  clear lenna
elseif InputSignal == 3
  load theater; X = theater(50:305, 10:265);
  clear theater
elseif InputSignal == 4
  load boats; X = boats;
  clear boats
elseif InputSignal == 5
  load birthday; X = birthday;
  clear birthday
else   
  error('You did not specify any valid InputSignal')
end

% The input signal is normalized between 0 and 1
X = X/255;
N = length(X);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Blurring function
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if Blur == 1
  q = 9;
  hrow = [ ones(1,q) , zeros(1, N-q)];
  hrow = [ ones(1,q) , zeros(1, N-q)];
  hrow = hrow(:);  
  hcol = hrow;
  ConvMat = (hcol * hrow.') / q^2;
elseif Blur == 2
  q = floor(input('What is the blur size? Give a positive integer\n'));
  if q <= 0
    error('Enter a positive blur size')
  end
  hrow = [ ones(1,q) , zeros(1, N-q)];
  hrow = hrow(:);  
  hcol = hrow;
  ConvMat = (hcol * hrow.') / q^2;
elseif Blur == 3
  ConvMat = zeros(N,N);
  ConvMat(1:4,1:3) = 1;
  ConvMat(1:3,1:4) = 1;
  ConvMat(5,1:2) = 1;
  ConvMat(1:2,5) = 1;
  
  ConvMat(N-3:N,1:2) = 1;
  ConvMat(N-4:N,1:3) = 1;
  ConvMat(N-1:N,1:4) = 1;
  ConvMat(N,1:5) = 1;

  ConvMat(1:2,N-3:N) = 1;
  ConvMat(1:3,N-2:N) = 1;
  ConvMat(1:4,N-1:N) = 1;
  ConvMat(1:5,N) = 1;
  
  ConvMat(N,N-3:N) = 1;
  ConvMat(N-1,N-2:N) = 1;
  ConvMat(N-2,N-1:N) = 1;
  ConvMat(N-3,N) = 1;
  ConvMat = ConvMat/(pi*17);
elseif Blur == 4
  % The r used below is the same as that used in the Kalifa paper.
  r = 1;
  fftrow = [ones(1,N/4) (2^r) * (abs( 2*((N/4+1):(N/2))/(N)-1)).^r ] ; 
  sum1 = norm([ fftrow  fliplr(fftrow)]);
  fftrow = [ fftrow  fliplr(fftrow)];
  sum2 = norm(fftrow(:));
  fftrow = fftrow*sum1/sum2;  
  hrow   = real(ifft(fftrow));
  hrow = hrow(:);
  clear fftcol
  hcol = hrow;
  ConvMat = hcol * hrow.';
elseif Blur == 5
   q = 5;
   hrow = [ ones(1,q) , zeros(1, N-q)];
   hrow = hrow(:);  
   hcol = hrow;
   ConvMat = (hcol * hrow.') / q^2;
elseif Blur == 6
   q = 6;
   hrow = [ ones(1,q) , zeros(1, N-q)];
   hrow = hrow(:);  
   hcol = hrow;
   ConvMat = (hcol * hrow.') / q^2;
elseif Blur == 7
   q = 7;
   hrow = [ ones(1,q) , zeros(1, N-q)];
   hrow = hrow(:);  
   hcol = hrow;
   ConvMat = (hcol * hrow.') / q^2;
elseif Blur == 8
   q = 3;
   hrow = [ ones(1,q) , zeros(1, N-q)];
   hrow = hrow(:);  
   hcol = hrow;
   ConvMat = (hcol * hrow.') / q^2;
else
  error('You did not specify any valid Blur')
end
clear hcol hrow 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Noise Setup
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tempvar = real(ifft2(fft2(X) .* fft2(ConvMat)));
% orgsigma is the standard deviation of noise
OrgSigma = sqrt(norm(tempvar(:)-mean(tempvar(:)),2)^2 /(N^2*10^(BSNR/10)));
n = OrgSigma*randn(N,N);
clear tempvar
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Blurred and noisy observation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Y = real(ifft2(fft2(X) .* fft2(ConvMat))) + n ;
clear n
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Defining the wavelet systems for WaRD. Redundant wavelet transforms are
% used.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[tempW, tempS] = mrdwtcycle2D(Y,WaveletFiltType1,DecompLevels);
clear tempS

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Noise variance estimation using median estimator
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if Variance == 'Estimate'
  temp = abs(tempW(:,:,DecompLevels,3));
  sigma = median(temp(:))/.67;
else
  sigma = OrgSigma;
end
clear temp tempW

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Power spectral density estimation
% This might not be robust for certain images
%
% Reference: A.D. Hillery and R.T. Chin, ``Iterative {W}iener filters
% for image restoration,'' {\em IEEE Trans. Signal Processing}, vol.~39,
% pp.~1892--1899,  Aug. 1991.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if PSD == 'Estimate'
  % Note: We add an offset to the PSD estimate to stabilize it at frequencies
  % where the convolution operator response goes to zero.
  Iterations = 10;
  delta = 5*10^-7;
  disp(char(strcat('Using an adhoc parameter delta = ',num2str(delta),{' '},'in the PSD estimation for Wiener filter. For best Wiener filtering results, try different values for delta. Ideally, delta = 0.')));
  Pyy = abs(fft2(Y)).^2;
  G = fft2(ConvMat);
  Phh = abs(G).^2;
  % Quick estimate to begin with
  sqrdtempnorm = ((norm(Y(:)-mean(Y(:)))^2 -N*(N-1)*sigma^2)/(norm(ConvMat(:),1))^2); % using minkowski's ineq
  tempreg = sigma^2*N*N/(sqrdtempnorm);
  GenericRegInv = (conj(G)) ./ ((abs(G).^2 ) + tempreg );
  Pxx = abs(fft2(Y) .* GenericRegInv).^2;
  Pxx = Pyy/((norm(ConvMat(:),1))^2);
  for i = 1:Iterations
    M = (conj(G) .* Pxx .* fft2(Y)) ./ (Phh .* Pxx + N^2*(sigma^2));
    PxxY = (Pxx .* N^2*sigma^2) ./ (Phh.*Pxx + N^2*(sigma^2));
    Pxx = PxxY + abs(M).^2;
    NormPxx(i) = norm(Pxx(:));
  end
    Pxx = (Pxx+delta*norm(Pxx(:)))*norm(Pxx(:))/norm(Pxx(:) +delta*norm(Pxx(:)));
  %Pxx = (Pxx+delta*max(Pyy(:))/(max(Phh(:)))) *norm(Pxx(:))/norm(Pxx(:) + delta*max(Pyy(:))/(max(Phh(:))));
  Weight = sqrt(Pxx);
else 
  G = fft2(ConvMat);  
  Pxx = abs(fft2(X)).^2;
  Weight = sqrt(Pxx);
end
clear Pyy Pxx Phh M PxxY i delta 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Wiener Estimate
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
WienerFilterEst = (conj(G) .* (Weight.^2)) ./ ((abs(G).^2 .* (Weight.^2)) + N^2*(sigma^2));
WienerEstimate  = real(ifft2(WienerFilterEst .* fft2(Y)));
xwiener = WienerEstimate;
clear WienerFilterEst Weight

if RegSearch == 'Y'
  regscalefactor = 10;
  RegMultFactor = 1;
  RegParam = [3.3*10.^(-1) 6.6*10.^(-1) 10.^(-0) 3.3*10.^(0) 6.6*10.^(-0)]*RegMultFactor;
  % normalization wrt samples and noise variance
  RegParam = RegParam(:); 
end
