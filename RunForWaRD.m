% Given the Observation, the noise variance, the PSD of the input, the
% blurring function and the wavelet system setup (defined in ExptSetup and
% ForWaRDSetup), this program now runs the actual ForWaRD algorithm to perform
% the deconvolution.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Wiener Estimate is already calculated in ForWaRD Setup
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Regularized inversion in ForWaRD. Fourier Processing.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

RegInvfilter = (conj(G)) ./ ((abs(G).^2 ) + regparam);
xford = real(ifft2(RegInvfilter .* fft2(Y)));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Wavelet domian wiener filtering 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Take the RDWT of the xford
[xfordw, xfords] = mrdwtcycle2D(xford,WaveletFiltType1,DecompLevels);

% Calculate the variances of the noise at different levels
tempfilt = real(ifft2(RegInvfilter));
[tempw,temps] = mrdwtcycle2D(tempfilt,WaveletFiltType1,DecompLevels);
for lev = 1:DecompLevels
  for sub = 1:3
    noisevarw(lev,sub) = matrixnorm(tempw(:,:,lev,sub),2)*sigma;
  end
end
noisevars = matrixnorm(temps,2)*sigma;
clear tempw temps 

% Step 1 in WDWF: Get the hard thresholded estimate

for lev = 1:DecompLevels
  for sub = 1:3
    tempw(:,:,lev,sub)=HardTh(xfordw(:,:,lev,sub),...
	ThreshFactor*noisevarw(lev,sub)); %thresh difft wavelet subbands
  end
end
temps = HardTh(xfords,ThreshFactor*noisevars); % thresh scaling subband
clear xfordw xfords
xref = mirdwtcycle2D(tempw,temps,WaveletFiltType1,DecompLevels);
clear tempw temps

% Use the estimate to perform wavelet domain wiener filtering.

% Calculate the variances of the noise at different levels
[tempw,temps] = mrdwtcycle2D(tempfilt,WaveletFiltType2,DecompLevels);
for lev = 1:DecompLevels
  for sub = 1:3
    noisevarw(lev,sub) = matrixnorm(tempw(:,:,lev,sub),2)*sigma;
  end
end
noisevars = matrixnorm(temps,2)*sigma;
clear tempw temps tempfilt

[xfordw, xfords] = mrdwtcycle2D(xford,WaveletFiltType2,DecompLevels);
[xrefw, xrefs] = mrdwtcycle2D(xref,WaveletFiltType2,DecompLevels);

% Do the wavelet domain wiener filtering

for lev = 1:DecompLevels
  for sub = 1:3
    tempw(:,:,lev,sub)=xfordw(:,:,lev,sub) .*...
	((xrefw(:,:,lev,sub).^2)./...
	(xrefw(:,:,lev,sub).^2 + noisevarw(lev,sub)^2)); 
  end
end
temps =xfords .* ((xrefs.^2)./(xrefs.^2 + noisevars^2)); 
clear xford xfordw xfords xrefw xrefs
xward = mirdwtcycle2D(tempw,temps,WaveletFiltType2,DecompLevels);
clear tempw temps 






