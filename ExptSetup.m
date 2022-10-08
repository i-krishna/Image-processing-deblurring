
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                            %
%                       ForWaRD Algorithm Implementation                     %
%                              (Image deblurring)                            %
%                                                                            %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% This code specifies the experimental setup required to demonstrate the
% performance of Fourier-Wavelet Regularized Deconvolution (ForWaRD)
% algorithm proposed by Ramesh Neelamani, Hyeokho Choi and Richard Baraniuk
% from Rice University.
% 
% Reference: R. Neelamani, H. Choi, and R. G. Baraniuk, ``ForWaRD: 
% Fourier-Wavelet Regularized Deconvolution for Ill-Conditioned Systems,''
% Submitted to IEEE Transactions on Image Processing (October 2001). 
% Check http://www.dsp.rice.edu/publications/ to
% download the reference.
%
% Your comments would be highly appreciated. Please send them to 
% neelsh@ece.edu
%
% This program uses some routines from the  Rice Wavelet toolbox.
%
% Copyright: All software, documentation, and related files in this 
%            distribution are Copyright (c) 2000 Rice University
%
% Permission is granted for use and non-profit distribution providing that this
% notice be clearly maintained. The right to distribute any portion for profit
% or as part of any commercial product is specifically reserved for the author.

%  addpath 'C:\Users\krishna\Desktop\masters-chse\4th sem - 5 subjects\seminar-Embedded Image Processing\Seminar - Forward algorithm\rwt-when u get error for forward1d due to windows new version\rwt-master\bin' 
%  addpath 'C:\Users\krishna\Desktop\masters-chse\4th sem - 5 subjects\seminar-Embedded Image Processing\Seminar - Forward algorithm\forward\WaRD-ver2.0\routines'

addpath(genpath('C:\Users\krishna\Desktop\masters-chse\4th sem - 5 subjects\seminar-Embedded Image Processing\Seminar - Forward algorithm'))
% Which Input signal would you like? Choose the corresponding number
% (1) Cameraman (default)
% (2) Lenna 
% (3) Theater
% (4) Boats
% (5) Birthday
% The images are grayscale images
InputSignal = 1;

% Which blurring function would you like?
% (1)  Box car blur (9x9 blur) (default)
% (2)  Adjustable Box car blur (you specify the amount of blur size later)
% (3)  Circular blur with radius 7
% (4)  Low pass filter used by Kalifa et al. in 
%      "Image Deconvolution in Mirror Wavelet bases", ICIP 98, pg. 565-569,
%      Chicago, 1998.
% (5)  Box car blur (5x5 blur) (default)
% (6)  Box car blur (6x6 blur) (default)
% (7)  Box car blur (7x7 blur) (default)
% (8)  Box car blur (3x3 blur) (default)

Blur = 1;

% What blurred-signal-to-noise ratio (in dB) would you like?
% If y = x*h + n is the deconvolution model, the BSNR = |norm(x*h)/norm(n)|^2
% default is 40 dB (* denotes convolution)

BSNR = 40;

% Two wavelet basis are required for implementing the wavelet-domain wiener
% filter (WDWF) used to perform the estimation. See daubcqf.m to set the
% type wavelet filter. 

WaveletFiltType1 = daubcqf(6,'min'); 
WaveletFiltType2 = daubcqf(2,'min'); 

% What threshold should be used in for hard thresholding (1st stage in WDWF)?
% default threshold is 3*sigma, sigma is standard deviation of noise.

ThreshFactor = 3.0;

% How many levels of decomposition do you want?
% The maximum possible is log(N), however 4 levels are often sufficient.

DecompLevels = 3;

% Do you want to estimate the power spectral density (PSD) of the input
% signal from the observation or do you want to use the original signal to
% get the PSD? The PSD is required for the Wiener filter.
% Choices are PSD = 'Estimate'
%             PSD = 'Original'
% NOTE: For arbitrary images, the estimation 

PSD = 'Estimate';

% Do you want to estimate the noise variance of the input
% signal from the observation or do you want to use the original signal to
% get the PSD? The noise variance is required for the Wiener filter and the 
% WaRD algorithm.
% Choices are Variance = 'Estimate'
%             Variance = 'Original'
%

Variance = 'Estimate';

% Do you want the regularization parameter to be automatically estimated?
% Choices are RegSearch = 'Y' (recommended, but about 5 times slower) 
%             RegSearch = 'N' (default)

RegSearch = 'N';

% If RegSearch = 'N', then you need to specify a regularization parameter.
% The choice is normalized with respect to the noise variance and average 
% signal power. Choice of RegParam is ignored if RegSearch = 'Y'.

RegParam = 1;

% This runs the experiment.
ForWaRDSetup

% This searches the regularization parameter.
RegParamSetup

% We are now in a position to run the ForWaRD algorithm.
RunForWaRD

% We are now in a position to display the results
DispResults













