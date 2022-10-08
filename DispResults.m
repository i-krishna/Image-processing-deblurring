
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Display the results of the ForWaRD algorithm
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

temp = Y;
if Blur == 1
  % the observation is shifted by (4,4) towards the top-left corner to minize
  % the ISNR.
  minshifti= 4;
  minshiftj= 4;
elseif (Blur == 2 | Blur == 5 | Blur == 6 |   Blur == 7| Blur == 8)    
  for shifti = floor(q/2):ceil(q/2)
    for shiftj = floor(q/2):ceil(q/2)
      convmat = 0* ConvMat;
      convmat(N-shifti+1,N-shiftj+1) = 1;  
      Y1 = real(ifft2(fft2(temp) .* fft2(convmat)));
      cost(shifti,shiftj) = norm(X(:)-Y1(:))^2;
    end    
  end
  tempcost = cost(floor(q/2):ceil(q/2),floor(q/2):ceil(q/2));
  [minshifti, minshiftj] = find(cost == min(tempcost(:)));
elseif Blur == 3
  minshifti= 1;
  minshiftj= 1;
else Blur == 4
  maxshifti = 2;
  maxshiftj = 2;
  minshifti = 2;
  minshiftj = 2;
  while (minshifti == maxshifti)|(minshiftj == maxshiftj)
    maxshifti = maxshifti + 3;
    maxshiftj = maxshiftj + 3;
    for shifti = 1:maxshifti
      for shiftj = 1:maxshiftj
	convmat = 0* ConvMat;
	convmat(N-shifti+1,N-shiftj+1) = 1;  
	Y1 = real(ifft2(fft2(temp) .* fft2(convmat)));
	cost(shifti,shiftj) = norm(X(:)-Y1(:))^2;
      end    
    end
    [minshifti, minshiftj] = find(cost == min(cost(:)));
  end  
end

convmat = 0* ConvMat;
convmat(N-minshifti+1,N-minshiftj+1) = 1;  
Y1 = real(ifft2(fft2(temp) .* fft2(convmat)));
clear temp convmat

figure(1)
subplot(2,2,2)
imagesc(X)
colormap('gray')
imrange = get(gca,'CLim');
imagesc(Y1)
axis('square')
rmaxis
title(strcat('Noisy image, SNR = ', num2str(BSNR),'dB'))

subplot(2,2,1)
imagesc(X, imrange)
axis('square')
rmaxis
title(strcat('Original image'))

subplot(2,2,3)
imagesc(WienerEstimate, imrange)
axis('square')
rmaxis
ISNRWiener = 10*log10(norm(X(:)-Y1(:))^2 /norm(WienerEstimate(:) - X(:))^2 );
% rounding off to two digits
ISNRWiener = round(100*ISNRWiener)/100;
SNRWiener = 10*log10(norm(X(:))^2 /norm(WienerEstimate(:) - X(:))^2 );
SNRWiener = round(100*SNRWiener)/100;
title(strcat('Wiener, ISNR =',num2str(ISNRWiener),'dB, SNR =',num2str(SNRWiener),'dB'))

subplot(2,2,4)
imagesc(xward, imrange)
axis('square')
rmaxis
ISNRWDWF = 10*log10(norm(X(:) - Y1(:))^2 /norm(xward(:) - X(:))^2 );
ISNRWDWF = round(100*ISNRWDWF)/100;
SNRWDWF = 10*log10(norm(X(:))^2 /norm(xward(:) - X(:))^2 );
SNRWDWF = round(100*SNRWDWF)/100;
title(strcat('ForWaRD, ISNR = ', num2str(ISNRWDWF), 'dB, SNR = ', num2str(SNRWDWF),'dB'))


