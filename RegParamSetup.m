% This program chooses an appropriate Tikhonov regularization parameter for
% the Fourier regularized inversion step. 

sqrdtempnorm = ((norm(Y(:)-mean(Y(:)))^2 -N*(N-1)*sigma^2)/(norm(ConvMat(:),1))^2); % using
                                                     % minkowski's ineq
tempreg = sigma^2*N*N/(sqrdtempnorm);

if RegSearch == 'N'
    % Nothing to do
    %regparam = RegParam*sigma^2*N*N/(norm(Y(:) -mean(Y(:)))^2);
    regparam = RegParam*sigma^2*N*N/(sqrdtempnorm);
    % regparam = 0.00010133; % for 40 db bsnr case
    disp(' ')
    disp(['Using a pre-set regularization parameter. For automatic estimation, set RegSearch = Y in ExptSetup.m file (recommended).'])
    disp(' ')
else 
  clear Actualmse Obsnorm Resultnorm

  GenericRegInv = (conj(G)) ./ ((abs(G).^2 ) + tempreg );
  GenericObs = real(ifft2(fft2(Y) .* GenericRegInv));
  for regcounter = 1:length(RegParam)
    regparam = RegParam(regcounter)*tempreg;  
    RunForWaRD
    % xward is the forward estimate
    Resultnorm(regcounter) = norm(xward(:));
    %disp(strcat('Testing regularization parameter number', {' '}, num2str(regcounter),{' '},'of', {' '}, num2str(length(RegParam))))
    Actualmse(regcounter) = norm(xward(:) - X(:));
    temp = real(ifft2(fft2(xward) .* fft2(ConvMat) .* GenericRegInv));
    Obsnorm(regcounter) = norm(temp(:) - GenericObs(:));
  end  
  clear GenericObs  GenericRegInv temp
  ActualmsedB = -20*log10(Actualmse/(norm(X(:))));
  [MaxMSE,ind] = max(ActualmsedB);
  MaxMSEReg = RegParam(ind);
  [a,ObsInd] = min(Obsnorm);
  if ((ObsInd == 1)| (ObsInd == length(RegParam)))
    disp('Increasing Reg Parameter range')
    if (ObsInd == 1)
      RegParam = RegParam/regscalefactor;
      disp((RegParam(:)).')
      RegParamSetup
    elseif (ObsInd == length(RegParam))
      RegParam = RegParam*regscalefactor;
      disp((RegParam(:)).')
      RegParamSetup      
    end
  else
    ActualmsedB = -20*log10(Actualmse/(norm(X(:))));
    [a,ind] = max(ActualmsedB);
    if MaxMSE < a
      MaxMSE = a;
      MaxMSEReg = RegParam(ind);
    end
    disp(strcat('Best possible performance is =',num2str(MaxMSE)));
    disp(strcat('Best Reg Parameter is =',num2str(MaxMSEReg)));
    
    [a,ObsInd] = min(Obsnorm);
    disp(strcat('Chosen Reg Parameter is =',num2str(RegParam(ObsInd))));
    disp(strcat('Chosen performance from Observation norm is =',num2str(ActualmsedB(ObsInd))));
    regparam = RegParam(ObsInd)*tempreg;
  end
end

