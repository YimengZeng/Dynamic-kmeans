function W = ComputeWaveletCorrelations(Y)


% Add path to Matlab Wavelet correlation toolbox
addpath(genpath('/mnt/musk2/home/ksupekar/kaustubh/ASDvsTDNetworkAnalysis/Analysis_1_19ASD_19TD/Scripts/Step3-WaveletAnalysis/MatlabWaveletCorrelations/wmtsa-matlab-0.2.6'))

data = Y;
wavelet = 'la8';
J0 = 3;
boundary = 'reflection';

for i = 1 : size(data,2)
    for j = 1 : size(data,2)
        [WJtX, VJ1tX] = modwt(data(4:end,i), wavelet, J0, boundary);
        [WJtY, VJ1tY] = modwt(data(4:end,j), wavelet, J0, boundary);
        [wcor, CI_wcor] = modwt_wcor(WJtX, WJtY);
        correlationScaleOne(i,j) = wcor(1);
        correlationScaleTwo(i,j) = wcor(2);
        correlationScaleThree(i,j) = wcor(3);
    end
end
W = correlationScaleThree;

