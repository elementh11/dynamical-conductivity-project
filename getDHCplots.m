
alat = 9.5632;
a1 = [   1.000000   0.000000   0.000000 ] * alat;
a2 = [   -0.500000   0.866025   0.000000 ] * alat;  
a3 = [   0.000000   0.000000   4.209200 ] * alat;  

vol = abs(dot(cross(a1,a2),a3));
b1 = 2 * pi * cross(a2,a3) / vol;
b2 = 2 * pi * cross(a3,a1) / vol;
b3 = 2 * pi * cross(a1,a2) / vol;

prefacSH = norm(b1) * norm(b2) / (2*pi);
prefacAH = norm(b1) * norm(b2) / (2*pi)^2;

load kfile.mat;
weight = datak.weightlist(2);%<--- assume equal weight for simplicity



load hallcond.mat;

dhckwn = data.dynHallCond;
sdhckwn = data.dynSpinHallCond;
dimfreq = size(dhckwn, 2);
numk = size(dhckwn, 1);
dim = size(dhckwn, 3);

dhc = zeros(1, dimfreq);
sdhc = zeros(1, dimfreq);
freqc = 0;
for freq = 1: dimfreq
    freqc = freqc + 1;
    dhckn = squeeze(dhckwn(:, freq, :));
    dhc(freqc) = sum(sum(dhckn)) * weight * prefacAH;
    
    sdhckn = squeeze(sdhckwn(:, freq, :));
    sdhc(freqc) = sum(sum(sdhckn)) * weight * prefacSH;
end
datap.dhc = dhc;
datap.sdhc = sdhc;
save('dhc-test','datap','-v7.3');







