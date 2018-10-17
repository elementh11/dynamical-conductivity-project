% ---- k, n, and \omega resolved (spin) Hall conductivity

ef = -0.6547;
T = 0;
kb = 1.38e-23;  qelec = 1.6e-19; kt = kb / qelec * T;
smallDelta = 1e-10;

freqlist = linspace(0, 1, 200);
dimf = length(freqlist);

alat = 9.5632;
a1 = [   1.000000   0.000000   0.000000 ] * alat;
a2 = [   -0.500000   0.866025   0.000000 ] * alat;
a3 = [   0.000000   0.000000   4.209200 ] * alat;

vol = abs(dot(cross(a1,a2),a3));
b1 = 2 * pi * cross(a2,a3) / vol;
b2 = 2 * pi * cross(a3,a1) / vol;
b3 = 2 * pi * cross(a1,a2) / vol;


load hfile.mat;
matrices = datah.matrices;
nrpts = datah.nrpts;
dim = datah.num_wann;

load kfile.mat;
kpoints = datak.kpoints;
nk = size(kpoints,1);

sx = zeros(dim);
sy = zeros(dim);
sz = zeros(dim);
for c = 1:dim/2
    upind = 2*c - 1;
    dnind = 2*c;
    sz(upind,upind) = +1;
    sz(dnind,dnind) = -1;
    sx(upind,dnind) = 1;
    sx(dnind,upind) = 1;
    sy(upind,dnind) = -1i;
    sy(dnind,upind) = +1i;
end


energy = zeros(nk, dim);
dhc = zeros(nk, dimf, dim);
dshc = zeros(nk, dimf, dim);
for kc = 1:nk
    k = kpoints(kc,1:3);
    realk = k(1)*b1 + k(2)*b2 + k(3)*b3;
    hamWk = zeros(dim);
    vop_Wx = zeros(dim);
    vop_Wy = zeros(dim);
    for counter = 1:nrpts
        matrix = matrices(counter);
        delta = matrix.disp;
        realdisp = delta(1) * a1 + delta(2) * a2 + delta(3) * a3;
        ham = matrix.ham;
        hcontr = (ham * exp(1i* sum(conj(realk).*realdisp))* matrix.deg);
        vxcontr = 1i * realdisp(1) * hcontr;
        vycontr = 1i * realdisp(2) * hcontr;
        hamWk = hamWk + hcontr;
        vop_Wx = vop_Wx + vxcontr;
        vop_Wy = vop_Wy + vycontr;
    end
    hamWk = .5 * (hamWk + hamWk');
    [Uk, ek] = eig(hamWk); ek = diag(ek);
    sop_Wx = (vop_Wx * sz + sz * vop_Wx) / 2;
    
    vop_x = Uk' * vop_Wx * Uk;
    vop_y = Uk' * vop_Wy * Uk;
    sop_x = Uk' * sop_Wx * Uk;
    
    freqc = 0;
    for freq = freqlist
        freqc = freqc + 1;
        for n = 1:dim
            hc = 0;
            shc = 0;
            fnk = 1 / (exp((ek(n) - ef)/kt) + 1);
            for m = 1:dim
                if abs(ek(n) - ek(m)) > 1e-2
                    fmk = 1 / (exp((ek(m) - ef)/kt) + 1);
                    hc = hc + (fnk - fmk) * imag(vop_x(n,m) * vop_y(m,n))/((ek(n)-ek(m))^2 - (freq+1i*smallDelta)^2);
                    shc = shc + (fnk - fmk) * imag(sop_x(n,m) * vop_y(m,n))/((ek(n)-ek(m))^2 - (freq+1i*smallDelta)^2);
                end
            end
            dhc(kc, freqc, n) = hc;
            dshc(kc, freqc, n) = shc;
            energy(kc,n) = ek(n);
        end
    end
    if mod(kc, 5000)==0
        prog = fopen('progress.txt','a');
        fprintf(prog,'calculating for k=%d/%d \n',kc,nk);
        fclose(prog);
    end
    
end
data.energy = energy;
data.dynHallCond = dhc;
data.dynSpinHallCond = dshc;
save('hallcond','data');


