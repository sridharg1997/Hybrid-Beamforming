Nt = 64;
NtRF = 4;

Nr = 16;
NrRF = 4;

rng(4096);
c = 3e8;
fc = 28e9;
lambda = c/fc;
txarray = phased.PartitionedArray(...
    'Array',phased.URA([sqrt(Nt) sqrt(Nt)],lambda/2),...
    'SubarraySelection',ones(NtRF,Nt),'SubarraySteering','Custom');
rxarray = phased.PartitionedArray(...
    'Array',phased.URA([sqrt(Nr) sqrt(Nr)],lambda/2),...
    'SubarraySelection',ones(NrRF,Nr),'SubarraySteering','Custom');

fc = 28e9;    % carrier in Hz
dlySpread = 50e-9;  % delay spread in seconds

% Channel Parameters

Ncl = 6;
Nray = 8;
Nscatter = Nray*Ncl;


snr_param = -40:5:0;
Nsnr = numel(snr_param);
Ns_param = [1 2 4];
NNs = numel(Ns_param);


Ropt = zeros(Nsnr,NNs);
Rhyb = zeros(Nsnr,NNs);
Niter = 1000;
%%
for m = 1:Nsnr
    snr = db2pow(snr_param(m));
    for n = 1:Niter
        % Channel realization
        [H,At,Ar] = channel_generation(Nt,Nr,Ncl,Nray);
        
        for k = 1:NNs
            Ns = Ns_param(k);
            % Compute optimal weights and its spectral efficiency
            [Fopt,Wopt] = OptimalWeights(H,Ns,1/snr);
            Ropt(m,k) = Ropt(m,k)+ SpectralEfficiency(H,Fopt,Wopt,Ns,snr);
            
            % Compute hybrid weights and its spectral efficiency
            [Fbb,Frf,Wbb,Wrf] = OMPWeights(H,Ns,NtRF,At,NrRF,Ar,1/snr);
            Rhyb(m,k) = Rhyb(m,k) + SpectralEfficiency(H,Frf*Fbb,Wrf*Wbb,Ns,snr);
        end
    end
end

Ropt = Ropt/Niter;
Rhyb = Rhyb/Niter;
figure();
plot(snr_param,Ropt(:,1),'--sr',snr_param,Rhyb(:,1),'-sr',...
    snr_param,Ropt(:,2),'--b',snr_param,Rhyb(:,2),'-b', ...
    snr_param,Ropt(:,3),'--g',snr_param,Rhyb(:,3),'-g');
   %
xlabel('SNR (dB)');
ylabel('Spectral Efficiency (bits/s/Hz');
legend('Ns=1 Optimal', 'Ns=1 Hybrid','Ns=2 Optimal','Ns=2 Hybrid', ...
    'Ns = 4 Optimal ','Ns = 4 Hybrid','Location','best');
grid on;


function R = SpectralEfficiency(H,F,W,Ns,snr)
A = W(:,1:Ns)'*H*F(:,1:Ns);
R = log2(det(eye(Ns)+snr/Ns*(real(W(:,1:Ns)'*W(:,1:Ns))\real(A'*A))));

end



