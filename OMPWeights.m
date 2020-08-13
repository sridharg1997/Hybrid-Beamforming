function [Fbb,Frf,Wbb,Wrf] = OMPWeights(Hin,Ns,NtRF,At,NrRF,Ar,noisevar)


H = Hin;
[~,~,v] = svd(H);
Fopt = v(:,1:Ns);
[Nr,Nt] = size(H);

Frf = complex(zeros(Nt,NtRF));
Fres = Fopt;
for m = 1:NtRF
    Psi = At'*Fres;
    [~,k] = max(diag(Psi*Psi'));
    Frf(:,m) = At(:,k);
    Fbb = ((Frf(:,1:m))'*Frf(:,1:m))\Frf(:,1:m)'*Fopt;
    temp = Fopt-Frf(:,1:m)*Fbb;
    Fres = temp/norm(temp,'fro');
end
Fbb = sqrt(Ns)*Fbb/norm(Frf*Fbb,'fro');

Wmmse = ((Fbb'*Frf'*(H'*H)*Frf*Fbb+noisevar*Ns*eye(Ns))\Fbb'*Frf'*H')';

Wrf = complex(zeros(Nr,NrRF));
Wres = Wmmse;
Ess = (1/Ns)*eye(Ns);
Eyy = (H*Frf*Fbb)*Ess*(H*Frf*Fbb)' + noisevar*eye(Nr);

for m = 1:NrRF
    Psi = Ar'*(Eyy*Wres);
    [~,k] = max(diag(Psi*Psi'));
    Wrf(:,m) = Ar(:,k);
    Wbb = (Wrf(:,1:m)'*Eyy*Wrf(:,1:m))\(Wrf(:,1:m)'*Eyy*Wmmse);
    temp = Wmmse-Wrf(:,1:m)*Wbb;
    Wres = temp/norm(temp,'fro');
end
