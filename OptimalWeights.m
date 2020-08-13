function [Fopt,Wopt] = OptimalWeights(Hin,Ns,noisevar)

H = Hin;
[~,~,v] = svd(H);
Fopt = v(:,1:Ns);


Wopt = (((H*Fopt)'*(H*Fopt) + Ns*noisevar*eye(Ns))\((H*Fopt)'))';



