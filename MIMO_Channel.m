function [H,AT,AR] = MIMO_Channel(N_t,N_r,N_c,N_ray)

%input: the numbers of transmit antennas and receive antennas
%output: the realized channel, codebook for vrf and wrf of OMP method 



E_aoa = 2*pi* rand(N_c,1);                                % The mean value of the cluster, obey the uniform distribution of (0,2*pi)

sigma_aoa = 5*pi/180;                                    % The angle is expanded to 5бу, converted to radians, which is the standard deviation
b = sigma_aoa/sqrt(2);                                    % Find the corresponding b according to the standard deviation, the scale parameter
a = rand(N_c,N_ray)-0.5;                                  % Generate a uniformly distributed random number sequence in the interval (-0.5, 0.5);
aoa = repmat(E_aoa,1,N_ray)-b*sign(a).*log(1-2*abs(a));   % Generate a random number sequence that conforms to the Laplace distribution 
                                                          % (each column represents the angle of each cluster)
aoa = sin(aoa);

%AOD

E_aod = 2*pi* rand(N_c, 1);                               %The mean value of the cluster, obey the uniform distribution of (0,2*pi)
sigma_aod = 5*pi/180;                                    %The angle is expanded to 5бу, converted to radians, which is the standard deviation
b = sigma_aod/sqrt(2);                                    % Find the corresponding b according to the standard deviation, the scale parameter
a = rand(N_c,N_ray)-0.5;                                  %Generate a uniformly distributed random number sequence in the interval (-0.5, 0.5);
aod = repmat(E_aod,1, N_ray)-b*sign(a).*log(1-2*abs(a));   %Generate a random number column conforming to the Laplace distribution 
                                                            %(each row represents the angle of each cluster))
aod = sin(aod);


signature_t = 1i*pi*[0:(N_t-1)]';
                         
signature_r =  1i*pi*[0:(N_r-1)]';
                                                         % Prepare for the next signature construction

H_ray = zeros(N_r, N_t, N_c, N_ray);
H_cl = zeros(N_r, N_t, N_c);


for i= 1: N_c
    for m = 1: N_ray
        H_ray(:,:,i,m)=complex(randn(1),randn(1))/sqrt(2)*exp((aoa(i,m)*signature_r))*exp((aod(i,m)*signature_t))'/sqrt(N_t*N_r); 
           
     end
end

H_cl = sum(H_ray, 4);
     
H(:,:) = sqrt(N_t*N_r/N_c/N_ray)*sum(H_cl,3);   
    
    
aod = aod(:).';
aoa = aoa(:).';
A = kron(aod,signature_t);
AT = 1/sqrt(N_t)*exp(A);
A = kron(aoa,signature_r);
AR = 1/sqrt(N_r)*exp(A);
   
