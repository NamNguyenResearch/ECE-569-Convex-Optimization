function w_MRT = function_MRT(H,D)
%Calculates the maximum ratio transmission (MRT) beamforming vectors 

%Number of users
Kr = size(H,1);

%Number of antennas
N = size(H,2);

if nargin<2
    D = repmat( eye(N), [1 1 Kr]);
end


%==========================================================================
%Pre-allocation of MRT beamforming
w_MRT = zeros(size(H'));

%MRT beamforming
for k = 1:Kr
    channelvector = (H(k,:)*D(:,:,k))';             %Useful channel
    w_MRT(:,k) = channelvector/norm(channelvector); %Normalization 
end