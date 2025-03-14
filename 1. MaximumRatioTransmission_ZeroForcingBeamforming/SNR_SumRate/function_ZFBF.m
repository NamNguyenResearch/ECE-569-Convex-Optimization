function w_ZFBF = function_ZFBF(H,D)
%Calculates the zero-forcing beamforming (ZFBF) vectors

%Number of users
Kr = size(H,1);

%Number of antennas
N = size(H,2);

if nargin<2
    D = repmat( eye(N), [1 1 Kr]);
end


%==========================================================================
%Pre-allocation of MRT beamforming
w_ZFBF = zeros(size(H'));

%ZFBF beamforming
for k = 1:Kr
    effectiveChannel = (H*D(:,:,k))';                                           %Effective channels
    channelInversion = effectiveChannel/(effectiveChannel'*effectiveChannel);   %Channel inversion
    w_ZFBF(:,k) = channelInversion(:,k)/norm(channelInversion(:,k));            %Normalization of zero-forcing direction
end