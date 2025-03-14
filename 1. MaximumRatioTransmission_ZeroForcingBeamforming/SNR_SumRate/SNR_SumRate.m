%==========================================================================
close all;
clear all;

%Simulation parameters
NumberAntennas = [4 12];            %Number of transmit antennas
K = 4;                              %Number of users
PdB = -10:1:20;                     %Power 
P = 10.^(PdB/10); 
ChannelRelization = 2;              %Number of channel realizations 
channelVariances = [1 1 1 1];
weights = [1 1 1 1]'; ones(K,1);    %User weights


%==========================================================================
%Beamforming Techniques
sumRateMRT = zeros(length(P),ChannelRelization,length(NumberAntennas));
sumRateZFBF = zeros(length(P),ChannelRelization,length(NumberAntennas));

for n = 1:length(NumberAntennas)
    N = NumberAntennas(n);
    
    %Rayleigh fading channel
    RayleighChannel = (randn(K,N,ChannelRelization)+1i*randn(K,N,ChannelRelization))/sqrt(2);
    
    for m = 1:ChannelRelization
        %Channel matrix
        H = repmat(sqrt(channelVariances)',[1 N]) .* RayleighChannel(:,:,m);
        
        %Normalized beamforming vectors of ZFBF
        w_ZFBF = function_ZFBF(H);
        
        %Normalized beamforming vectors  MRT
        w_MRT = function_MRT(H);
        
        for i = 1:length(P)
            %Power allocation 
            rhos = diag(abs(H*w_MRT).^2)';
            powerAllocation_MRT = function_PowerAllocation(rhos,P(i),weights);
            
            rhos = diag(abs(H*w_ZFBF).^2)';
            powerAllocationw_ZFBF = function_PowerAllocation(rhos,P(i),weights);
            
            %Sum rate
            W = kron(sqrt(powerAllocation_MRT),ones(N,1)).*w_MRT;
            
            channelGains = abs(H*W).^2;
            signalGains = diag(channelGains);
            interferenceGains = sum(channelGains,2)-signalGains;
            
            rates = log2(1+signalGains./(interferenceGains+1));
            sumRate_MRT(i,m,n) = weights'*rates;  
            
            W = kron(sqrt(powerAllocationw_ZFBF),ones(N,1)).*w_ZFBF;
            
            channelGains = abs(H*W).^2;
            signalGains = diag(channelGains);
            interferenceGains = sum(channelGains,2)-signalGains;
            
            rates = log2(1+signalGains./(interferenceGains+1));
            sumRate_ZFBF(i,m,n) = weights'*rates;              
        end    
    end  
end


%==========================================================================
%Simulation results
for n = 1:length(NumberAntennas)
    figure; hold on; grid on
    plot(PdB,mean(sumRate_ZFBF(:,:,n),2),'b-o','LineWidth',1);
    plot(PdB,mean(sumRate_MRT(:,:,n),2),'r-*','LineWidth',1);

    title(['N = ' num2str(NumberAntennas(n)) ' antennas']);
    xlabel('Average SNR [dB]')
    ylabel('Average Sum Rate [bit/channel use]');
    legend('ZFBF','MRT');
end