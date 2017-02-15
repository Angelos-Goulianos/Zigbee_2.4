
%This script implements an end to end  Zigbee transmitter-receiver under AWGN noise at 2.4 GHz
%considering the half-sinusoidal pulse shaping as defined in the standard
%It also implements a coherent matched filtering receiver and then computes
%the Packet Error Probability (PER) and the Bit Error probability (BER) for
%various packet sizes, packet numbers and SNR values.

            
                                                                        %PN_sequence table
                                                                          
                                                                         
Table=[1 1 0 1 1 0 0 1 1 1 0 0 0 0 1 1 0 1 0 1 0 0 1 0 0 0 1 0 1 1 1 0; 
       1 1 1 0 1 1 0 1 1 0 0 1 1 1 0 0 0 0 1 1 0 1 0 1 0 0 1 0 0 0 1 0;
       0 0 1 0 1 1 1 0 1 1 0 1 1 0 0 1 1 1 0 0 0 0 1 1 0 1 0 1 0 0 1 0;
       0 0 1 0 0 0 1 0 1 1 1 0 1 1 0 1 1 0 0 1 1 1 0 0 0 0 1 1 0 1 0 1;
       0 1 0 1 0 0 1 0 0 0 1 0 1 1 1 0 1 1 0 1 1 0 0 1 1 1 0 0 0 0 1 1;
       0 0 1 1 0 1 0 1 0 0 1 0 0 0 1 0 1 1 1 0 1 1 0 1 1 0 0 1 1 1 0 0;
       1 1 0 0 0 0 1 1 0 1 0 1 0 0 1 0 0 0 1 0 1 1 1 0 1 1 0 1 1 0 0 1;
       1 0 0 1 1 1 0 0 0 0 1 1 0 1 0 1 0 0 1 0 0 0 1 0 1 1 1 0 1 1 0 1;
       1 0 0 0 1 1 0 0 1 0 0 1 0 1 1 0 0 0 0 0 0 1 1 1 0 1 1 1 1 0 1 1;
       1 0 1 1 1 0 0 0 1 1 0 0 1 0 0 1 0 1 1 0 0 0 0 0 0 1 1 1 0 1 1 1;
       0 1 1 1 1 0 1 1 1 0 0 0 1 1 0 0 1 0 0 1 0 1 1 0 0 0 0 0 0 1 1 1;
       0 1 1 1 0 1 1 1 1 0 1 1 1 0 0 0 1 1 0 0 1 0 0 1 0 1 1 0 0 0 0 0;
       0 0 0 0 0 1 1 1 0 1 1 1 1 0 1 1 1 0 0 0 1 1 0 0 1 0 0 1 0 1 1 0;
       0 1 1 0 0 0 0 0 0 1 1 1 0 1 1 1 1 0 1 1 1 0 0 0 1 1 0 0 1 0 0 1;
       1 0 0 1 0 1 1 0 0 0 0 0 0 1 1 1 0 1 1 1 1 0 1 1 1 0 0 0 1 1 0 0;
       1 1 0 0 1 0 0 1 0 1 1 0 0 0 0 0 0 1 1 1 0 1 1 1 1 0 1 1 1 0 0 0];
   
re_table=transpose(Table); % transpose table for despreading
re_table(re_table==0)=-1;
 
quads=[0 0 0 0;1 0 0 0;0 1 0 0;1 1 0 0;0 0 1 0;1 0 1 0;0 1 1 0;1 1 1 0;0 0 0 1; % symbols table
      1 0 0 1;0 1 0 1;1 1 0 1;0 0 1 1;1 0 1 1;0 1 1 1;1 1 1 1];

N_bits=No_bytes*8; % bits in a packet                                                                        
data=rand(1,No_packets*N_bits)>0.5; % generate the data for all packets
packets=reshape(data,N_bits,No_packets); % each column is a packet
% No_bytes=1.25*10^5; % 10^6 bits
 
 for k=1:length(SNR)
     %msg=reshape(msg,4,[]); % 4 bits of a message consist one symbol                                                                      
     f=SNR(k) % mmonitor progress of SNR loop
                                                                           
 for   j=1:No_packets                                                                          
       packet=packets(:,j);
       msg=reshape(packet,4,[]); % reshape each packets in quads for dsss processing
       
       
                                                                           % Bit to symbol matching
                                                                                   
     for i=1:size(msg,2);
         symbol=msg(:,i);
         symbol=symbol';% from column to row

         symbol_flipped=fliplr(symbol); % in binary represtenation b0 is the last number on the left, whereas in the standard b0 is the first number. so in order matlab to convert to decimal i flip
   
     %Convert to decimal for pn-sequence matching'
         binary_symbol=num2str(symbol_flipped);
    
                                                                           % Symbol to chip matching
                                     
         decimal_symbol=bin2dec(binary_symbol);
    
         tx_sequence=Table(decimal_symbol+1,:);
    
 
                                                    % NRZ  pulses, I-Q branch construction and half-sinusoidal pulse shaping    
                                                                                                                             
         Tc=1/(2*10^6); % Chip period
         N_samples=4; % samples per bit for pulse shaping
         bipolar_sequence=2*tx_sequence-1; % bit 0 -> -1, bit 1 -> 1 
         IBranch=bipolar_sequence(1:2:end); % In-phase component
         QBranch=bipolar_sequence(2:2:end); % Quadrature component
         t=Tc/4:Tc/4:2*Tc; % 4 samples per bit period, thus 8 samples per 2-bit periods
         p=sin(pi*t/(2*Tc)); % half-sine pulse shape
         shaped_IBranch=kron(IBranch,p);
         shaped_QBranch=kron(QBranch,p);
                                                                 % O-QPSK modulation (Offset QPSK)/ transmitted signal 
                                                       
         final_IBranch=[shaped_IBranch zeros(1,N_samples)]; % add zeros to have the same dimensions with the Q-branch
         final_QBranch=[zeros(1,N_samples) shaped_QBranch]; % delay addition for O-QPSK modulation
         tx_waveform =(final_IBranch + 1i*final_QBranch); 
      
                                                                         % Additive white Gaussian Noise
                                                                         
         rand('state', sum(100*clock));                                                                  
         S=10*log10(mean(abs(tx_waveform).^2)); % signal power  in dB                                         
         noise = 1/sqrt(2)*[randn(1,length(tx_waveform)) + 1i*randn(1,length(tx_waveform))]; % white gaussian noise samples of 0dB variance 
         noise_power=S-SNR(k); % noise variance-dB (noise power)
         noise_volt=10^(noise_power/20);
         tx_waveform=tx_waveform +noise_volt*noise; % noisy signal
      
                                    
                                                                   % Demodulation (half-sinusoidal matched filter)       
                                                                                                                                   
         RxIbranch= real(tx_waveform);
         RxQbranch=imag(tx_waveform);
         odd_bits=conv(RxIbranch,p); % I-Branch matched filter
         even_bits=conv(RxQbranch,p); % Q-Branch matched filter
         
                                                                       % Peak-Sampling and decision function
                                                                                
        rx_bits=zeros(1,length(tx_sequence));
        rx_bits(1:2:end)=odd_bits(2*N_samples:2*N_samples:end-2*N_samples);
        rx_bits(2:2:end)=even_bits(3*N_samples:2*N_samples:end-N_samples);
        find_ones=find(rx_bits>0);
        find_zeros=find(rx_bits<0);
        rx_bits(find_ones)=1;
        rx_bits(find_zeros)=-1;
     
                                                                     % Despreading and symbol extraction
                                                                            
        despread=rx_bits*re_table; % despread (multiplication and integration)
        index=find(despread==max(despread));
        index_chosen=index(1);
        symbol_rx(:,i)=quads(index_chosen,:); % LUT among chips and symbols (extract transmitted symbol)
      
    end % end symbols
 
 packet_tx=msg(:)';
 packet_rx=symbol_rx(:)';
 pe(j)=1-isequal(packet_tx,packet_rx); % if equal, no error so  pe(j)=0, if not pe(j)=1
 nErr(j)=size(find(packet_tx-packet_rx),2); % bits in error within a packet
 
 end % end packets
                                                              
                                                                         % PER and BER calculation
                                                                       
 PER(k)=sum(pe)/No_packets; 
 BER(k)=sum(nErr)/(No_packets*N_bits);
                                                                              
 end % end SNR         
      
    
    
    

    
    
    
                                                                    

