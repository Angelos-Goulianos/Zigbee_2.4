

% This is a Zigbee fh end to end transmitter receiver which implements a
% sub-optimum  receiver that  compares the received signal to the OQPSK
% reference constellation, according to MATLAB comm. toolbox objects.


Table=[1 1 0 1 1 0 0 1 1 1 0 0 0 0 1 1 0 1 0 1 0 0 1 0 0 0 1 0 1 1 1 0;     % PN sequence
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
  
hmod=comm.OQPSKModulator('BitInput',true); % construct OQPSK mod object
hDemod=comm.OQPSKDemodulator('BitOutput',true);% construct OQPSK demod object

No_bytes=No_bytes+0.5; % add  1 more symbol due to oqpsk matlab algorithm: to be removed in the end so the total bits are No_bytes
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
         tx_sequence=tx_sequence.';
                                      
                                                                 % O-QPSK modulation (Offset QPSK)/ transmitted signal 
                                                                 
         modData=step(hmod,tx_sequence); % modulate data                                                        
                                                       
         
                                                                         % Additive white Gaussian Noise
                                                                         
         rand('state', sum(100*clock));                                                                  
         S=10*log10(mean(abs(modData).^2)); % signal power  in dB                                         
         noise = 1/sqrt(2)*[randn(1,length(modData)) + 1i*randn(1,length(modData))]; % white gaussian noise samples of 0dB variance 
         noise_power=S-SNR(k); % noise variance-dB (noise power)
         noise_volt=10^(noise_power/20);
         tx_waveform=modData +noise_volt*noise.'; % noisy signal
      
                                    
                                                                                 % Demodulation      
        rx_bits=step(hDemod,tx_waveform);
        rx_bits=circshift(rx_bits,[-2 0]); % 2 bits lagging due to the matlab oqpsk demod. algorithm                                                                                                                           
         
     
                                                                        % Despreading and symbol extraction
        rx_seq=rx_bits.';
        rx_seq(rx_seq==0)=-1;
        despread=rx_seq*re_table; % despread (multiplication and integration)
        index=find(despread==max(despread));
        index_chosen=index(1);
        symbol_rx(:,i)=quads(index_chosen,:); % LUT among chips and symbols (extract transmitted symbol)
      
    end % end symbols
 
 packet_tx=msg(:)';
 packet_tx=packet_tx(1:end-4); %the last symbol is not compared due to matlab oqpsk alg.
 packet_rx=symbol_rx(:)';
 packet_rx=packet_rx(1:end-4); %same at the receiver which has been shifted before and the first bits have become last and do not count
 pe(j)=1-isequal(packet_tx,packet_rx); % if equal, no error so  pe(j)=0, if not pe(j)=1
 nErr(j)=size(find(packet_tx-packet_rx),2); % bits in error within a packet
 
 end % end packets
                                                                             % PER and BER calculation                                                                      
 PER(k)=sum(pe)/No_packets; 
 BER(k)=sum(nErr)/(No_packets*N_bits);
                                                                              
 end % end SNR         
