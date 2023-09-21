function [p, state_AB] = generate_fm0_sig2(tx_code, fm0_samps, state_AB)
    clk = 1;
    
    if nargin < 3
        state_AB = [1, 1];
    end
    
    p = [];
    
    samps = floor(fm0_samps);
    resid = fm0_samps - samps;
    resid_tot = 0;
    flag = 0;
    
    for k=1:length(tx_code)
       for m = 1:2
          if clk == 1
             state_AB = [(~state_AB(2)) & 1, xor(tx_code(k), (state_AB(2) & 1)) & 1];
          end
          gen_sig = ((clk&1) * state_AB(1) + ((~clk)&1)*state_AB(2)) & 1;
          % make DC center at 0
          gen_sig = (gen_sig - 0.5);
          if (m == 1) % deal with odd size sample
             
              if resid_tot >=0.999
                  samps = samps+1;
                  resid_tot = resid_tot-1;
                  flag = 1;
              else
                  resid_tot = resid_tot + resid;
              end
              
              p = [p, gen_sig*ones(1, floor(samps/2))];
              
              if flag == 1
                samps = samps-1;
                flag = 0;
              end
          else
              if resid_tot >=0.999
                  samps = samps+1;
                  resid_tot = resid_tot-1;
                  flag = 1;
              else
                  resid_tot = resid_tot + resid;
              end
              p = [p, gen_sig*ones(1, ceil(samps/2))]; 
              
              if flag == 1
                samps = samps-1;
                flag = 0;
              end
              
          end
          clk = ~clk;
       end
       
    end

    p = 2*p;
end