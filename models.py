
condlif_sig = '''
              dV/dt = (El-V + ge*(Ee-V) + gi*(Ei-V))/tau 
                      + sigma * xi / (tau **.5): volt
              Vt : volt 
              dge /dt = -ge/tau_e : 1
              dgi /dt = -gi/tau_i : 1

              Asum : 1
              
              spk_count : 1
              h_ip : Hz (constant)
              '''

# refractory period???

nrnEE_thrshld = 'V > Vt'

nrnEE_reset = '''
              V = Vr_e
              spk_count = spk_count + 1
              '''

synEE_mod = '''
            a : 1
            active : integer

            dApre  /dt = -Apre/taupre  : 1
            dApost /dt = -Apost/taupost : 1

            Asum_post = a : 1 (summed)         
            insert_P : 1 (shared) 
            '''

# synEE_mod = '''
#             a : 1
#             dApre  /dt = -Apre/taupre  : 1 (event-driven)
#             dApost /dt = -Apost/taupost : 1 (event-driven)
#             '''

synEE_pre = '''
            ge_post += a
            Apre += Aplus
            a = clip(a+Apost, 0, amax)
            '''

synEE_post = '''
             Apost+= Aminus
             a = clip(a+Apre, 0, amax)
             '''


synEE_scaling = '''
                a = a*(ATotalMax/Asum_post)
                '''

intrinsic_mod = '''
                Vt = Vt + eta_ip*(spk_count/it_dt - h_ip)
                spk_count = 0
                '''

# rand() == uniform(0,1)
#strct_mod = ''
strct_mod = '''
            r = rand()
            should_stay_active = (a > prn_thrshld)
            should_become_active = (r < insert_P)
            was_active_before = active
            active = int(active==1) * int(should_stay_active) \
                     + int(active==0) * int(should_become_active)
            a = a*int(was_active_before==1)*int(active==1) \
                + a_insert*int(was_active_before==0)*int(active==1)
             '''
