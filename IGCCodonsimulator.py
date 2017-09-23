# -*- coding: utf-8 -*-   

#This is my final course project of CSC530.
#It is a simulation study of interlocus gene conversion (IGC)
# between paralogs together with point mutation process during
# evolution on a single branch. Simulated shotgun sequencing reads are 
#then simulated using the simulated sequence and then re-assembled. 
#It aims at understanding the influence of IGC in sequencing identifiability of the paralogs. 
#Xiang Ji
#xji3@ncsu.edu

from CodonGeneconFunc import *
from copy import deepcopy

class OneBranchIGCCodonSimulator:
    def __init__(self, blen, num_exon,
                 x_exon, x_IGC,  log_file, div_file,
                 initial_seq = None, 
                 intron_model = 'HKY', exon_model = 'MG94',
                 has_intron = True, num_paralog = 2):
        
        
        self.blen         = blen          # branch length
        self.num_exon     = num_exon      # number of site in the exon-intron-exon sequence exonic region (same at both ends)
        self.x_exon       = x_exon        # parameter values for exon model
        self.x_IGC        = x_IGC         # parameter values for IGC process
                                          # IGC_init_G, IGC_geo_q, IGC_threshold
        self.log_file     = log_file      # log file position
        self.div_file     = div_file      # divergence file position

        self.has_intron   = has_intron    # whether consider intronic region in the simulation, default value is True
        self.intron_model = intron_model  # intronic region substitution model
        self.exon_model   = exon_model    # exonic region substitution model
        self.num_paralog  = num_paralog   # number of paralogs in the simulation, currently 2
        self.initial_seq  = initial_seq   # initial sequences of paralogs, will simulate if not given
        self.current_seq  = initial_seq
        # list of [paralog site state list]
        # paralog site state list = [[exon 1], [intron], [exon 2]]

        IGC_parameters     = self.unpack_x_IGC()
        self.IGC_init_g    = IGC_parameters[0]
        self.IGC_geo_q     = IGC_parameters[1]
        self.IGC_threshold = IGC_parameters[2]  # IGC tract similarity threshold control

        self.point_mut_count   = 0        # total number of point mutation in the simulation
        self.IGC_total_count   = 0        # total number of IGC events in the simulation
        self.IGC_failure_count = 0        # total number of failed IGC events in the simulation
        self.IGC_contribution_count = 0   # total number of sites experiencing IGC
        self.IGC_tract_length  = 0        # used to track each IGC event tract length
        self.IGC_change_sites  = 0        # used to track changes due to IGC event

        self.exon_mut_Q   = None
        self.intron_mut_Q = None

        self.exon_distn   = None
        self.intron_distn = None

        self.mut_poisson_rate = None  # list of lists of poisson rates in exons and introns
        self.mut_total_rate   = None  # total rate of point mutation
        self.IGC_total_rate   = None  # total rate of IGC 
        self.IGC_J            = None  # IGC "transition" matrix

        # Constants for Sequence operations
        bases = 'tcag'.upper()
        codons = [a+b+c for a in bases for b in bases for c in bases]
        amino_acids = 'FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG'
        
        self.nt_to_state    = {a:i for (i, a) in enumerate('ACGT')}
        self.state_to_nt    = {i:a for (i, a) in enumerate('ACGT')}
        self.codon_table    = dict(zip(codons, amino_acids))
        self.codon_nonstop  = [a for a in self.codon_table.keys() if not self.codon_table[a]=='*']
        self.codon_to_state = {a.upper() : i for (i, a) in enumerate(self.codon_nonstop)}
        self.state_to_codon = {i:a.upper() for (i, a) in enumerate(self.codon_nonstop)}

        self.initialize()


    def initialize(self):
        if self.exon_model == 'MG94':
            # get stationary nucleotide distribution of codon model of exonic region
            pi_codon = self.unpack_x_exon()[0]
            distn_codon = [ reduce(mul, [pi_codon['ACGT'.index(b)]  for b in codon], 1) for codon in self.codon_nonstop ]#codon 先验
            distn_codon = np.array(distn_codon) / sum(distn_codon)#codon 先验分布
            self.exon_distn = distn_codon
            self.exon_mut_Q = self.get_MG94()

        self.get_IGC_rate()
        if self.initial_seq == None:
            self.initial_seq = self.sim_initial_seq()
        self.current_seq = self.convert_seq_to_list(self.initial_seq)
        self.get_mutation_rate()
        with open(self.log_file, 'w+') as f:
            f.write('Initial seq: \n' + '\n'.join(self.convert_list_to_seq()) + '\n')
            f.write('Exon model: ' + self.exon_model + ' Exon nsite: ' + str(self.num_exon) + ' blen: ' + str(self.blen) + '\n')
            f.write('x_exon: ' + ', '.join([str(parameter) for parameter in self.x_exon]) + '\n')
            f.write('x_IGC: ' + ', '.join([str(parameter) for parameter in self.x_IGC]) + '\n')

        with open(self.div_file, 'w+') as f:
            f.write('# time \t divergence \t # of mutations \t # IGC chances \t IGC length \t # IGC failures \t Changes due to IGC  \t % change due to IGC \n')

        
    def sim_initial_seq(self):
        sim_seq = []
        exon_1 = draw_from_distribution(self.exon_distn, self.num_exon, self.codon_nonstop)
            
        for i in range(self.num_paralog):
            seq = deepcopy(exon_1)
            sim_seq.append(''.join(seq))
        return sim_seq
        

    def unpack_x_IGC(self, log = False):
        if log:
            x_IGC = np.exp(self.x_IGC)
        else:
            x_IGC = self.x_IGC
        assert(len(x_IGC) == 3)
        return x_IGC

    def unpack_x_exon(self, log = False):
        if log:
            x_exon = np.exp(self.x_exon)
        else:
            x_exon = self.x_exon
            
        if self.exon_model == 'MG94':
            assert(len(self.x_exon) == 5)
            # %AG, %A, %C, kappa, omega
            pi_a = x_exon[0] * x_exon[1]
            pi_c = (1 - x_exon[0]) * x_exon[2]
            pi_g = x_exon[0] * (1 - x_exon[1])
            pi_t = (1 - x_exon[0]) * (1 - x_exon[2])
            pi = [pi_a, pi_c, pi_g, pi_t]

            kappa = x_exon[3]
            omega = x_exon[4]
            return [pi, kappa, omega]


    def get_MG94(self):
        Qbasic = np.zeros((61, 61), dtype = float)
        pi, kappa, omega = self.unpack_x_exon()
        for ca in self.codon_nonstop:
            for cb in self.codon_nonstop:
                if ca == cb:
                    continue
                Qbasic[self.codon_to_state[ca], self.codon_to_state[cb]] = get_MG94BasicRate(ca, cb, pi, kappa, omega, self.codon_table)
        expected_rate = np.dot(self.exon_distn, Qbasic.sum(axis = 1))
        Qbasic = Qbasic / expected_rate
        return Qbasic

    def get_HKY(self):
        pi, kappa = self.unpack_x_intron()
        Qbasic = np.array([
            [0, 1.0, kappa, 1.0],
            [1.0, 0, 1.0, kappa],
            [kappa, 1.0, 0, 1.0],
            [1.0, kappa, 1.0, 0],
            ]) * np.array(pi)
        expected_rate = np.dot(self.intron_distn, Qbasic.sum(axis = 1))
        Qbasic = Qbasic / expected_rate
        return Qbasic

    def get_mutation_rate(self):
        #print self.current_seq
        poisson_rate_sum = 0.0
        exon_diag_rates = self.exon_mut_Q.sum(axis = 1)
#        intron_diag_rates = self.intron_mut_Q.sum(axis = 1)

        seq_rate_list = []
        for i in range(self.num_paralog):
            seq_rate = []
            exon_rate = [exon_diag_rates[self.codon_to_state[self.current_seq[i][j]]] for j in range(self.num_exon)]
            seq_rate.extend(exon_rate)
            seq_rate_list.extend(seq_rate)

        self.mut_poisson_rate = seq_rate_list
        self.mut_total_rate = sum(seq_rate_list)

    def get_IGC_rate(self):
        # Now prepare "rate matrix" for the gene conversion event
        # row number is the starting position, column number is the ending position
        # This matrix only need to prepare once and then stay constant
        if self.exon_model == 'MG94':
            S = self.num_exon
        J = np.zeros((S, S))  # copying my formula part from Alex's script
        l = np.arange(S+1, dtype=float)
        p = (self.IGC_geo_q) * (1 - self.IGC_geo_q)**l  # geometric distribution
        for i in range(S):
            for j in range(i, S):
                l = j - i  # This is the number of failures in Geometric distribution
                # start and end inside the gene with length l <= S - 2
                if 0 < i <= j < S - 1:
                    J[i, j] += p[l]
                # start before gene, end within gene with length l <= S - 1
                if 0 == i <= j < S - 1:
                    J[i, j] += p[l] * 1.0 / self.IGC_geo_q
                # start within gene, end after gene with length l <= S - 1
                if 0 < i <= j == S - 1:
                    J[i, j] += p[l] * 1.0 / self.IGC_geo_q
                # start before gene, end after gene with length l = S
                if 0 == i <= j == S - 1:
                    J[i, j] += p[l] * 1.0 / self.IGC_geo_q ** 2                
        # Now assign initiation rate using tau parameter
        init_rate = self.IGC_init_g
        J *= init_rate
        self.IGC_total_rate = J.sum()
        self.IGC_J = J

    def sim_one_branch(self, starting_seq, blen):
        current_seq = deepcopy(starting_seq)
        cummulate_time = 0.0

        self.point_mut_count   = 0
        self.IGC_total_count   = 0
        self.IGC_failure_count = 0
        self.IGC_contribution_count = 0

        while(cummulate_time < blen):
            # Now sample exponential distributed waiting time for next event
            # point mutation or IGC event
            ## not only recording divergence this time
            
            div_info = [str(cummulate_time), str(1.0 - self.compare_similarity()), str(self.point_mut_count),
                        str(self.IGC_total_count), str(self.IGC_tract_length), str(self.IGC_failure_count), str(self.IGC_change_sites)]
            if self.point_mut_count + self.IGC_change_sites == 0:
                div_info.append(str(0.0))
            else:
                div_info.append(str((self.IGC_change_sites + 0.0) / (self.IGC_change_sites + self.point_mut_count + 0.0)))
                
            self.add_div(div_info)
            total_rate = self.mut_total_rate + 2 * self.IGC_total_rate#第一项所有位点rate之和，第二项所有位点IGCrate之和
            if total_rate == 0.0:
                break
            cummulate_time += np.random.exponential(1.0 / total_rate)

            if cummulate_time > blen:
                break
            else:
                event = draw_from_distribution(np.array([self.mut_total_rate, self.IGC_total_rate, self.IGC_total_rate]) / total_rate,
                                               1, range(3))
                if event == 0:
                    # It's a point mutation event
                    self.get_one_point_mutation()
                elif event == 1:
                    # It's an IGC event
                    # copy from paralog 0 to paralog 1
                    self.get_one_IGC_event(template_paralog = 0, target_paralog = 1)
                elif event == 2:
                    # It's an IGC event
                    # copy from paralog 1 to paralog 0
                    self.get_one_IGC_event(template_paralog = 1, target_paralog = 0)

        # intermidiate step, don't need to add final sequence
        # self.add_final_seq()

    def add_final_seq(self):
        with open(self.log_file, 'a') as f:
            f.write('Final seq: \n')
            f.write('\n'.join(self.convert_list_to_seq()) + '\n')
            f.write(str(self.compare_similarity()) + '\n')


    def get_one_point_mutation(self):
        mut_pos = draw_from_distribution(self.mut_poisson_rate / sum(self.mut_poisson_rate), 1, range(len(self.mut_poisson_rate)))
        #print len(self.mut_poisson_rate)
        # mut_pos is where point mutation happens
        mut_paralog = mut_pos / (self.num_exon)
        
        mut_pos = mut_pos - mut_paralog * (self.num_exon)

        exon_pos = mut_pos
        old_seq = self.current_seq[mut_paralog][exon_pos]
        prob = np.array(self.exon_mut_Q[self.codon_to_state[self.current_seq[mut_paralog][exon_pos]],:])
        new_state = draw_from_distribution(prob / sum(prob), 1, range(len(prob)))
        self.current_seq[mut_paralog][exon_pos] = self.codon_nonstop[new_state]
        new_seq = self.current_seq[mut_paralog][exon_pos]

        self.point_mut_count += 1

        # Now updaate mutation rate
        self.get_mutation_rate()
        mut_info = ['Point mutation at paralog ' + str(mut_paralog) + ' position ' + str(mut_pos) + ' from', old_seq, 'to', new_seq, str(self.compare_similarity())]
        self.add_log(mut_info)
##        print('  '.join(self.convert_list_to_seq()) + '  ' + str(self.compare_similarity()))

    def get_one_IGC_event(self, template_paralog, target_paralog):
        IGC_J_reshape = np.array(self.IGC_J).reshape(-1)
        IGC_pos = draw_from_distribution(IGC_J_reshape / sum(IGC_J_reshape), 1, range(len(IGC_J_reshape)))
        IGC_start = IGC_pos / self.IGC_J.shape[1]   # now in unit as codon 
        IGC_stop  = IGC_pos - self.IGC_J.shape[1] * IGC_start   # now in unit as codon

        # print IGC_pos, IGC_start, IGC_stop, template_paralog, target_paralog
        self.IGC_copy(IGC_start, IGC_stop, template_paralog, target_paralog)
        # Now updaate mutation rate
        self.get_mutation_rate()
        
    def IGC_copy(self, start, stop, template_paralog, target_paralog):  # now start and stop positions are in unit as codon not nucleotide
        if self.exon_model == 'MG94':
            exon_seq_len = self.num_exon

        nuc_start = 3 * start  # transform from codon into nucleotide as unit
        nuc_stop  = 3 * (stop + 1)

        paralog_seq_list = self.convert_list_to_seq()
        # now compare similarity between the two tracts
        tract_1 = paralog_seq_list[template_paralog][nuc_start:nuc_stop]
        tract_2 = paralog_seq_list[target_paralog][nuc_start:nuc_stop]
        compare_result = [tract_1[i] == tract_2[i] for i in range(len(tract_1))]
        
        self.IGC_change_sites += nuc_stop - nuc_start - sum(compare_result)  # New, track changes due to IGC rather than point mutation
        
        similarity = (sum(compare_result) + 0.0) / (len(tract_1) + 0.0)
        if similarity > self.IGC_threshold:
            paralog_seq_list[target_paralog] = paralog_seq_list[target_paralog][:nuc_start] + paralog_seq_list[template_paralog][nuc_start:nuc_stop] + paralog_seq_list[target_paralog][nuc_stop:]

            new_paralog_seq = self.convert_seq_to_list(paralog_seq_list)#            [self.convert_seq_to_list(paralog_seq) for paralog_seq in paralog_seq_list]
            self.current_seq = new_paralog_seq
        else:
            self.IGC_failure_count += 1#IGC failure due to trying to copy from a bad model
        self.IGC_total_count += 1
        self.IGC_contribution_count += (nuc_stop - nuc_start) * 2
        self.IGC_tract_length = (nuc_stop - nuc_start)

        # move add_log into here        
        IGC_info = ['IGC ' + str(template_paralog) + ' -> ' + str(target_paralog) + ' : start', str(nuc_start + 1), 'stop', str(nuc_stop), str(self.compare_similarity())
                    + ' from ' + tract_2 + ' -> ' + tract_1]
        self.add_log(IGC_info)

        
    def convert_list_to_seq(self):
        paralog_seq_list = []
        for i in range(self.num_paralog):
            paralog_seq_list.append(''.join(self.current_seq[i]))

        return (deepcopy(paralog_seq_list))
        
    def convert_seq_to_list(self, seq):
        exon_1 = [seq[0][(3*event):(3*event + 3)] for event in range(self.num_exon)]
        exon_2 = [seq[1][(3*event):(3*event + 3)] for event in range(self.num_exon)]
        return([exon_1, exon_2])        
        
    def compare_similarity(self):
        assert(self.num_paralog == 2)
        paralog_seq_list = self.convert_list_to_seq()
        assert(len(paralog_seq_list) == 2)
        assert(len(paralog_seq_list[0]) == len(paralog_seq_list[1]))
        compare_result = [paralog_seq_list[0][i] == paralog_seq_list[1][i] for i in range(len(paralog_seq_list[0]))]
        similarity = (sum(compare_result) + 0.0) / (len(paralog_seq_list[0]) + 0.0)
        return similarity

    def add_log(self, info):
        with open(self.log_file, 'a') as f:
            f.write('\t'.join(info) + '\n')
            #f.write('\n'.join(self.convert_list_to_seq()) + '\n')

    def add_div(self, info):
        with open(self.div_file, 'a') as f:
            f.write('\t'.join(info) + '\n')

def draw_from_distribution(prob, size, values):
    prob = np.array(prob) / sum(prob)
    bins = np.add.accumulate(prob)
    if size == 1:
        return values[np.digitize(np.random.random_sample(size), bins)[0]]
    else:
        return [values[i] for i in np.digitize(np.random.random_sample(size), bins)]
 

            
if __name__ == '__main__':
    blen = 2.5
    #blen = 200.0
    num_exon = 100
    IGC_g = 0.07
    IGC_q = 0.025
    IGC_threshold = 0.8

    IGC_geo = 3.0
    tau = 1.409408
    IGC_geo_codon = IGC_geo / 3.0
    IGC_init = tau / IGC_geo_codon
    x_IGC = [IGC_init, 1.0 / IGC_geo_codon, IGC_threshold]  # These values vary for the simulation study

    
    # x_exon and x_intron from MG94_YBR191W_YPL079W_nonclock_save.txt
    x_exon = [0.49249355302375575, 0.60985035555249456, 0.42155795722934408, 8.1662933909645563, 0.092804167727196338]
#    x_IGC = [IGC_g, IGC_q, IGC_threshold]  # These values vary for the simulation study
    
    pi = [0.3003473684188488, 0.21394338106802796, 0.19214618460490698, 0.2935630659082163]
    div_limit = 1 - sum([i **2 for i in pi])
    for replicate in range(0, 1):
        log_file = './logs/3rd_log_g_' + str(IGC_g) + '_q_' + str(IGC_q) + '_threshold_' + str(IGC_threshold) + '_rep_' + str(replicate) + '.log'
        div_file = './logs/3rd_div_g_' + str(IGC_g) + '_q_' + str(IGC_q) + '_threshold_' + str(IGC_threshold) + '_rep_' + str(replicate) + '.log'
        test = OneBranchIGCCodonSimulator(blen, num_exon, x_exon, x_IGC, log_file, div_file)

        self = test

    test.sim_one_branch(test.initial_seq, blen)
    
##        try:
##            test.sim_one_branch(test.initial_seq, blen)
##        except:
##            print "Failed"
##            test.add_final_seq()
