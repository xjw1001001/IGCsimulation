# IGC Simulation on a Tree
# Imports simulation on a branch from my CSC530 project with modification
# Xiang Ji
# xji3@ncsu.edu

# I know I really shouldn't do this, but forgive me for being lazy...


# Unlike the one branch simulator, this one is quite simple
# consider only one exon with codon model

from IGCCodonsimulator import OneBranchIGCCodonSimulator, draw_from_distribution
from CodonGeneconFunc import *
import cPickle

class TreeIGCCodonSimulator:
    def __init__(self, num_exon, newick_tree, paralog, seq_file, log_file,
                 x_exon, x_IGC, log_folder, div_folder, seed_number):
        self.newicktree   = newick_tree
        self.num_exon     = num_exon
        self.pair         = paralog
        self.seq_file     = seq_file
        self.log_file     = log_file
        self.seed_number  = seed_number 
        self.edge_to_blen = None
        self.node_to_num  = None
        self.num_to_node  = None
        self.edge_list    = None
        self.outgroup     = [('N0', 'kluyveri')]  # I am so lazy

        # Node sequence
        self.node_to_sequence = None
        self.node_to_sim      = {}

        # OneBranchIGCSimulator Related Parameters
        self.x_exon       = x_exon        # parameter values for exon model
        self.x_IGC        = x_IGC
        self.Model   = 'MG94'
        self.num_paralog  = 2
        self.log_folder   = log_folder
        self.div_folder   = div_folder

        self.distn         = None
        self.mut_Q         = None

        self.OneBranchSimulator = None

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
        if self.Model == 'MG94':
            self.pair_to_state  = {pair:i for i, pair in enumerate(product(self.codon_nonstop, repeat = 2))}
        self.state_to_pair  = {self.pair_to_state[pair]:pair for pair in self.pair_to_state}

        # Total event counts
        self.total_mut         = 0  # Total mutation events
        self.total_IGC         = 0  # Total IGC events
        self.total_IGC_sites   = 0  # Total IGC affecting sites
        self.total_IGC_changes = 0  # Total changes due to IGC rather than point mutation


        self.initiate()
        


    def initiate(self):
        self.get_tree()
        self.unpack_x()
        self.set_seed()
##        # If the seed file exists, set numpy's random seed state according to the seed file
##        seed_file = self.log_file.replace('.log', '_seed.log')
##        if os.path.isfile(seed_file):
##            prng = cPickle.load(open(seed_file, 'r'))
##            np.random.set_state(prng.get_state())
##        else:
##            prng = np.random.RandomState()
##            seed_file = self.log_file.replace('.log', '_seed.log')
##            cPickle.dump(prng, open(seed_file, 'w+'))

    def set_seed(self):
        assert(type(self.seed_number) == int)
        np.random.seed(self.seed_number)

    def unpack_x(self):
        if self.Model == 'MG94':
            # get stationary nucleotide distribution of codon model of exonic region
            pi_codon = self.unpack_x_exon()[0]
            distn_codon = [ reduce(mul, [pi_codon['ACGT'.index(b)]  for b in codon], 1) for codon in self.codon_nonstop ]
            distn_codon = np.array(distn_codon) / sum(distn_codon)
            self.distn = distn_codon
            self.mut_Q = self.get_MG94()
            
        self.node_to_sequence = {node:[] for node in self.node_to_num.keys()}

    def unpack_x_exon(self, log = False):
        if log:
            x_exon = np.exp(self.x_exon)
        else:
            x_exon = self.x_exon
            
        if self.Model == 'MG94':
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
        expected_rate = np.dot(self.distn, Qbasic.sum(axis = 1))
        Qbasic = Qbasic / expected_rate
        return Qbasic
    
    def get_tree(self):
        tree = Phylo.read( self.newicktree, "newick")
        #set node number for nonterminal nodes and specify root node
        numInternalNode = 0
        for clade in tree.get_nonterminals():
            clade.name = 'N' + str(numInternalNode)
            numInternalNode += 1
        tree_phy = tree.as_phyloxml(rooted = 'True')
        tree_nx = Phylo.to_networkx(tree_phy)

        triples = ((u.name, v.name, d['weight']) for (u, v, d) in tree_nx.edges(data = True)) # data = True to have the blen as 'weight'
        T = nx.DiGraph()
        edge_to_blen = {}
        for va, vb, blen in triples:
            edge = (va, vb)
            T.add_edge(*edge)
            edge_to_blen[edge] = blen

        self.edge_to_blen = edge_to_blen

        # Now assign node_to_num
        leaves = set(v for v, degree in T.degree().items() if degree == 1)
        self.leaves = list(leaves)
        internal_nodes = set(list(T)).difference(leaves)
        node_names = list(internal_nodes) + list(leaves)
        self.node_to_num = {n:i for i, n in enumerate(node_names)}
        self.num_to_node = {self.node_to_num[i]:i for i in self.node_to_num}

        # Prepare for generating self.tree so that it has same order as the self.x_process
        nEdge = len(self.edge_to_blen)  # number of edges
        l = nEdge / 2 + 1               # number of leaves
        k = l - 1   # number of internal nodes. The notation here is inconsistent with Alex's for trying to match my notes.

        leaf_branch = [edge for edge in self.edge_to_blen.keys() if edge[0][0] == 'N' and str.isdigit(edge[0][1:]) and not str.isdigit(edge[1][1:])]
        out_group_branch = [edge for edge in leaf_branch if edge[0] == 'N0' and not str.isdigit(edge[1][1:])] [0]
        internal_branch = [x for x in self.edge_to_blen.keys() if not x in leaf_branch]
        assert(len(internal_branch) == k-1)  # check if number of internal branch is one less than number of internal nodes

        leaf_branch.sort(key = lambda node: int(node[0][1:]))  # sort the list by the first node number in increasing order
        internal_branch.sort(key = lambda node: int(node[0][1:]))  # sort the list by the first node number in increasing order
        edge_list = []
        for i in range(len(internal_branch)):
            edge_list.append(internal_branch[i])
            edge_list.append(leaf_branch[i])
        for j in range(len(leaf_branch[i + 1:])):
            edge_list.append(leaf_branch[i + 1 + j])

        self.edge_list = edge_list

    def unpack_x_rates(self, x_rates):
        for edge_iter in range(len(self.edge_list)):
            edge = self.edge_list[edge_iter]
            self.edge_to_blen[edge] = x_rates[edge_iter]

    def sim(self):
        self.sim_root()
        for edge in self.edge_list:
            
            #print edge, self.edge_to_blen[edge]

            # Now need to adapt branchsim
            # create an instance for each branch
            blen = self.edge_to_blen[edge]
            num_exon = self.num_exon
            x_exon = self.x_exon
            x_IGC  = deepcopy(self.x_IGC)
            log_file = self.log_folder + '_'.join(edge) + '_log.log'
            div_file = self.div_folder + '_'.join(edge) + '_div.log'
            starting_seq = self.node_to_sequence[edge[0]]

            if edge in self.outgroup:
                x_IGC[0] = 0.0
                        
            self.OneBranchSimulator = OneBranchIGCCodonSimulator(blen = blen, num_exon = num_exon,
                                               x_exon = x_exon, x_IGC = x_IGC,
                                               log_file = log_file, div_file = div_file, initial_seq = starting_seq)
            
            blen = self.edge_to_blen[edge] 

            self.OneBranchSimulator.sim_one_branch(starting_seq, blen)
            self.node_to_sequence[edge[1]] = self.OneBranchSimulator.convert_list_to_seq()

            self.total_mut += self.OneBranchSimulator.point_mut_count
            self.total_IGC += self.OneBranchSimulator.IGC_total_count
            self.total_IGC_sites += self.OneBranchSimulator.IGC_contribution_count
            self.total_IGC_changes += self.OneBranchSimulator.IGC_change_sites

        self.output_seq()
        self.get_log()
        self.write_log()




    def sim_root(self):
        if self.Model == 'MG94':
            seq = draw_from_distribution(self.distn, self.num_exon, self.codon_nonstop)

        #self.node_to_sequence['N0'] = np.array([self.pair_to_state[(i, i)] for i in seq])
        self.node_to_sequence['N0'] = [''.join(seq), ''.join(seq)]
        self.node_to_sim['N0'] = [self.node_to_sequence['N0'], 0, 0]

    def output_seq(self):
        with open(self.seq_file, 'w+') as f:
            for node in self.node_to_sequence:
                if not node in self.outgroup[0]:
                    for paralog_counter in range(self.num_paralog):
                        paralog = self.pair[paralog_counter]
                        f.write('>' + node + paralog + '\n')
                        f.write(self.node_to_sequence[node][paralog_counter] + '\n')
                else:  # only observe one paralog of the outgroup species
                    paralog = self.pair[0]
                    f.write('>' + node + paralog + '\n')
                    f.write(self.node_to_sequence[node][paralog_counter] + '\n')        

    def get_log(self):
        with open(self.log_file, 'w+') as f:
            f.write('Model: ' + self.Model + '  nSites: ' + str(self.num_exon) + ' TreeFile: ' + self.newicktree + '\n')
            f.write('x_exon: ' + ', '.join([str(parameter) for parameter in self.x_exon]) + '\n')
            f.write('x_IGC: ' + ', '.join([str(parameter) for parameter in self.x_IGC]) + '\n')
            f.write('\n'.join(['_'.join(edge) + ': ' + str(self.edge_to_blen[edge]) for edge in self.edge_list]) + '\n')
            f.write('Total point mutation: ' + str(self.total_mut) + '  Total IGC: ' + str(self.total_IGC)
                    + '  Total IGC affecting sites: ' + str(self.total_IGC_sites)
                    + ' Total changes due to IGC: ' + str(self.total_IGC_changes) + '\n')
            f.write('% change due to IGC: ' + str((self.total_IGC_changes + 0.0) / (self.total_mut + self.total_IGC_changes + 0.0)) + '\n')
            
    def write_log(self):

        label = ['%IGC', '#mut', '#IGC']
        summary = np.matrix([(self.total_IGC_changes + 0.0) / (self.total_mut + self.total_IGC_changes + 0.0),
                   self.total_mut, self.total_IGC])
        short_log_file = self.log_file.replace('.log', '_short.log')
        footer = ' '.join(label)
        np.savetxt(open(short_log_file, 'w+'), summary.T, delimiter = ' ', footer = footer)


if __name__ == '__main__':
    paralog1 = 'YDR418W'
    paralog2 = 'YEL054C'

    paralog = [paralog1, paralog2]
    newicktree = './YeastTree.newick'
    sim_num = 1
    num_exon = 163
    log_folder = './sim1/log/'
    div_folder = './sim1/div/'
    seed_number = 27606
    

    tau = 1.409408

    
    IGC_geo  = 3.0 / 3.0
    IGC_init = tau / IGC_geo 
    IGC_threshold = -0.1
    x_IGC = [IGC_init, 1.0 / IGC_geo, IGC_threshold]  # These values vary for the simulation study

    seq_file = './sim' + str(sim_num) + '/' + '_'.join(paralog) + '_MG94_geo_' + str(IGC_geo) + '_Sim_' + str(sim_num) + '.fasta'
    log_file = './sim' + str(sim_num) + '/' + '_'.join(paralog) + '_MG94_geo_' + str(IGC_geo) + '_Sim_' + str(sim_num) + '.log'

    save_file = './save/MG94_' + '_'.join(paralog) + '_nonclock_save.txt'

    #x_exon = [0.49249355302375575, 0.60985035555249456, 0.42155795722934408, 8.1662933909645563, 0.092804167727196338]
    pi_a = 2.983653651297249465e-01
    pi_c = 2.095729139806868369e-01
    pi_g = 1.871536409219771990e-01
    pi_t = 3.049080799676109899e-01
    kappa = 8.404333642032199236
    omega = 7.618529705823594289e-02
    x_exon = [pi_a + pi_g, pi_a / (pi_a + pi_g), pi_c / (pi_c + pi_t), kappa, omega]  # parameter from MG94_YDR418W_YEL054C_nonclock_summary.txt

    save_file = './save/MG94_' + '_'.join(paralog) + '_nonclock_save.txt'
    x_rates = np.exp([-3.925914311698581738, -1.533949335440421891, -1.564219363086546410, -1.762095498829083562, -3.660826292734839171, -3.626559649929371076,
               -3.438639957077882059, -2.460890888502840657, -3.690966112631306473, -2.870638256407853639, -2.844813362532052192, -3.822236386160084542])

    test = TreeIGCCodonSimulator(num_exon, newicktree, paralog, seq_file, log_file, x_exon, x_IGC, log_folder, div_folder, seed_number)
    test.unpack_x_rates(x_rates)
    self = test

    test.sim_root()
    self.sim()

##    for edge in self.edge_list:
##        print edge
##        if edge in self.outgroup:
##            rate_mat = self.mut_Q
##        
##        # Now need to adapt branchsim
##        # create an instance for each branch
##        blen = self.edge_to_blen[edge]
##        num_exon = self.num_exon
##        num_intron = 0
##        x_exon = self.x_exon
##        x_IGC  = self.x_IGC
##        x_intron = self.x_exon
##        log_file = self.log_folder + 'log.log'
##        div_file = self.div_folder + 'div.log'
##        starting_seq = self.node_to_sequence[edge[0]]
##        
##        branch_sim = OneBranchIGCSimulator(blen = blen, num_exon = num_exon, num_intron = num_intron,
##                                           x_exon = x_exon, x_IGC = x_IGC, x_intron = x_intron,
##                                           log_file = log_file, div_file = div_file, initial_seq = starting_seq)
##
##
##        print blen
##        blen = self.edge_to_blen[edge]
##        self = branch_sim
##
##        branch_sim.sim_one_branch(starting_seq, blen)
##
##    self = branch_sim
