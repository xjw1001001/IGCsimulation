import os
import subprocess
from Bio import Seq, SeqIO, AlignIO, Phylo
from cStringIO import StringIO
from Bio.Phylo.PAML import codeml, baseml
import numpy as np
import networkx as nx


def initialize(paralog, out_path = './output/', alignment_path = '../MafftAlignment/'):
    if not os.path.isdir(out_path + '_'.join(paralog)):
        os.mkdir(out_path + '_'.join(paralog))

    input_alignment = out_path + '_'.join(paralog) + '/' + '_'.join(paralog) + '_input.fasta' 
    if not os.path.isfile(input_alignment):
        subprocess.check_output(['cp', input_alignment.replace(out_path, alignment_path), input_alignment])

    input_tree = out_path + '_'.join(paralog) + '/' + '_'.join(paralog) + '_tree.newick'
    old_paml_tree_path = '/Users/xji3/Genconv_Copy/NewClusterPackRun/NewPairsAlignment/'
    if not os.path.isfile(input_tree):
        subprocess.check_output(['cp', input_tree.replace(out_path, old_paml_tree_path), input_tree])


def run_paml(wk_dir, ctl_file, codeml_dir = '/Users/xji3/Downloads/paml4.8/bin/codeml'):
    codeml_cmd = [codeml_dir, ctl_file.replace(wk_dir, './')]
    os.chdir(wk_dir)
    print(codeml_cmd)
    subprocess.check_output(codeml_cmd)
    

def prepare_ctl(tree_loc, seq_loc, out_file, ctl_loc):
    fixed_stuff = ['noisy = 9', 'verbose = 1', 'runmode = 0', 'seqtype = 1',
                   'CodonFreq = 4', 'estFreq = 1', 'ndata = 1', 'clock = 0', 'aaDist = 0',
                   'model = 0', 'NSsites = 0', 'icode = 0', 'Mgene = 0',
                   'fix_kappa = 0', 'kappa = 2', 'fix_omega = 0', 'omega = 0.4',
                   'fix_alpha = 1', 'alpha = 0.', 'Malpha = 0', 'ncatG = 1',
                   'getSE = 0', 'RateAncestor =0', 'Small_Diff = .5e-6',
                   'cleandata = 0', 'fix_blength = 0', 'method = 0']
    with open(ctl_loc, 'w+') as f:
        f.write('seqfile = ' + seq_loc + '\n')
        f.write('treefile = ' + tree_loc + '\n')
        f.write('outfile = ' + out_file + '\n')
        for stuff in fixed_stuff:
            f.write(stuff + '\n')


# Copy from my code for reading newick tree
def get_tree(tree_file, name_tree):
    tree = Phylo.read( open(tree_file, 'r'), "newick")
    tree_name = Phylo.read( open(name_tree, 'r'), "newick")
    #set node number for nonterminal nodes and specify root node
    numInternalNode = 0
    for clade in tree.get_nonterminals():
        clade.name = 'N' + str(numInternalNode)
        clade.branch_length = clade.confidence
        numInternalNode += 1

    
    for clade_iter in range(len(tree.get_terminals())):
        clade = tree.get_terminals()[clade_iter]
        clade.branch_length = clade.confidence
        clade.name = tree_name.get_terminals()[clade_iter].name
    tree_phy = tree.as_phyloxml(rooted = 'True')
    tree_nx = Phylo.to_networkx(tree_phy)


    triples = ((u.name, v.name, d['weight']) for (u, v, d) in tree_nx.edges(data = True)) # data = True to have the blen as 'weight'
    T = nx.DiGraph()
    edge_to_blen = {}
    for va, vb, blen in triples:
        edge = (va, vb)
        T.add_edge(*edge)
        edge_to_blen[edge] = blen

    edge_list = edge_to_blen.keys()
    edge_list.sort(key = lambda node: int(node[0][1:]))

    return edge_to_blen, edge_list

    

if __name__ == '__main__':

    
##    finished_list = []


    IGC_geo_list = [3.0, 10.0, 50.0, 100.0, 500.0]
    tree_loc = '/Users/xji3/GitFolders/IGCCodonSimulation/YDR418W_YEL054C_tree.newick'

    IGC_geo_list = [500.0]
    name_tree = '/Users/xji3/GitFolders/IGCCodonSimulation/YDR418W_YEL054C.newick'

    for IGC_geo in IGC_geo_list:
        label = ['ll', 'kappa', 'omega']
        header = []
        summary_mat = []
        for sim_num in range(100):
            wk_dir = '/Users/xji3/GitFolders/IGCCodonSimulation/YDR418W_YEL054C/IGCgeo_' + str(IGC_geo) + '/sim_' + str(sim_num) + '/'
            seq_loc = wk_dir + 'YDR418W_YEL054C_MG94_geo_' + str(IGC_geo) + '_Sim_' + str(sim_num) + '.fasta'
            ctl_loc = wk_dir + 'geo_' + str(IGC_geo) + '_Sim_' + str(sim_num) + '_codeml.ctl'
            out_file = wk_dir + 'unrooted_MG94_geo_' + str(IGC_geo) + '_Sim_' + str(sim_num) + '_codeml_output.txt'
            prepare_ctl(tree_loc, seq_loc, out_file, ctl_loc)
            run_paml(wk_dir, ctl_loc)

            if os.path.isfile(out_file):
                codeml_result = codeml.read(out_file)
                tree_file = out_file.replace('codeml_output.txt', 'codeml_tree_est.newick')
                with open(tree_file, 'w+') as f:
                    f.write(codeml_result['NSsites'][0]['tree'] + '\n')

                edge_to_blen, edge_list = get_tree(tree_file, name_tree)
                summary = [codeml_result['NSsites'][0]['lnL'],
                           codeml_result['NSsites'][0]['parameters']['kappa'],
                           codeml_result['NSsites'][0]['parameters']['omega']]
                summary.extend([edge_to_blen[edge] for edge in edge_list])
                summary_mat.append(summary)
                header.append('geo_' + str(IGC_geo) + '_sim_' + str(sim_num))

        
        label.extend(['_'.join(edge) for edge in edge_list])
        print len(header), len(label)
        footer = ' '.join(label)
        header = ' '.join(header)
        np.savetxt(open('./geo_' + str(IGC_geo) + '_paml_unrooted_summary.txt', 'w+'), np.matrix(summary_mat).T, delimiter = ' ', footer = footer, header = header)
            


##    #pairs = pairs[0:2]
##    for pair in pairs:
##        codeml_result = codeml.read('/Users/xji3/GitFolders/Genconv/PAMLCheck/output/' + '_'.join(pair) + '/' + '_'.join(pair) + '_codeml_result.txt')
##        summary_mat.append([codeml_result['NSsites'][0]['tree length'],
##                            codeml_result['NSsites'][0]['lnL']])
##        finished_list.append(pair)
##
##    header = ' '.join(['_'.join(pair) for pair in finished_list])  # column labels
##    np.savetxt(open('/Users/xji3/GitFolders/Genconv/PAMLCheck/paml_summary.txt', 'w+'), np.matrix(summary_mat).T, delimiter = ' ', footer = footer, header = header)
