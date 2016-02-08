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

def Seperate_codeml_result(codeml_output_file, new_files):
    if os.path.isfile(codeml_output_file):
        with open(codeml_output_file, 'r') as f:
            all_lines = f.readlines()
        tree_start = [i for i in range(len(all_lines)) if all_lines[i][:6] == 'TREE #']
    assert(len(new_files) == len(tree_start))
    basic_lines = all_lines[:tree_start[0]]
    for fiter in range(len(new_files)):
        new_file = new_files[fiter]
        if fiter == len(new_files) - 1:
            tree_result = all_lines[tree_start[fiter]:]
        else:
            tree_result = all_lines[tree_start[fiter]:tree_start[fiter + 1]]
        with open(new_file, 'w+') as f:
            f.writelines(basic_lines)
            f.writelines(tree_result)
            

if __name__ == '__main__':

    
##    finished_list = []


##    IGC_geo_list = [3.0, 10.0, 50.0, 100.0, 500.0]
##    tree_loc = '/Users/xji3/GitFolders/IGCCodonSimulation/YDR418W_YEL054C_tree.newick'
##
##    IGC_geo_list = [500.0]
##    name_tree = '/Users/xji3/GitFolders/IGCCodonSimulation/YDR418W_YEL054C.newick'

    IGC_geo_list = [3.0, 10.0, 50.0, 100.0, 500.0]
    tree_loc = '/Users/xji3/GitFolders/IGCCodonSimulation/YDR418W_YEL054C_tree.newick'

    IGC_geo_list = [3.0]
    name_tree_1st = '/Users/xji3/GitFolders/IGCCodonSimulation/YDR418W_YEL054C_1st.newick'
    name_tree_2nd = '/Users/xji3/GitFolders/IGCCodonSimulation/YDR418W_YEL054C_2nd.newick'

    for IGC_geo in IGC_geo_list:
        label = ['ll', 'kappa', 'omega']
        header = []
        summary_mat = []
        summary_mat_2 = []
        for sim_num in range(100):
            #wk_dir = '/Users/xji3/GitFolders/IGCCodonSimulation/YDR418W_YEL054C/IGCgeo_' + str(IGC_geo) + '/sim_' + str(sim_num) + '/'
            wk_dir = '/Users/xji3/GitFolders/IGCCodonSimulation/YDR418W_YEL054C_10Tau/IGCgeo_' + str(IGC_geo) + '/sim_' + str(sim_num) + '/'
            seq_loc = wk_dir + 'YDR418W_YEL054C_MG94_geo_' + str(IGC_geo) + '_Sim_' + str(sim_num) + '.fasta'
            ctl_loc = wk_dir + 'geo_' + str(IGC_geo) + '_Sim_' + str(sim_num) + '_codeml.ctl'
            out_file = wk_dir + 'unrooted_MG94_geo_' + str(IGC_geo) + '_Sim_' + str(sim_num) + '_codeml_output.txt'
##            prepare_ctl(tree_loc, seq_loc, out_file, ctl_loc)
##            run_paml(wk_dir, ctl_loc, "/Users/xji3/Software/paml4.8/bin/codeml")
            out_tree1_file = out_file.replace('_output.txt', '_tree1_output.txt')
            out_tree2_file = out_file.replace('_output.txt', '_tree2_output.txt')
            out_tree_files = [out_tree1_file, out_tree2_file]
            #Seperate_codeml_result(out_file, out_tree_files)

            if os.path.isfile(out_tree1_file):
                codeml_result = codeml.read(out_tree1_file)
                tree1_file = out_file.replace('codeml_output.txt', 'codeml_tree1_est.newick')
                with open(tree1_file, 'w+') as f:
                    f.write(codeml_result['NSsites'][0]['tree'] + '\n')

                edge_to_blen, edge_list_1 = get_tree(tree1_file, name_tree_1st)
                summary = [codeml_result['NSsites'][0]['lnL'],
                           codeml_result['NSsites'][0]['parameters']['kappa'],
                           codeml_result['NSsites'][0]['parameters']['omega']]
                summary.extend([edge_to_blen[edge] for edge in edge_list_1])
                summary_mat.append(summary)
                header.append('geo_' + str(IGC_geo) + '_sim_' + str(sim_num))

            edge_to_blen = None                
            if os.path.isfile(out_tree2_file):
                codeml_result = codeml.read(out_tree2_file)
                tree2_file = out_file.replace('codeml_output.txt', 'codeml_tree2_est.newick')

                with open(out_tree2_file, 'r') as f:
                    all_lines = f.readlines()
                tree_line = [i for i in all_lines if i[:17] == '(kluyveriYDR418W:'][0]

                with open(tree2_file, 'w+') as f:
                    f.write(tree_line)

                edge_to_blen, edge_list_2 = get_tree(tree2_file, name_tree_2nd)
                summary = [codeml_result['NSsites'][0]['lnL'],
                           codeml_result['NSsites'][0]['parameters']['kappa'],
                           codeml_result['NSsites'][0]['parameters']['omega']]
                summary.extend([edge_to_blen[edge] for edge in edge_list_2])
                summary_mat_2.append(summary)
                #header.append('geo_' + str(IGC_geo) + '_sim_' + str(sim_num))
        
        label.extend(['_'.join(edge) for edge in edge_list_1])
        print len(header), len(label)
        footer = ' '.join(label)
        header = ' '.join(header)
        np.savetxt(open('./geo_' + str(IGC_geo) + '_10Tau_paml_unrooted_1stTree_summary.txt', 'w+'), np.matrix(summary_mat).T, delimiter = ' ', footer = footer, header = header)
        label = ['ll', 'kappa', 'omega']
        label.extend(['_'.join(edge) for edge in edge_list_2])
        footer = ' '.join(label)
        np.savetxt(open('./geo_' + str(IGC_geo) + '_10Tau_paml_unrooted_2ndTree_summary.txt', 'w+'), np.matrix(summary_mat_2).T, delimiter = ' ', footer = footer, header = header)


            
        label = ['ll', 'kappa', 'omega']
        header = []
        summary_mat = []
        summary_mat_2 = []
        for sim_num in range(100):
            #wk_dir = '/Users/xji3/GitFolders/IGCCodonSimulation/YDR418W_YEL054C/IGCgeo_' + str(IGC_geo) + '/sim_' + str(sim_num) + '/'
            wk_dir = '/Users/xji3/GitFolders/IGCCodonSimulation/YDR418W_YEL054C_estimatedTau/IGCgeo_' + str(IGC_geo) + '/sim_' + str(sim_num) + '/'
            seq_loc = wk_dir + 'YDR418W_YEL054C_MG94_geo_' + str(IGC_geo) + '_Sim_' + str(sim_num) + '.fasta'
            ctl_loc = wk_dir + 'geo_' + str(IGC_geo) + '_Sim_' + str(sim_num) + '_codeml.ctl'
            out_file = wk_dir + 'unrooted_MG94_geo_' + str(IGC_geo) + '_Sim_' + str(sim_num) + '_codeml_output.txt'
##            prepare_ctl(tree_loc, seq_loc, out_file, ctl_loc)
##            run_paml(wk_dir, ctl_loc, "/Users/xji3/Software/paml4.8/bin/codeml")
            out_tree1_file = out_file.replace('_output.txt', '_tree1_output.txt')
            out_tree2_file = out_file.replace('_output.txt', '_tree2_output.txt')
            out_tree_files = [out_tree1_file, out_tree2_file]
            Seperate_codeml_result(out_file, out_tree_files)

            if os.path.isfile(out_tree1_file):
                codeml_result = codeml.read(out_tree1_file)
                tree1_file = out_file.replace('codeml_output.txt', 'codeml_tree1_est.newick')
                with open(tree1_file, 'w+') as f:
                    f.write(codeml_result['NSsites'][0]['tree'] + '\n')

                edge_to_blen, edge_list_1 = get_tree(tree1_file, name_tree_1st)
                summary = [codeml_result['NSsites'][0]['lnL'],
                           codeml_result['NSsites'][0]['parameters']['kappa'],
                           codeml_result['NSsites'][0]['parameters']['omega']]
                summary.extend([edge_to_blen[edge] for edge in edge_list_1])
                summary_mat.append(summary)
                header.append('geo_' + str(IGC_geo) + '_sim_' + str(sim_num))

            edge_to_blen = None                
            if os.path.isfile(out_tree2_file):
                codeml_result = codeml.read(out_tree2_file)
                tree2_file = out_file.replace('codeml_output.txt', 'codeml_tree2_est.newick')

                with open(out_tree2_file, 'r') as f:
                    all_lines = f.readlines()
                tree_line = [i for i in all_lines if i[:17] == '(kluyveriYDR418W:'][0]

                with open(tree2_file, 'w+') as f:
                    f.write(tree_line)

                edge_to_blen, edge_list_2 = get_tree(tree2_file, name_tree_2nd)
                summary = [codeml_result['NSsites'][0]['lnL'],
                           codeml_result['NSsites'][0]['parameters']['kappa'],
                           codeml_result['NSsites'][0]['parameters']['omega']]
                summary.extend([edge_to_blen[edge] for edge in edge_list_2])
                summary_mat_2.append(summary)
                #header.append('geo_' + str(IGC_geo) + '_sim_' + str(sim_num))
        
        label.extend(['_'.join(edge) for edge in edge_list_1])
        print len(header), len(label)
        footer = ' '.join(label)
        header = ' '.join(header)
        np.savetxt(open('./geo_' + str(IGC_geo) + '_estimatedTau_paml_unrooted_1stTree_summary.txt', 'w+'), np.matrix(summary_mat).T, delimiter = ' ', footer = footer, header = header)
        label = ['ll', 'kappa', 'omega']
        label.extend(['_'.join(edge) for edge in edge_list_2])
        footer = ' '.join(label)
        np.savetxt(open('./geo_' + str(IGC_geo) + '_estimatedTau_paml_unrooted_2ndTree_summary.txt', 'w+'), np.matrix(summary_mat_2).T, delimiter = ' ', footer = footer, header = header)
  

##    #pairs = pairs[0:2]
##    for pair in pairs:
##        codeml_result = codeml.read('/Users/xji3/GitFolders/Genconv/PAMLCheck/output/' + '_'.join(pair) + '/' + '_'.join(pair) + '_codeml_result.txt')
##        summary_mat.append([codeml_result['NSsites'][0]['tree length'],
##                            codeml_result['NSsites'][0]['lnL']])
##        finished_list.append(pair)
##
##    header = ' '.join(['_'.join(pair) for pair in finished_list])  # column labels
##    np.savetxt(open('/Users/xji3/GitFolders/Genconv/PAMLCheck/paml_summary.txt', 'w+'), np.matrix(summary_mat).T, delimiter = ' ', footer = footer, header = header)
