import os
import subprocess
from Bio import Seq, SeqIO, AlignIO, Phylo
from cStringIO import StringIO
from Bio.Phylo.PAML import codeml, baseml
import numpy as np
import networkx as nx
from PamlCheck import run_paml, prepare_ctl, get_tree, Seperate_codeml_result

def ReverseFastaOrder(input_fasta_file, output_fasta_flie):
    if os.path.isfile(input_fasta_file):
        with open(input_fasta_file, 'r') as f:
            all_lines = f.readlines()
            
        with open(output_fasta_file, 'w+') as g:
            for num_line in list(reversed(range(len(all_lines)))):
                if all_lines[num_line][0] == '>':
                    g.write(all_lines[num_line])
                    g.write(all_lines[num_line + 1])


if __name__ == '__main__':
##    input_fasta_file = '/Users/Xiang/GitFolders/IGCCodonSimulation/YDR418W_YEL054C_estimatedTau/IGCgeo_3.0/sim_0/YDR418W_YEL054C_MG94_geo_3.0_Sim_0.fasta'
##    output_fasta_file = '/Users/Xiang/GitFolders/IGCCodonSimulation/YDR418W_YEL054C_estimatedTau/IGCgeo_3.0/sim_0/YDR418W_YEL054C_MG94_geo_3.0_Sim_0_reversed.fasta'
##    ReverseFastaOrder(input_fasta_file, output_fasta_file)
    #tree_loc = '/Users/xji3/GitFolders/IGCCodonSimulation/YDR418W_YEL054C_tree.newick'
    tree_list = ['/Users/Xiang/GitFolders/IGCCodonSimulation/YDR418W_YEL054C_local_tree.newick',
                 '/Users/Xiang/GitFolders/IGCCodonSimulation/YDR418W_YEL054C_local_tree2.newick',
                 '/Users/Xiang/GitFolders/IGCCodonSimulation/YDR418W_YEL054C_local_tree3.newick']

    nametree_list = ['/Users/Xiang/GitFolders/IGCCodonSimulation/YDR418W_YEL054C_local_nametree.newick',
                 '/Users/Xiang/GitFolders/IGCCodonSimulation/YDR418W_YEL054C_local_nametree2.newick',
                 '/Users/Xiang/GitFolders/IGCCodonSimulation/YDR418W_YEL054C_local_nametree3.newick']
    IGC_geo_list = [3.0]

    local_tree_num = 1
    tree_loc = tree_list[local_tree_num - 1]
    name_tree = nametree_list[local_tree_num - 1]
            
    for IGC_geo in IGC_geo_list:
        label = ['ll', 'kappa', 'omega']
        header = []
        summary_mat = []
        for sim_num in range(100):
            #wk_dir = '/Users/xji3/GitFolders/IGCCodonSimulation/YDR418W_YEL054C/IGCgeo_' + str(IGC_geo) + '/sim_' + str(sim_num) + '/'
            wk_dir = '/Users/Xiang/GitFolders/IGCCodonSimulation/YDR418W_YEL054C_estimatedTau/IGCgeo_' + str(IGC_geo) + '/sim_' + str(sim_num) + '/'
            seq_loc = wk_dir + 'YDR418W_YEL054C_MG94_geo_' + str(IGC_geo) + '_Sim_' + str(sim_num) + '.fasta'
            ctl_loc = wk_dir + 'geo_' + str(IGC_geo) + '_Sim_' + str(sim_num) + '_localTree_' + str(local_tree_num) + '_codeml.ctl'
            out_file = wk_dir + 'unrooted_MG94_geo_' + str(IGC_geo) + '_Sim_' + str(sim_num) + '_localTree_' + str(local_tree_num) + '_codeml_output.txt'
            prepare_ctl(tree_loc, seq_loc, out_file, ctl_loc)
            run_paml(wk_dir, ctl_loc, "/Users/Xiang/Software/paml4.8/bin/codeml")
            
            if os.path.isfile(out_file):
                codeml_result = codeml.read(out_file)
                tree_file = out_file.replace('codeml_output.txt', 'codeml_est.newick')
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
        np.savetxt(open('./geo_' + str(IGC_geo) + '_estimatedTau_paml_unrooted_localTree_' + str(local_tree_num) +'_summary.txt', 'w+'), np.matrix(summary_mat).T, delimiter = ' ', footer = footer, header = header)
     
                
