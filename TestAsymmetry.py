import os, subprocess

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
    input_fasta_file = '/Users/Xiang/GitFolders/IGCCodonSimulation/YDR418W_YEL054C_estimatedTau/IGCgeo_3.0/sim_0/YDR418W_YEL054C_MG94_geo_3.0_Sim_0.fasta'
    output_fasta_file = '/Users/Xiang/GitFolders/IGCCodonSimulation/YDR418W_YEL054C_estimatedTau/IGCgeo_3.0/sim_0/YDR418W_YEL054C_MG94_geo_3.0_Sim_0_reversed.fasta'
    ReverseFastaOrder(input_fasta_file, output_fasta_file)
    
                
