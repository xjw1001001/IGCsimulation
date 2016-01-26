from Rewrite_CodonGeneconv import ReCodonGeneconv
import os
import argparse

def main(args):
    paralog1 = 'YDR418W'
    paralog2 = 'YEL054C'
    paralog = [paralog1, paralog2]
    newicktree = './YeastTree.newick'

    IGC_geo = args.IGC_geo
    sim_num = args.sim_num

    alignment_file = './' + '_'.join(paralog) + '/IGCgeo_' + str(IGC_geo) + '/sim_' + str(sim_num) + '/' + '_'.join(paralog) + '_MG94_geo_' + str(IGC_geo) + '_Sim_' + str(sim_num) + '.fasta'
    save_name = './SimulationSave/' + '_'.join(paralog) + '/IGCgeo_' + str(IGC_geo) + '/' + '_'.join(paralog) + '_MG94_geo_'  + str(IGC_geo) + '_Sim_' + str(sim_num) + '_save.txt'
    summary_name = './SimulationSummary/' + '_'.join(paralog) + '/IGCgeo_' + str(IGC_geo) + '/' + '_'.join(paralog) + '_MG94_geo_'  + str(IGC_geo) + '_Sim_' + str(sim_num) + '_summary.txt'

    # generate folder if not exist
    if not os.path.isdir('./SimulationSave/' + '_'.join(paralog) + '/'):
        os.mkdir('./SimulationSave/' + '_'.join(paralog) + '/')

    if not os.path.isdir('./SimulationSave/' + '_'.join(paralog) + '/IGCgeo_' + str(IGC_geo) + '/'):
        os.mkdir('./SimulationSave/' + '_'.join(paralog) + '/IGCgeo_' + str(IGC_geo) + '/')

    if not os.path.isdir('./SimulationSummary/' + '_'.join(paralog) + '/'):
        os.mkdir('./SimulationSummary/' + '_'.join(paralog) + '/')

    if not os.path.isdir('./SimulationSummary/' + '_'.join(paralog) + '/IGCgeo_' + str(IGC_geo) + '/'):
        os.mkdir('./SimulationSummary/' + '_'.join(paralog) + '/IGCgeo_' + str(IGC_geo) + '/')
        
    test = ReCodonGeneconv( newicktree, alignment_file, paralog, Model = 'MG94', Force = None, clock = False, save_name = save_name)

    test.get_mle(False, True, 0, 'BFGS')
    test.get_individual_summary(summary_path = './SimulationSummary/' + '_'.join(paralog) + '/', file_name = summary_name)
    
    


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--Geo', dest = 'IGC_geo', help = 'IGC tract length parameter')
    parser.add_argument('--sim_num', dest = 'sim_num', help = 'Simulation number')
    
    main(parser.parse_args())    
