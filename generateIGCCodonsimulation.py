from IGCCodonSimTree import *
import os


# This file is used to generate simulation data
# Xiang Ji
# xji3@ncsu.edu

# again, I am too lazy and I know it...

if __name__ == '__main__':
    # constants in the simulation
    # mostly from  inference with IGC expansion see MG94_YDR418W_YEL054C_nonclock_summary.txt
    paralog1 = 'YDR418W'
    paralog2 = 'YEL054C'
    paralog = [paralog1, paralog2]
    newicktree = './YeastTree.newick'
    num_exon = 163
    tau = 1.409408
    IGC_threshold = -0.1

    IGC_geo_list = [3.0, 10.0, 50.0, 100.0, 500.0]
    IGC_geo_list = [3.0]
    for IGC_geo in IGC_geo_list:
        IGC_geo_codon = IGC_geo / 3.0
        IGC_init = tau / IGC_geo_codon
        x_IGC = [IGC_init, 1.0 / IGC_geo_codon, IGC_threshold]  # These values vary for the simulation study
        
        #sim_num = 1
        for sim_num in range(99, 100):
            log_folder = './' + '_'.join(paralog) + '/IGCgeo_' + str(IGC_geo) + '/sim_' + str(sim_num) + '/log/'
            div_folder = './' + '_'.join(paralog) + '/IGCgeo_' + str(IGC_geo) + '/sim_' + str(sim_num) + '/div/'

            # Now create folder if they don't exist
            
            if not os.path.isdir('./' + '_'.join(paralog) + '/IGCgeo_' + str(IGC_geo) + '/'):
                os.mkdir('./' + '_'.join(paralog) + '/IGCgeo_' + str(IGC_geo) + '/')

            if not os.path.isdir('./' + '_'.join(paralog) + '/IGCgeo_' + str(IGC_geo) + '/sim_' + str(sim_num) + '/'):
                os.mkdir('./' + '_'.join(paralog) + '/IGCgeo_' + str(IGC_geo) + '/sim_' + str(sim_num) + '/')

            if not os.path.isdir(log_folder):
                os.mkdir(log_folder)

            if not os.path.isdir(div_folder):
                os.mkdir(div_folder)

            
            seq_file = './' + '_'.join(paralog) + '/IGCgeo_' + str(IGC_geo) + '/sim_' + str(sim_num) + '/' + '_'.join(paralog) + '_MG94_geo_' + str(IGC_geo) + '_Sim_' + str(sim_num) + '.fasta'
            log_file = './' + '_'.join(paralog) + '/IGCgeo_' + str(IGC_geo) + '/sim_' + str(sim_num) + '/' + '_'.join(paralog) + '_MG94_geo_' + str(IGC_geo) + '_Sim_' + str(sim_num) + '.log'

            save_file = './save/MG94_' + '_'.join(paralog) + '_nonclock_save.txt'

            #x_exon = [0.49249355302375575, 0.60985035555249456, 0.42155795722934408, 8.1662933909645563, 0.092804167727196338]
            pi_a = 2.983653651297249465e-01
            pi_c = 2.095729139806868369e-01
            pi_g = 1.871536409219771990e-01
            pi_t = 3.049080799676109899e-01
            kappa = 8.404333642032199236
            #omega = 7.618529705823594289e-02
            omega = 1.0
            x_exon = [pi_a + pi_g, pi_a / (pi_a + pi_g), pi_c / (pi_c + pi_t), kappa, omega]  # parameter from MG94_YDR418W_YEL054C_nonclock_summary.txt

            save_file = './save/MG94_' + '_'.join(paralog) + '_nonclock_save.txt'
            x_rates = np.exp([-3.925914311698581738, -1.533949335440421891, -1.564219363086546410, -1.762095498829083562, -3.660826292734839171, -3.626559649929371076,
                       -3.438639957077882059, -2.460890888502840657, -3.690966112631306473, -2.870638256407853639, -2.844813362532052192, -3.822236386160084542])

            test = TreeIGCCodonSimulator(num_exon, newicktree, paralog, seq_file, log_file, x_exon, x_IGC, log_folder, div_folder)
            test.unpack_x_rates(x_rates)
##            try:
##                test.sim()
##            except:
##                print 'failed at sim  ' + str(sim_num) + '  IGC_geo = ' + str(IGC_geo)
##                test.write_log()
##        

