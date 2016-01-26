import os

if __name__ == '__main__':
    paralog1 = 'YDR418W'
    paralog2 = 'YEL054C'
    paralog = [paralog1, paralog2]

    IGC_geo_list = [1.0, 10.0, 50.0, 100.0, 500.0]
    sim_num_list = range(100)

    sh_line = 'sbatch -o IGCSim-%j.out --mail-type=FAIL --mail-user=xji3@ncsu.edu ./ShFiles/'

    for IGC_geo in IGC_geo_list:
        IGC_geo_sh_file = './' + '_'.join(paralog) + '_IGCgeo_' + str(IGC_geo) + '.sh'
        with open(IGC_geo_sh_file, 'w+') as f:
            f.write('#!/bin/bash' + '\n')
            for num_sim in sim_num_list:
                f.write(sh_line + '_'.join(paralog) + '_IGCgeo_' + str(IGC_geo) + '_sim_' + str(num_sim) + '.sh \n')
                with open('./ShFiles/' + '_'.join(paralog) + '_IGCgeo_' + str(IGC_geo) + '_sim_' + str(num_sim) + '.sh', 'w+') as g:
                    g.write('#!/bin/bash' + '\n')
                    g.write('python RunSimulation.py --Geo ' + str(IGC_geo) + ' --sim_num ' + str(num_sim) + '\n')


##    sh_line = 'sbatch -o IGCSim-%j.out --mail-type=FAIL --mail-user=xji3@ncsu.edu ./ShFiles_multiply3/'
##
##    for IGC_geo in IGC_geo_list:
##        IGC_geo_sh_file = './' + '_'.join(paralog) + '_IGCgeo_' + str(IGC_geo) + '_multiply3.sh'
##        with open(IGC_geo_sh_file, 'w+') as f:
##            f.write('#!/bin/bash' + '\n')
##            for num_sim in sim_num_list:
##                f.write(sh_line + '_'.join(paralog) + '_IGCgeo_' + str(IGC_geo) + '_sim_' + str(num_sim) + '_multiply3.sh \n')
##                with open('./ShFiles_multiply3/' + '_'.join(paralog) + '_IGCgeo_' + str(IGC_geo) + '_sim_' + str(num_sim) + '_multiply3.sh', 'w+') as g:
##                    g.write('#!/bin/bash' + '\n')
##                    g.write('python RunSimulation_multiply3.py --Geo ' + str(IGC_geo) + ' --sim_num ' + str(num_sim) + '\n')
