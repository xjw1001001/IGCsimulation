#
# Code to parse the interesting bits of a PAML results file that contains
# statistics for several trees. Note, this approach contains some hacks that
# are probably specific to multiple-tree CODEML files, Biopython has a module
# for handling other PAML outputs, and it probably a better starting point for
# 'normal' (single tree) input files.
#
# TODO - should clean up handlng of resifues - use @property decorator and
# ._functions to simplify writing/reading.
#

##https://gist.github.com/dwinter/2851069


import sys
import re
from collections import defaultdict
from scipy.stats import chisqprob

class PAMLResult(object):
    """Represent the results from one PAMl model for one tree """

    def __init__(self, lnL, residues =None):
        self.lnL = lnL
        self.residues = residues

    def __repr__(self):
        return 'PAMLResult(lnL= {0})'.format(self.lnL)


class PAMLTree(dict):
    """Hold PAML restuls for one tree across several models"""

    def __init__(self, *args):
        dict.__init__(self, args)


    def _add_result(self, model, res):
        self[model] = res


    def _calculate_LRTs(self):
        """Run likelihood ratio test if there are enough results  """
        if all( [m in self.keys() for m in [1,2]] ):
            D = -2 * self[1].lnL + 2 * self[2].lnL
            pval = chisqprob(D,2)
            self.LRT_m1m2 = (D, pval)

        if all( [m in self.keys() for m in [7,8]] ):
            D = -2 * self[7].lnL + 2 * self[8].lnL
            pval = chisqprob(D,2)
            self.LRT_m7m8 = (D, pval)


class PAMLParser(object):
    """Parse an entire CODEML result file
    Collects data from all models and trees analysed and runs LRTs for model
    comparions as well as collating residues found to be under selection using
    Bayes Emperical Bayes.
    Usage
    result_file = PAMLParser("my_run.out")
    result_file.write("my_run", LRT=True, residues=False)
    will write files:
    'my_run_LRTs.csv', 'my_run_res_m2.csv' and 'my_run_resm8.csv'
    """

    def __init__(self, fname):
        self.fname = fname
        self.read()


    def _get_residues(self, model):
        """ """
        d = defaultdict(list)
        for tree, mod in self.combined.items():
            for line in mod[model].residues:
                d[line[0]].append(line[2])
        return d


    def read(self):
        self.combined = defaultdict(PAMLTree)
        lnL_pattern = re.compile('np:\s+\d+\):\s+(-\d+\.\d+)')

        mods = open(self.fname).read().split("Model ")[1:]
        for m in mods:
            m_name = int(m.split(':', 1)[0])
            trees = m.split('TREE ')[1:]
            for t in trees:
                t_name = int(t.split(':')[0][-1])
                ln_match = lnL_pattern.search(t)
                t_ln = float(ln_match.group(1))
                if "Bayes Empirical Bayes" in t:
                    beb = t.split("Bayes Empirical Bayes (BEB) analysis")[1]
                    table = beb.split("SE for w")[1].split('The grid')[0]
                    #index, default, posterior w _ SE  for residues
                    residues = [tuple(l.split()) for l in table.split("\n") if l]
                else:
                    residues = None
                res = PAMLResult(lnL=t_ln, residues=residues)
                self.combined[t_name]._add_result(m_name, res)

        #all done, if there are any LRT-tests to run we can do them now
        [p._calculate_LRTs() for p in self.combined.values()]

    def write(self, file_stem, LRT = True, residues= True):
        """ """
        if LRT:
            out = open(file_stem + "_LRTs.csv", "w")
            out.write("tree,comp,D, p-val\n")
            t = 0
            for tree, results in self.combined.items():
                try:
                    D, pval = results.LRT_m1m2
                    out.write("{0}, m1am2a, {1}, {2}\n".format(
                                                      tree,round(D,3), pval))
                    t += 1
                except(AttributeError):
                    pass #no results for that comparisom

                try:
                    D, pval = results.LRT_m7m8
                    out.write("{0}, m7m8, {1}, {2}\n".format(
                                                     tree,round(D,3), pval))
                    t += 1
                except(AttributeError):
                    pass #no results for that comparisom

            out.close()
            print "wrote data for {0} LRTs".format(t)

        if residues:

            try:
                m8 = self._get_residues(8)
                out = open(file_stem +"_res_m8.csv","w")
                header = "res, " + ",".join(
                        ["t" + str(i) for i in range(1, len(m8.values()[1])+1)]
                )
                out.write(header + "\n")
                counter = 0
                for res, vals in m8.items():
                    line = "{0},{1}\n".format(
                                    res, ",".join(f.strip("*") for f in vals))
                    out.write(line)
                    counter += 1
                print "wrote data for {0} residues using model 8".format(counter)

            except KeyError:
                pass

            try:
                m2 = self._get_residues(2)
                out = open(file_stem +"_res_m2.csv","w")
                header = "res, " + ",".join(
                        ["t" + str(i) for i in range(1, len(m2.values()[1])+1)]
                )
                out.write(header + "\n")
                counter = 0
                for res, vals in m2.items():
                    line = "{0},{1}\n".format(
                                    res, ",".join(f.strip("*") for f in vals))
                    out.write(line)
                    counter += 1
                print "wrote data for {0} residues using model 2a".format(counter)

            except KeyError:
                pass #not finised/not results for this model


def main():
    """Get the information form one file, given as a commandline argument """
    result = PAMLParser(sys.argv[1])
    result.write(sys.argv[1] + '_result')

if __name__ == "__main__":
    main()
