import sys, os
sys.path.insert(0, 'msprime')
import msprime

growth_rate = 0.3
n_gens = 17

des = [
        msprime.PopulationParametersChange(
            time=n_gens+1, initial_size=1000, growth_rate=0)
      ]

pcs = [
        msprime.PopulationConfiguration(
            sample_size=100000, initial_size=1.4e6, growth_rate=0.3)
      ]

dd = msprime.DemographyDebugger(population_configurations=pcs, demographic_events=des)
dd.print_history()
# sys.exit()

ts = msprime.simulate(model='dtwf', length=3e8, population_configurations=pcs,
        recombination_rate=1e-8, mutation_rate=1e-8, demographic_events=des)

ts.dump(os.path.expanduser('~/project/pedigree_msp/results/wf_100Ksamples_3M.h5'))

import IPython; IPython.embed()
