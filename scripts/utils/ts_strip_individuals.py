import os
import sys
sys.path.insert(0, '../msprime')
import msprime

if len(sys.argv) != 3:
    print("Usage: python ts_strip_individuals.py ts_file out_file")
    sys.exit()

ts_file = os.path.expanduser(sys.argv[1])
out_file = os.path.expanduser(sys.argv[2])

ts = msprime.load(ts_file)

if ts.num_individuals > 0:
    ll_tables = ts.dump_tables().asdict()
    ll_tables['individuals'] = msprime.tskit.IndividualTable().asdict()
    ll_tables['nodes']['individual'].fill(-1)
    ts_noinds = msprime.tskit.tables.TableCollection.fromdict(
            ll_tables).tree_sequence()
    ts_noinds.dump(out_file)
else:
    ts.dump(out_file)
