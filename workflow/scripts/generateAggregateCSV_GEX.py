#!/usr/bin/env python3

import glob, sys, os

def main():
    files = glob.glob('*/outs/molecule_info.h5')
    aggrFile = open('AggregatedDatasets.csv', 'w')

    aggrFile.write('sample_id,molecule_h5\n')
    files.sort()
    for i in files:
        aggrFile.write('%s,%s\n' % (i.split('/')[0], os.path.abspath(i)))
    aggrFile.close()

if __name__ == "__main__":
    main()
