#!/usr/bin/env python3

import argparse, sys, os, re

class SmartFormatter(argparse.HelpFormatter):

    def _split_lines(self, text, width):
        if text.startswith('R|'):
            return text[2:].splitlines()
        # this is the RawTextHelpFormatter._split_lines
        return argparse.HelpFormatter._split_lines(self, text, width)

from argparse import ArgumentParser

def main(raw_args=None):
    parser = argparse.ArgumentParser(description="Create libraries file for each sample given input csv", formatter_class=SmartFormatter)
    parser.add_argument('file_name', metavar='libraries.csv', action = "store",
        type=str, help="CSV file containing the sample name, fastq path, and the library type for each file. Expected column order in file: \n(Sample) Name, Flowcell, Sample, Library Type")
    parser.add_argument('fastqs', metavar="outs/fastq_path/HLGW3DRXX",
        nargs='?', action="store", type=str,
        help="Full path to FASTQ files, used to help fill in new library files. If missing will just use values in the libraries.csv file. Multiple paths should be comma delimited")

    args = parser.parse_args(raw_args)
    if args.fastqs != None:
        fastqs = args.fastqs.split(',')

    with open(args.file_name) as f:
        headers = next(f).strip().split(',')
        #print(headers)
        try:
            name_index = [i for i in range(len(headers)) if headers[i].casefold() == 'Name'.casefold()][0] 
            flowcell_index = [i for i in range(len(headers)) if headers[i].casefold() == 'Flowcell'.casefold()][0]
            sample_index = [i for i in range(len(headers)) if headers[i].casefold() == 'Sample'.casefold()][0]
            library_index = [i for i in range(len(headers)) if headers[i].casefold() == 'Type'.casefold()][0]
        except:
            print("File headers could not be parsed. Please check that they match the expected format of: Name,Flowcell,Sample,Type\n")
        samples = dict()
        for line in f:
            line = line.strip().split(',')
            if line[name_index] in samples:
                samples[line[name_index]].append([line[flowcell_index], line[sample_index], line[library_index]])
            else:
                samples[line[name_index]] = [[line[flowcell_index], line[sample_index], line[library_index]]]
    for sample in samples:
        text = []
        for values in samples[sample]:
            if args.fastqs != None:
                runs = [path for path in fastqs if values[0].rstrip('/') in path]
                if len(runs) != 1:
#                    runs = [run for run in runs if values[1] in run]
                    runs = [run for run in runs for j in run.split(os.sep) if len(re.findall(values[1] + r'$', j)) > 0]
                if len(runs) != 1:
                    sys.exit("Problems finding unique match for %s in %s" % (values[0], args.fastqs))
                else:
                    print([runs[0], values[1], values[2]])
                    text.append(",".join([runs[0], values[1], values[2]]))
            else:
                text.append(",".join(values))

        with open('%s_libraries.csv' % sample, 'w') as f:
            f.write('fastqs,sample,library_type\n')
            f.write('\n'.join(text))

    # print(samples)


if __name__ == '__main__':
    main()
