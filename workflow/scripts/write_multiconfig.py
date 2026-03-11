#!/usr/bin/env /mnt/nasapps/development/python/3.7.1/bin/python

import argparse, csv, re

def cellLimits(s):
    try:
        #print(s)
        #print(re.split(seps, s))
        situp = []
        return(tuple(s.split(',')))
        #    situp.append(si.split(','))
        return situp
    except:
        raise argparse.ArgumentTypeError("Coordinates must be given divided by commas and space, dot, or semicolon e.g.: 'x,y k,l,m'")

def main(raw_args=None):
    parser = argparse.ArgumentParser(description="""Help to set up and run the single cell multi pipeline""")
    parser.add_argument("-o", "--output", metavar="output.csv",
        action = "store", type=str, required=True,
        help="Output file name")
    parser.add_argument("-r", "--ref", metavar="refdata-gex-GRCh38-2020-A",
        action = "store", type=str, required=True,
        help="Path to reference")
    parser.add_argument("-l", "--lib", metavar="libraries.csv",
        action = "store", type=str, required=True,
        help="Path to libraries file to create the config file based off of")
    parser.add_argument("--cellranger", metavar="8.0.0",
        action = "store", type=str, required=True,
        help="Version of CellRanger that is being run to handle changes in flags")
    parser.add_argument("--forcecells", metavar = 3000, default=None,
        nargs='?', action = "store", type=int,
        help="Cell count")
    parser.add_argument("--multiplexforcecells", metavar = "sample1,3000 sample2,5000", default=None,
        nargs='*', action = "store", type=cellLimits,
        help="Cell count")
    parser.add_argument("--cmoref", metavar="cmo_ref.csv",
        nargs='?', action = "store", type=str,
        help="Path to cmo reference file if applicable. If not included but cmosample provided will use default 10x CMOs")
    parser.add_argument("--cmosample", metavar="cmo_sample.csv",
        nargs='?', action = "store", type=str,
        help="Path to cmo file defining sample per hashtag if applicable. If not included but cmoref provided will use CMO names from cmoref")
    parser.add_argument("--htosample", metavar="hto_sample.csv",
        nargs='?', action = "store", type=str,
        help="Path to hto file defining sample per hashtag if applicable. HTO sequence information should be included as part of feature.csv")
    parser.add_argument("--ocmsample",  metavar="ocm_sample.csv", 
        nargs='?', action = "store", type=str,
        help="Path to ocm file defining sample per on-chip hashtag if applicable.")
    parser.add_argument("--probesample",  metavar="probe_sample.csv",
        nargs='?', action = "store", type=str,
        help="Path to probe barcode file defining sample per probe barcode if applicable.")
    parser.add_argument("--feature", metavar="feature.csv",
        nargs='?', action = "store", type=str,
        help="Path to feature barcode reference file if applicable")
    parser.add_argument("--vdjref", metavar="refdata-cellranger-vdj-GRCh38-alts-ensembl-5.0.0",
        nargs='?', action = "store", type=str,
        help="Path to vdj reference")
    parser.add_argument("--probeset", metavar="Chromium_Human_Transcriptome_Probe_Set_v1.0.1_GRCh38-2020-A.csv",
        nargs='?', action = "store", type=str,
        help="Path to probe barcode reference")
    parser.add_argument("--innerprimer", metavar="enrichment_primers.txt",
        nargs='?', action = "store", type=str,
        help="Text file containing one primer per line if non-10x inner enrichment primers were used")
    parser.add_argument("-f", "--force", action="store_true",
        help="Use force-cells flag instead of expect-cells")
    parser.add_argument("-i", "--exclude_introns", action="store_true",
        help="Set include-introns flag to false")
    parser.add_argument("-b", "--create_bam", action="store_true",
        help="Set flag to for creating bam file")


    args = parser.parse_args(raw_args)
    print(args)

    with open(args.output, 'w', newline="") as csvfile:
        spamwriter = csv.writer(csvfile, delimiter=',')
        spamwriter.writerow(['[gene-expression]'])
        spamwriter.writerow(['reference', args.ref])
        if args.forcecells != None:
            spamwriter.writerow(['force-cells', args.forcecells])
#        else:
#            spamwriter.writerow(['expect-cells', args.cell])
        if args.probeset != None:
            spamwriter.writerow(['probe-set', args.probeset])
        if args.cmoref != None:
            spamwriter.writerow(['cmo-set', args.cmoref])
        if args.exclude_introns:
            spamwriter.writerow(['include-introns', 'false'])
        if args.create_bam:
            if int(args.cellranger.split('.')[0]) >= 8:
                spamwriter.writerow(['create-bam', 'true'])
            else:
                spamwriter.writerow(['no-bam', 'false'])
        else:
            if int(args.cellranger.split('.')[0]) >= 8:
                spamwriter.writerow(['create-bam', 'false'])
            else:
                spamwriter.writerow(['no-bam', 'true'])

        if args.feature != None:
            spamwriter.writerow([])
            spamwriter.writerow(['[feature]'])
            spamwriter.writerow(['reference', args.feature])

        if args.vdjref != None:
            spamwriter.writerow([])
            spamwriter.writerow(['[vdj]'])
            spamwriter.writerow(['reference', args.vdjref])
            if args.innerprimer != None:
                spamwriter.writerow(['inner-enrichment-primers', args.innerprimer])

        spamwriter.writerow([])
        spamwriter.writerow(['[libraries]'])
        spamwriter.writerow(['fastq_id', 'fastqs', 'lanes', 'feature_types'])
        with open(args.lib, 'r') as lib:
            line = next(lib)
            for line in lib:
                line = line.strip().split(',')
                spamwriter.writerow([line[1], line[0], 'Any', line[2]])

        if args.cmoref != None or args.cmosample != None or args.htosample != None or args.ocmsample != None or args.probesample != None:
            spamwriter.writerow([])
            spamwriter.writerow(['[samples]'])
            if args.htosample != None:
                headers = ['sample_id', 'hashtag_ids', 'description']
            elif args.ocmsample != None:
                headers = ['sample_id', 'ocm_barcode_ids', 'description']
            elif args.probesample != None:
                headers = ['sample_id', 'probe_barcode_ids', 'description']
            else:
                headers = ['sample_id', 'cmo_ids', 'description']
            if args.multiplexforcecells != None:
                headers += ['force_cells']
            spamwriter.writerow(headers)
            if args.cmosample == None and args.cmoref != None:
                with open(args.cmoref, 'r') as lib:
                    line = next(lib)
                    index = 1
                    for line in lib:
                        line = line.strip().split(',')
                        row = ['HTO_%s' % index, line[0], '']
                        if args.multiplexforcecells != None:
                            if any([row[0] in i for i in args.multiplexforcecells]) == 1:
                                row.append([i[1] for i in args.multiplexforcecells if row[0] in i][0])
                            else:
                                row.append('')
                        spamwriter.writerow(row)
                        index += 1
            else:
                notempty = [i for i in [args.cmosample, args.htosample, args.ocmsample, args.probesample] if i != None]
                if len(notempty) == 1:
                    with open(notempty[0], 'r') as lib:
                        line = next(lib)
                        for line in lib:
                            if len(line.strip().split(',')) == 3:
                                row = line.strip().split(',')
                            else:
                                row = line.strip().split(',') + ['']
                            if args.multiplexforcecells != None:
                                if any([row[0] in i for i in args.multiplexforcecells]) == 1:
                                    row.append([i[1] for i in args.multiplexforcecells if row[0] in i][0])
                            spamwriter.writerow(row)
                else:
                    raise RuntimeError('More than one demultiplexing sample flag used when calling function. Only one can be processed in one run')



if __name__ == '__main__':
    main()
