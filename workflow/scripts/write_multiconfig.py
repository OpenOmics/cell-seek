#!/usr/bin/env /mnt/nasapps/development/python/3.7.1/bin/python

import argparse, csv

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
#    parser.add_argument("--cell", metavar = 3000, default=3000,
#        nargs='?', action = "store", type=int,
#        help="Cell count")
    parser.add_argument("--cmoref", metavar="cmo_ref.csv",
        nargs='?', action = "store", type=str,
        help="Path to cmo reference file if applicable. If not included but cmosample provided will use default 10x CMOs")
    parser.add_argument("--cmosample", metavar="cmo_sample.csv",
        nargs='?', action = "store", type=str,
        help="Path to cmo file defining sample per hashtag if applicable. If not included but cmoref provided will use CMO names from cmoref")
    parser.add_argument("--feature", metavar="feature.csv",
        nargs='?', action = "store", type=str,
        help="Path to feature barcode reference file if applicable")
    parser.add_argument("--vdjref", metavar="refdata-cellranger-vdj-GRCh38-alts-ensembl-5.0.0",
        nargs='?', action = "store", type=str,
        help="Path to vdj reference")
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
#        if args.force:
#            spamwriter.writerow(['force-cells', args.cell])
#        else:
#            spamwriter.writerow(['expect-cells', args.cell])
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

        if args.cmoref != None or args.cmosample != None:
            spamwriter.writerow([])
            spamwriter.writerow(['[samples]'])
            spamwriter.writerow(['sample_id', 'cmo_ids', 'description'])
            if args.cmosample == None:
                with open(args.cmoref, 'r') as lib:
                    line = next(lib)
                    index = 1
                    for line in lib:
                        line = line.strip().split(',')
                        spamwriter.writerow(['HTO_%s' % index, line[0]])
                        index += 1
            else:
                with (open(args.cmosample, 'r')) as lib:
                    line = next(lib)
                    for line in lib:
                        spamwriter.writerow(line.strip().split(','))




if __name__ == '__main__':
    main()
