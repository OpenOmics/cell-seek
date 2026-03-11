#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Jul 31 2018

@author: Vicky

usage: python generateSummaryFiles.py metric_summary

"""


import glob, xlsxwriter, csv, sys, ntpath, os, re
from shutil import copyfile

metricsPath = 'finalreport/'
summaryPath = 'finalreport/summaries/'
def main(arg1='metric_summary'):
    createMetricsSummary(arg1)
    copyWebSummary()

def createMetricsSummary(arg1):
    try:
        os.makedirs(metricsPath)
    except OSError:
        if not os.path.isdir(metricsPath):
            raise


    files = glob.glob('./*/outs/per_sample_outs/*/metrics_summary.csv')
    files.sort()

    stats = dict()
    headers = []
    samples = list()
    for filename in files:
        if len(set([filename.split('/')[i] for i in [-5, -2]])) == 1:
            sample = filename.split('/')[-5]
        else:
            sample = '|'.join([filename.split('/')[i] for i in [-5, -2]])
        samples.append(sample)
        f = open(filename)
        temp = csv.reader(f, delimiter = ',')
        newheaders = []
        for row in temp:
            if row[0] == 'Cells' or (row[0] == 'Library' and (row[2] == 'Physical library ID' or row[1] == 'Multiplexing Capture' or row[2] == 'Fastq ID')):
                if row[2] == 'CMO Name':
                    header = ' '.join([row[1], row[-3], row[-2]])
                    newheaders.append(header)
                elif row[2] == 'Fastq ID':
                    header = ' '.join(['Fastq', row[-2]])
                    if header not in newheaders:
                        newheaders.append(header)
                else:
                    header = ' '.join([row[0], row[1], row[-2]])
                    newheaders.append(header)

                samplestats = stats.get(header, dict())
                if row[2] == 'Fastq ID':
                    samplestats[row[-3]] = row[-1]
                else:
                    samplestats[sample] = row[-1]
                stats[header] = samplestats
        if len(newheaders) > len(headers):
            if len(set(headers)-set(newheaders)) > 0:
                newheaders += list(set(headers)-set(newheaders))
            headers = newheaders
        f.close()

    samples.sort()

    workbook = xlsxwriter.Workbook(metricsPath + arg1+'.xlsx')
    workbook_clean = xlsxwriter.Workbook(metricsPath + arg1+'_clean.xlsx')
    #if len([i for i in headers if 'Multiplexing Capture' in i]) > 0:
    #    write_sheet(workbook, stats, [i for i in stats['Sample Antibody Capture Cells'] if stats['Sample Antibody Capture Cells'][i] != '0'], headers, "Sample")
    #else:
    write_sheet(workbook, stats, samples, headers, "Cells")
    write_sheet(workbook_clean, stats, samples, headers, "Cells", clean=True)
    write_sheet(workbook, stats, samples, headers, "Library")
    write_sheet(workbook_clean, stats, samples, headers, "Library")
    if len([i for i in headers if 'Multiplexing Capture' in i]) > 0:
        write_sheet(workbook, stats, samples, headers, "Multiplexing")
        write_sheet(workbook_clean, stats, samples, headers, "Multiplexing")

    #Pull out the FASTQ file names to provide as sample names
    write_sheet(workbook, stats, sorted(list(set([x for i in stats if 'Fastq' in i for x in stats[i].keys()]))), headers, "Fastq")
    write_sheet(workbook_clean, stats, sorted(list(set([x for i in stats if 'Fastq' in i for x in stats[i].keys()]))), headers, "Fastq")

    workbook.close()
    workbook_clean.close()

def write_sheet(workbook, stats, samples, headers, filter, clean=False):
    formatNum = workbook.add_format({'num_format': '#,##0'})
    formatFrac = workbook.add_format({'num_format': '#,##0.##'})
    formatPer = workbook.add_format({'num_format': '0.00%'})
    formatHead = workbook.add_format({'bold': True, 'italic': True, 'text_wrap': True, 'align': 'center'})

    worksheet = workbook.add_worksheet(filter)

    sheet_headers = [header for header in headers if header.split(' ')[0] == filter]
    print(sheet_headers)
    [worksheet.set_column(sheet_headers.index(header)+1, sheet_headers.index(header)+1, 10) for header in sheet_headers if ('Gene Expression' in header or 'cell-associated' in header.lower())]
    [header for header in sheet_headers if 'Gene Expression' in header or 'cell-associated' in header.lower()]
    [worksheet.set_column(sheet_headers.index(header)+1, sheet_headers.index(header)+1, 12) for header in sheet_headers if ('number of reads' in header.lower() or 'Multiplexing')]
    [sheet_headers.index(header)+1 for header in sheet_headers if ('number of reads' in header.lower() or 'Multiplexing')]

    row = 1
    printed = list()
    for sample in samples:
        if clean:
            if stats['Cells Gene Expression Cells'][sample] == '0':
                continue
        if '|' in sample:
            if filter in ['Library', 'Multiplexing']:
                if sample.split('|')[0] in printed:
                    continue
                printed.append(sample.split('|')[0])
        else:
            printed.append(sample)
        if filter in ['Library', 'Multiplexing']:
            worksheet.write(row, 0, printed[-1].replace('|', ' '))
        else:
            worksheet.write(row, 0, sample.replace('|', ' '))

        col = 1
        for category in sheet_headers:
            if sample in stats[category]:
                i = stats[category][sample]
                if '%' in i:
                    if ' ' in i:
                        worksheet.write(row, col, i)
                    elif '---' in i:
                        worksheet.write(row, col, i)
                    else:
                        worksheet.write(row, col, float(i.strip('%'))/100, formatPer)
                elif '.' in i:
                    if i.count('.') > 1:
                        worksheet.write(row, col, i.replace(',', ''))
                    else:
                        if float(i.replace(',','')) <= 1:
                            worksheet.write(row, col, float(i.replace(',','')), formatPer)
                        else:
                            worksheet.write(row, col, float(i.replace(',','')), formatFrac)
                elif len(re.findall('\D', ''.join(re.findall('[^,]', i)))) > 0:
                    worksheet.write(row, col, i)
                else:
                    worksheet.write(row, col, int(i.replace(',','')), formatNum)
            col += 1
        row += 1

    col = 1
    row = 0
    worksheet.write(0, 0, "Sample", formatHead)
    for i in sheet_headers:
        worksheet.write(row, col, ' '.join(i.split(' ')[1:]), formatHead)
        col += 1

def copyWebSummary():
    try:
        os.makedirs(summaryPath)
    except OSError:
        if not os.path.isdir(summaryPath):
            raise
    files = glob.glob('./*/outs/per_sample_outs/*/web_summary.html')
    for filename in files:
        if len(set([filename.split('/')[i] for i in [1, -2]])) == 1:
            new_name = filename.split('/')[1]
        else:
            new_name = '_'.join([filename.split('/')[i] for i in [1, -2]])
        copyfile(filename, '%s/%s_web_summary.html' % (summaryPath, new_name))

if __name__ == "__main__":
    if len(sys.argv) == 1:
        main()
    else:
        main(sys.argv[1])
