from sys import argv
import os
import csv

script, listdir, chdir, length, probability_threshold = argv


def ahmm2bed(file_posterior, probability_threshold):
    inp_posterior = open(listdir + file_posterior, 'r')
    out_name = file_posterior.split('/')[-1].split('.')[0] + '_ahmm2bed'
    out = open(chdir + out_name, 'w')

    probability_threshold = float(probability_threshold)
    ancestry_list = ['2,0', '1,1', '0,2']
    chrom_list = []
    ancestry_last = 0
    ancestry_state_last = ''
    chrom_last = ''
    chromStart = 1
    chromEnd = 1

    header = inp_posterior.readline()
    out.write('chrom' + '\t' + 'chromStart' + '\t' + 'chromEnd' + '\t' + 'ancestry' + '\n')

    for line in inp_posterior.readlines():
        chrom = line.split('\t')[0].strip()
        if chrom not in dic_chrom_end:
            if chromEnd != 1:
                if ancestry_last < probability_threshold:
                    out.write(str(dic_chrom_end[chrom_last]) + '\t' + 'N' + '\n')
                else:
                    out.write(str(chrom_last) + '\t' + str(chromStart) + '\t' + str(
                        dic_chrom_end[chrom_last]) + '\t' + ancestry_state_last + '\n')
            chromEnd = 1
            continue

        position = int(line.split('\t')[1].strip())
        line_ancestry_list = [float(i) for i in line.split('\t')[2:] if i.strip()]
        ancestry = max(line_ancestry_list)
        ancestry_state = ancestry_list[line_ancestry_list.index(ancestry)]

        if chrom not in chrom_list:
            if chromEnd != 1:
                if ancestry_last < probability_threshold:
                    out.write(str(dic_chrom_end[chrom_last]) + '\t' + 'N' + '\n')
                else:
                    out.write(str(chrom_last) + '\t' + str(chromStart) + '\t' + str(dic_chrom_end[chrom_last]) + '\t' + ancestry_state_last + '\n')
            chrom_list.append(chrom)
            chromStart = 1
            chromEnd = position

            if ancestry < probability_threshold:
                out.write(str(chrom) + '\t' + str(chromStart) + '\t')

        else:
            if ancestry_last < probability_threshold and ancestry >= probability_threshold:
                out.write(str(position) + '\t' + 'N' + '\n')
                chromStart = position
                chromEnd = position

            elif ancestry_last >= probability_threshold and ancestry < probability_threshold:
                    out.write(
                        str(chrom) + '\t' + str(chromStart) + '\t' + str(chromEnd) + '\t' + ancestry_state_last + '\n')
                    out.write(str(chrom) + '\t' + str(chromEnd) + '\t')
                    chromStart = chromEnd
                    chromEnd = position

            elif ancestry_last >= probability_threshold and ancestry >= probability_threshold:
                if ancestry_state == ancestry_state_last and position - chromEnd <= 4000:
                    pass
                    chromEnd = position
                elif ancestry_state != ancestry_state_last and position - chromEnd <= 4000:
                    middle_position = chromEnd + round((position - chromEnd) / 2)
                    out.write(
                        str(chrom) + '\t' + str(chromStart) + '\t' + str(middle_position) + '\t' + ancestry_state_last + '\n')
                    chromStart = middle_position
                    chromEnd = middle_position

                elif position - chromEnd > 4000:
                    out.write(
                        str(chrom) + '\t' + str(chromStart) + '\t' + str(chromEnd) + '\t' + ancestry_state_last + '\n')
                    out.write(str(chrom) + '\t' + str(chromEnd) + '\t' + str(position) + '\t' + 'NA' + '\n')
                    chromStart = position
                    chromEnd = position

        ancestry_last = ancestry
        ancestry_state_last = ancestry_state
        chrom_last = chrom

    if chromEnd != 1:
        if ancestry_last < probability_threshold:
            out.write(str(dic_chrom_end[chrom_last]) + '\t' + 'N' + '\n')
        else:
            out.write(str(chrom_last) + '\t' + str(chromStart) + '\t' + str(dic_chrom_end[chrom_last]) + '\t' + ancestry_state_last + '\n')

    inp_posterior.close()
    out.close()


dic_chrom_end = {}
with open(length, 'r') as inp_length:
    header = inp_length.readline()
    for line in inp_length.readlines():
        chrom = line.strip().split('\t')[0]
        chrom_end = line.strip().split('\t')[2]
        dic_chrom_end[chrom] = chrom_end
    print(dic_chrom_end)


file_list = os.listdir(listdir)
for file in file_list:
    if file.endswith(".posterior"):
        ahmm2bed(file, probability_threshold)

