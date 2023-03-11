#!/bin/python

input_file = './genome/bos_mature.fa'
output_file = './genome/all_bta_miR_family_info.txt'
with open(input_file,'r') as f:
    lines = f.readlines()
    for line in lines:
        if '>' in line:
            output_line = ''
            family = line.strip().split(' ')[-1]
            output_line += family
            output_line += '\t'
        else:
            output_seed = line[1:8]
            output_line += output_seed
            output_line += '\t'
            output_line += '9913\n'
            with open(output_file,'a+') as output:
                output.write(output_line)
                
# input_file = './mature.fa'
# output_file = './all_miR_info.txt'              
# with open(input_file,'r') as f:
    # lines = f.readlines()
    # for line in lines:
        # if '>' in line:
            # output_line = ''
            # family = line.strip().split(' ')[-1]
            # output_line += family
            # output_line += '\t'
        # else:
            # output_seed = line[1:8]
            # output_line += output_seed
            # output_line += '\t'
            # output_line += '9913\n'
            # with open(output_file,'a+') as output:
                # output.write(output_line)            