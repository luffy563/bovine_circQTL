# /home/liuhongfei/miniconda3/envs/py3/bin/python

import pandas as pd
import numpy as np
import os, time
from multiprocessing.dummy import Pool as ThreadPool

circRNA = './ABS_events.txt'
introns = './introns.gtf'

intron_info_convert = pd.read_csv('./alu_info.bed', sep='\t')
def extract_pos(pos, t):
    '''
    pos: tuple, (start,end)

    pos[0]
    Returns:
    start
    end
    '''
    
    if pos == 'nan':
        return 'NA'
    else:
        start = pos.split(', ')
        return start


p = ThreadPool(10)
result_start = []
result_end = []
idx = (intron_info_convert['Alu_number'] != '0.0') & (intron_info_convert['Alu_number'] !=0)
Alu_start = intron_info_convert[idx]['Alu_start']
Alu_end = intron_info_convert[idx]['Alu_end']
chrom = intron_info_convert[idx]['chrom']
strand = intron_info_convert[idx]['strand']
names = intron_info_convert[idx]['Alu_name']
n,c,s = [],[],[]
for pos,pos_1,ch,st,name in zip(Alu_start,Alu_end,chrom,strand,names):
    print(pos)
    
    res = p.apply_async(func=extract_pos, args=(pos, 'start'))
    res_1 = p.apply_async(func=extract_pos, args=(pos_1, 'start'))
    res_2 = p.apply_async(func=extract_pos, args=(name, 'start'))
    for j in range(len(res.get())):
        c.append(ch[3:]);s.append(st);n.append(res_2.get()[j])
        result_start.append(res.get()[j])
        result_end.append(res_1.get()[j])
p.close()
p.join()

# ouput the alu bed
alu_info = {'chrom':c,'start':result_start,'end':result_end,'name':n,'strand':s}
alu_info = pd.DataFrame(alu_info).loc[:,['chrom','start','end','name','strand']]
alu_info = alu_info.drop_duplicates()

alu_info.to_csv('./alu_info_sep.bed',sep='\t',index=0,columns=None)