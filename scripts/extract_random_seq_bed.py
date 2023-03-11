#!/bin/python3

# import optarse
import pandas as pd
import numpy as np

## config the attribute of random seq set
config = {
    "number":10,
    "meanlen":650.2962,
    "sdlen":573.4345,
    "minlen":6,
}

## randomly generate chrom and position according to the propability of chrom length and normal distribution respectively
genome_faindex = pd.read_csv('./genome/Bos_taurus.ARS-UCD1.2.dna.toplevel.fa.fai',sep='\t',header=None) # need fasta index file
genome_faindex.columns = ['CHROM','LENGTH','OFFSET','LINEBASES','LINEWIDTH']
random_len = np.random.normal(config['meanlen'],config['sdlen'],config['number'])
random_len[random_len<=config["minlen"]] = config["minlen"]

## set the chrom list and related propability list
chroms = list(range(1,30))+["X","MT"]
chrom_pos=list(genome_faindex['LENGTH'][0:31]/sum(genome_faindex['LENGTH'][0:31]))

## initialize output: random_seq_data (.bed)
random_seq_data = {
    "chrom":[],
    "start":[],
    "end":[],
    'name':[],
    'score':[],
    "strand":[],
}

## Main generate function
def main():
    n=1
    while n <= config['number']:
        # print(n)
        if n%100==0:
            print(n)
        chrom = np.random.choice(chroms,p=chrom_pos)
        max_len = int(genome_faindex[genome_faindex['CHROM']==chrom]['LENGTH'])
        # all_starts = list(range(max_len))
        ## accerlate the generate process 
        start = int(np.random.choice(np.array(range(int(max_len/10000))))*1e4+np.random.choice(list(range(max_len%10000))))
        # if random_len[n-1]<=6:
        end = start + int(random_len[n-1])
        random_seq_data['chrom'].append(chrom)
        random_seq_data['start'].append(start)
        random_seq_data['end'].append(end)
        random_seq_data['name'].append(".")
        random_seq_data['score'].append(".")
        random_seq_data['strand'].append("*")
        n+=1
    # convert to dataframe
    random_seq_dataframe = pd.DataFrame(random_seq_data)
    # export file
    random_seq_dataframe.to_csv("./random_seq.bed",sep="\t",index=False)
    
if __name__ == "__main__":
    main()