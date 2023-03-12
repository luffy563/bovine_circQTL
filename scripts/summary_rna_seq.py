import os
import venn as v
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import random


# # Define the color vector
# def randomcolor(n):
    # colorArr = ['1','2','3','4','5','6','7','8','9','A','B','C','D','E','F']
    # colors=[]
    # for i in range(n):
        # color=''
        # for j in range(6):
            # color += colorArr[random.randint(0,14)]
        # color = '#' + color
        # colors.append(color)
    # return colors

# colors = randomcolor(len(softlist))

# Import raw predict results of different circRNA detection software
def import_each_raw(line, i):
    if line == 'CIRCexplorer2':
        path = './' + line + '/' + i + '/' + 'circ_candidates_convert.bed'
        print(path)
        data = pd.read_table(path, names=['chrom', 'start', 'end', 'name', 'score', 'strand', 'unknown'])
        data_name = data['name'].str.split('/', expand=True)
        data[['name', 'n_reads', 'score']] = data_name
        data = data.astype(str)
        data['readcounts'] = pd.Series([int(x) for x in data['n_reads']])
        data = data[['chrom', 'start', 'end', 'strand','readcounts']]
        data = data[data['readcounts'] >= 2]

    elif line == 'CircRNAfinder':
        path = './' + line + '/' + i + '/' + 'circ_candidates_convert.bed'
        print(path)
        data = pd.read_table(path, names=['chrom', 'start', 'end', 'name', 'readcounts', 'strand'], index_col=False)
        data['readcounts'] = data['readcounts'].fillna(0)
        data = data.astype(str)
        data['readcounts'] = pd.Series([int(x) for x in data['readcounts']])
        data = data[['chrom', 'start', 'end', 'strand','readcounts']]
        data = data[data['readcounts'] >= 2]


    elif line == 'CIRI2':
        path = './' + line + '/' + i + '/' + 'circ_candidates_convert.bed'
        print(path)
        data = pd.read_table(path).iloc[:, 0:10]
        data.columns = ['chrom', 'start', 'end', 'name', 'readcounts', 'non_junction_reads', 'junction_ratio',
                        'circRNA_type','gene','strand']
        data_name = data['name'].str.split('/', expand=True)
        data[['name', 'n_reads']] = data_name
        data = data.astype(str)
        data['readcounts'] = pd.Series([int(x) for x in data['n_reads']])
        data = data[data['readcounts'] >= 2]

    elif line == 'CircMarker':
        path = './' + line + '/' + i + '/' + 'circ_candidates_convert.bed'
        print(path)
        data = pd.read_table(path, names=['chrom', 'start', 'end', 'name', 'strand', 'type', 'readcounts'],
                             index_col=False)
        data['readcounts'] = data['readcounts'].fillna(0)
        data = data.astype(str)
        data['readcounts'] = pd.Series([int(x) for x in data['readcounts']])
        data = data[['chrom', 'start', 'end','strand', 'readcounts']]
        data = data[data['readcounts'] >= 2]
    return data

# Read RNA-seq dataset
def import_data(working_dir,all_rna_seq_meta_file,pe_rna_seq_meta_file,softlist,quant):
    norm,hypo,yLabel = [],[],[]
    total_dataset = []
    for line,color in zip(softlist, colors):
        line = line.replace('\n', '')
        dataset_score = []
        if quant == True:
            ## Accurate quantity by CIRIquant
            # data_list = open('PE.txt')
            sra_list = pd.read_table(pe_rna_seq_meta_file, sep='\t', header=None)
            data_list_total = open(all_rna_seq_meta_file)
            sra_list_total = pd.read_table(all_rna_seq_meta_file, sep='\t', header=None)
            # print(len(sra_list))
            for i, index in zip(data_list_total, range(len(sra_list_total))):
                i = i.replace('\n', '')
                if i in np.array(sra_list).ravel().tolist():
                    path = os.path.join(working_dir, line, i, i + '_quant.gtf')
                    print("Loading circRNA detection calibrated results from ".format(path))
                    data = pd.read_table(path, names=['chrom', 'soft', 'circRNA', 'start', 'end', 'score', 'strand',
                                                      'unknown',
                                                      'info'], comment='#')
                    data_name = data['info'].str.split(';', expand=True)
                    data[['n_reads']] = (data_name.iloc[:, 2].str.split(' ', expand=True)).iloc[:, 2]

                    data = data[['chrom', 'start', 'end', 'strand','n_reads']]
                    data['start'] = data['start'] - 1
                    data = data.astype(str)
                    data['chrom'] = pd.Series([('chr' + x) for x in data['chrom']])
                    data['readcounts'] = pd.Series([float(x) for x in data['n_reads']])
                    data = data[data['readcounts'] >= 2]
                else:
                    data = import_each_raw(line,i)

                dataset_score.append(data)
                if index == 0:
                    # total_data = data
                    total_data = data[['chrom', 'start', 'end','strand']]

                else:
                    total_data = pd.merge(total_data, data[['chrom', 'start', 'end','strand']], how='outer',
                                          on=['chrom', 'start', 'end','strand'])
                    total_data = total_data.drop_duplicates(subset=['chrom', 'start', 'end','strand'], keep='first')

        else:
            data_list = open(all_rna_seq_meta_file)
            sra_list = pd.read_table(all_rna_seq_meta_file,sep='\t',header=None)
            # print(len(sra_list))
            # assert(len(sra_list) == 95)

            for i,index in zip(data_list,range(len(sra_list))):
                i = i.replace('\n', '')
                data = import_each_raw(line,i)
                if index == 0:
                    # total_data = data
                    total_data = data[['chrom', 'start', 'end','strand']]

                else:
                    total_data = pd.merge(total_data, data[['chrom', 'start', 'end','strand']], how='outer', on=['chrom', 'start', 'end','strand'])
                    total_data = total_data.drop_duplicates(subset=['chrom', 'start', 'end','strand'],keep='first')

                dataset_score.append(data)


        total_data = total_data.drop_duplicates(subset=['chrom', 'start', 'end','strand'],keep='first')
        total_data['index'] = total_data.index
        # add the readcounts of each sample
        data_list = open(all_rna_seq_meta_file)
        for i, j in zip(data_list, range(len(dataset_score))):
            i = i.replace('\n', '')
            total_data[i]=0
            temp = dataset_score[j]
            temp['index'] = temp.index
            temp_data = pd.merge(total_data, temp, how='inner', on=['chrom', 'start', 'end','strand'])
            idx = temp_data['index_x']
            total_data.loc[idx, i] = temp_data['readcounts'].to_list()
        total_data=total_data.fillna(0)
        total_dataset.append(total_data)
        yLabel.append(line)
    return total_dataset


def filtration(total_dataset,filt,all_rna_seq_meta_file,softlist):
    sra_list_total = pd.read_table(all_rna_seq_meta_file, sep='\t', header=None)
    sra_list = sra_list_total.iloc[:,0].to_list()
    dataset_filt = []
    for i in range(len(softlist)):
        if filt == False:
            # No filtration
            total_dataset_filt = total_dataset[i]
        else:
            # Filtration
            idx = (total_dataset[i].loc[:,sra_list].sum(axis=1) >= 10).to_list()
            total_dataset_filt = total_dataset[i][idx]

        dataset_filt.append(total_dataset_filt)
        if i == 0:
            total = total_dataset_filt[['chrom', 'start', 'end','strand']]
        else:
            total = pd.merge(total,total_dataset_filt[['chrom', 'start', 'end','strand']],on=['chrom', 'start', 'end','strand'],how='outer')

    total = total.drop_duplicates(subset=['chrom', 'start', 'end','strand'],keep='first')
    total['ID'] = range(len(total))
    norm =[]
    for i in range(len(softlist)):
        norm.append(pd.merge(total, dataset_filt[i],on=['chrom', 'start', 'end','strand'],how='inner'))

    # normsets = [set(norm[0]['ID']), set(norm[1]['ID']), set(norm[2]['ID'])]
    normsets=[]
    for i in range(len(norm)):
        normsets.append(set(norm[i]['ID']))
    return normsets,dataset_filt,total,norm



# Venn plot of dataset identified
def intersect(total,normsets,softlist):
    plt.figure(figsize=(12, 8), dpi=300)
    #
    d={}
    d.fromkeys(d.fromkeys(softlist))
    for line, i in zip(softlist, range(len(softlist))):
        line = line.replace('\n', '')
        print(i)
        print(line)
        d[line] = normsets[i]

    out = v.venn(d, fontsize=20)
    plt.show()
    # Get the intersection of high confidence cicRNA in different softwares
    soft_ID = []
    for l in range(len(softlist)):
        for j in range(l+1, len(softlist)):
            key = softlist[l]
            key_1 = softlist[j]
            soft_ID.append(d[key] & d[key_1])
    for i in range(len(soft_ID)):
        if i == 0:
            outer_soft_ID = soft_ID[i]
        else:
            outer_soft_ID = outer_soft_ID | soft_ID[i]

    inter_loc = total.merge(pd.DataFrame(list(outer_soft_ID),columns=['ID']),on='ID',how='inner')[['chrom', 'start', 'end','strand']].drop_duplicates(subset=['chrom', 'start', 'end','strand'],keep='first')
    return inter_loc

# Get the CPM of circRNA from each software
def normlizatin(dataset_filt,inter_loc):
    inter_loc=inter_loc.astype(str)
    # Normalization of the read counts for different software
    for i in range(len(dataset_filt)):
        norm_data = pd.merge(dataset_filt[i], inter_loc, on=['chrom', 'start', 'end','strand'], how='outer').drop_duplicates(subset=['chrom', 'start', 'end','strand'],keep='first')
        norm_data = pd.merge(norm_data, inter_loc, on=['chrom', 'start', 'end','strand'], how='inner').drop_duplicates(subset=['chrom', 'start', 'end','strand'],keep='first')
        norm_data = norm_data.fillna(0)
        name = inter_loc['chrom'] + ':' + inter_loc['start'] + '|' + inter_loc['end']
        norm_name = norm_data['chrom'] + ':' + norm_data['start'] + "|" + norm_data['end']
        norm_data.index = pd.Series(norm_name.to_list())
        norm_data = norm_data.loc[name.to_list()]
        if i == 0:
            # array = np.array(norm_data.loc[:, sra_list])
            array = np.array(norm_data.loc[:, sra_list])
        else:
            array = (array + np.array(norm_data.loc[:, sra_list])/ (i + 1))
        # array = array / (i + 1)
        norm_info = norm_data.iloc[:, 0:4]
        norm_info.reset_index(drop=True, inplace=True)
    array = (array*10**6)/np.sum(array,axis=0)
    log_array = np.log10(np.nansum(array,axis=0)+1)
    return array,log_array,norm_info


# Import annotation results from annotated circRNA
def import_anno(working_dir,all_rna_seq_meta_file,softlist):
    total_dataset_anno = []
    for line,color in zip(softlist, colors):
        line = line.replace('\n', '')
        dataset_score = []
        data_list = open(all_rna_seq_meta_file)
        sra_list = pd.read_table(all_rna_seq_meta_file,sep='\t',header=None)
        # print(len(sra_list))
        # assert(len(sra_list) == 94)

        for i,index in zip(data_list,range(len(sra_list))):
            i = i.replace('\n', '')
            path = os.path.join(working_dir, line, i + '/circularRNA_known.txt')
            print("Loading annotated result from ".format(path))
            if os.stat(path).st_size == 0:
                data = pd.DataFrame(columns=['chrom', 'start', 'end', 'name', 'none', 'strand', 'startExon_1', 'startExon', 'unknown',
                                'rea', 'Exons', 'Endexons', 'readcounts', 'type'])
                data['readcounts'] = 0
                dataset_score.append(data)
                continue
            data = pd.read_table(path).iloc[:, :14]
            data.columns = ['chrom', 'start', 'end', 'name', 'none', 'strand', 'startExon_1', 'startExon', 'unknown',
                            'rea', 'Exons', 'Endexons', 'readcounts', 'type']
            data['readcounts'] = data['readcounts'].fillna(0)
            data = data.astype(str)
            data['readcounts'] = pd.Series([int(x) for x in data['readcounts']])
            data = data[['chrom', 'start', 'end', 'strand','readcounts','type']]
            data = data[data['readcounts'] >= 2]
            dataset_score.append(data)
            if index == 0:
                # total_data = data
                total_data_anno = data[['chrom', 'start', 'end','strand','type']]

            else:
                total_data_anno = pd.merge(total_data_anno, data[['chrom', 'start', 'end','strand','type']], how='outer', on=['chrom', 'start', 'end','strand','type'])
                total_data_anno = total_data_anno.drop_duplicates(subset=['chrom', 'start', 'end','strand','type'], keep='first')

        total_data_anno = total_data_anno.drop_duplicates(subset=['chrom', 'start', 'end','strand','type'], keep='first')
        total_data_anno['index'] = total_data_anno.index
        # add the readcounts of each sample
        data_list = open(all_rna_seq_meta_file)
        for i, j in zip(data_list, range(len(dataset_score))):
            i = i.replace('\n', '')
            total_data_anno[i] = 0
            temp = dataset_score[j]
            temp['index'] = temp.index
            temp_data = pd.merge(total_data_anno, temp, how='inner', on=['chrom', 'start', 'end','strand','type'])
            idx = temp_data['index_x']
            total_data_anno.loc[idx, i] = temp_data['readcounts'].to_list()

        total_data_anno = total_data_anno.fillna(0)
        total_dataset_anno.append(total_data_anno)
    return total_dataset_anno

# Filtration of annotated circRNA results
def filtration_anno(total_dataset_anno,filt,softlist):
    # Filtration of circRNAs annotated
    dataset_anno = []
    for i in range(len(softlist)):
        if filt == False:
            # No filtration
            anno_dataset_filt = total_dataset_anno[i]
        else:
            # Filtration by total readscounts >= 10 in all samples
            idx = (total_dataset_anno[i].loc[:,sra_list].sum(axis=1) >= 10).to_list()
            anno_dataset_filt = total_dataset_anno[i][idx]

        dataset_anno.append(anno_dataset_filt)
        if i == 0:
            total_anno = anno_dataset_filt[['chrom', 'start', 'end','strand','type']]
        else:
            total_anno = pd.merge(total_anno,anno_dataset_filt[['chrom', 'start', 'end','strand','type']],on=['chrom', 'start', 'end','strand','type'],how='outer')
    # Define the ID of circRNA annotated and filted
    total_anno = total_anno.drop_duplicates(subset=['chrom', 'start', 'end','strand','type'],keep='first')
    total_anno['ID'] = range(len(total_anno))
    anno =[]
    for i in range(len(softlist)):
        # inter_anno = pd.merge(dataset_anno[i],inter_loc,on=['chrom', 'start', 'end'],how='inner')
        anno.append(pd.merge(total_anno, dataset_anno[i],on=['chrom', 'start', 'end','strand','type'],how='inner'))

    # Combine all annotated candidates by each soft
    annosets=[]
    for i in range(len(anno)):
        annosets.append(set(anno[i]['ID']))
    return dataset_anno,annosets,total_anno,anno


def intersect_anno(total_anno,annosets,softlist):
    # Annotated intersection
    d={}
    d.fromkeys(d.fromkeys(softlist))
    for line, i in zip(softlist, range(len(softlist))):
        line = line.replace('\n', '')
        print(i)
        print(line)
        d[line] = annosets[i]
    # Venn plot of circRNA annotated
    plt.figure(figsize=(12,8),dpi=300)
    out = v.venn(d, fontsize=20)
    plt.show()

    # Get the intersection of high confidence cicRNA in different software
    soft_ID = []
    for l in range(len(softlist)):
        for j in range(l+1, len(softlist)):
            print(l,j)
            key = softlist[l]
            key_1 = softlist[j]
            soft_ID.append(d[key] & d[key_1])
    for i in range(len(soft_ID)):
        if i == 0:
            outer_soft_ID = soft_ID[i]
        else:
            outer_soft_ID = outer_soft_ID.union(soft_ID[i])

    inter_loc = total_anno.merge(pd.DataFrame(list(outer_soft_ID),columns=['ID']),on='ID',how='inner')[['chrom', 'start', 'end','strand']].drop_duplicates(subset=['chrom', 'start', 'end','strand'],keep='first')
    return inter_loc

# # The cumulative curves of expression for CircRNA between before and after quantity calibration 
# def plot_cumulate_curve(dataSets,label,colors):
    # plt.figure(figsize=(12,8),dpi=300)
    # temp = []
    # for set in dataSets:
        # temp2 = []
        # for item in set:
            # if item != '':
                # temp2.append(float(item))
        # temp2.sort()
        # temp.append(temp2)
    # dataSets = temp

    # for set,label,color in zip(dataSets,label,colors):

        # plotDataset = [[], []]
        # count = len(set)
        # for i in range(count):
            # # plotDataset[0].append(float(set[i]))
            # plotDataset[0].append(float(set[i]))
            # plotDataset[1].append((i + 1) / count)
        # print(plotDataset)
        # plt.step(plotDataset[0], plotDataset[1], '-', color=color,linewidth=2,label=label)
        # # plt.legend()
        # plt.xticks(fontsize=18)
        # plt.yticks(fontsize=18)
        # plt.xlabel('${log_2}$(CPM+1)', fontdict={'fontsize': 25})
        # plt.ylabel('Frequency of circRNAs', fontdict={'fontsize': 25})
        # plt.legend(fontsize=18)
    # plt.savefig('./results/fig2-5C.pdf')
    # plt.show()
        
# Get the CPM of annotated circRNA results
def nomalization_anno(dataset_anno,inter_loc):
    # Nomalization of annotated circRNA with filtration
    for i in range(len(dataset_anno)):
        anno_data = pd.merge(dataset_anno[i], inter_loc, on=['chrom', 'start', 'end','strand'], how='outer').drop_duplicates(
            subset=['chrom', 'start', 'end','strand'], keep='first')
        anno_data = pd.merge(anno_data, inter_loc, on=['chrom', 'start', 'end','strand'], how='inner').drop_duplicates(
            subset=['chrom', 'start', 'end','strand'], keep='first')
        anno_data = anno_data.fillna(0)
        name=inter_loc['chrom']+':'+inter_loc['start']+ '|' +inter_loc['end']
        anno_name = anno_data['chrom'] + ':' + anno_data['start'] + "|" + anno_data['end']
        anno_data.index = pd.Series(anno_name.to_list())
        anno_data = anno_data.loc[name.to_list()]
        if i == 0:
            anno_array = np.array(anno_data.loc[:, sra_list])
        else:
            anno_array = anno_array + np.array(anno_data.loc[:, sra_list])/ (i + 1)
        # anno_array = anno_array / (i + 1)
        anno_info = anno_data.iloc[:, 0:4]
        anno_info.reset_index(drop=True, inplace=True)
    anno_array = (anno_array*10**6)/np.sum(anno_array,axis=0)
    return anno_array,anno_info

# Importing the basic info of circRNAs detection software
def get_software_list(soft_list_file):
    softlist,typelist=[],[]
    with open(soft_list_file, 'r') as f:
        line = f.readline().replace('\n', '')
        while line != '':
            softlist.append(line)
            if line == 'CircMarker' or line == 'CircDBG':
                typelist.append('k-mer')
            else:
                typelist.append('readsmapping')
            line = f.readline().replace('\n', '')
    return softlist,typelist
        
if __name__ == "__main__":
    # Parameters setting
    # CIRIquant
    compare=True
    quant=False
    filt=True
    working_dir = './'
    all_rna_seq_meta_file = 'SRA_list.txt'
    pe_rna_seq_meta_file = 'PE.txt'
    softlist,typelist = get_software_list('software_list_quant.txt')

    if compare == True:
        total_dataset = import_data(working_dir,all_rna_seq_meta_file,pe_rna_seq_meta_file,softlist,quant=False)
        total_dataset_quant = import_data(working_dir,all_rna_seq_meta_file,pe_rna_seq_meta_file,softlist,quant=True)

        # Filtration on normsets 
        normsets, dataset_filt, total,norm = filtration(total_dataset,filt,all_rna_seq_meta_file,softlist)
        quantsets, quantdataset_filt, quanttotal,anno = filtration(total_dataset_quant, filt,all_rna_seq_meta_file,softlist)
        
        
        inter_loc = intersect(total, normsets,softlist)
        quantinter_loc = intersect(quanttotal, quantsets,softlist)

        array,log_array,norm_info = normlizatin(dataset_filt, inter_loc)
        anno_array,log_anno_array,anno_info = normlizatin(quantdataset_filt, quantinter_loc)
    else:
        
        total_dataset = import_data(working_dir,all_rna_seq_meta_file,pe_rna_seq_meta_file,softlist,quant=False)
        # Import annotated datasets
        total_dataset_anno = import_anno(working_dir,all_rna_seq_meta_file,softlist)
        
        # Filtration on normsets
        normsets, dataset_filt, total, norm = filtration(total_dataset, filt,all_rna_seq_meta_file,softlist)
        dataset_anno,annosets,total_anno,anno=filtration_anno(total_dataset_anno,filt,all_rna_seq_meta_file,softlist)
        
        ###################
        #   Figure S1A-D  #
        ###################
        # Venn plot of circRNA detection results from different software
        inter_loc = intersect(total, normsets,softlist)
        anno_inter_loc = intersect_anno(total_anno, annosets,softlist)
        
        # Nomalization of circRNA quantity results
        array, log_array, norm_info = normlizatin(dataset_filt, inter_loc)
        anno_array, anno_info = nomalization_anno(dataset_anno, anno_inter_loc)

    # log2-scaled expression (log2(mean+1))
    log_array = np.log2(np.nanmean(array,axis=0)+1)
    log_anno_array = np.log2(np.nanmean(anno_array,axis=0)+1)

    log_anno_array_total = np.log2(anno_array+1)
    log_anno_array_total = pd.DataFrame(log_anno_array_total)
    # anno_array.columns = anno[0].columns[6:101]
    if compare == True:
        log_anno_array_total.columns = anno[0].columns[6:]
    else:
        log_anno_array_total.columns = anno[0].columns[7:]
    anno_reads_matrix = pd.concat([log_anno_array_total,anno_info],axis=1)
    # Writing to the information of circRNA annotated
    anno_info = pd.DataFrame(anno_info)
    circ_name = []
    for i in range(len(anno_info.iloc[:,0])):
        name = anno_info.iloc[:,0][i]+':'+str(anno_info.iloc[:,1][i])+'-'+str(anno_info.iloc[:,2][i])
        circ_name.append(name)

    anno_info_convert = pd.concat([pd.Series(circ_name),anno_info],axis=1)
    anno_info_convert.to_csv('anno_info_strand.txt',sep='\t',index=False)
    # Writing to the expression matrix of circRNA annotated
    anno_reads_matrix.to_csv('log_anno_array_hfilt_3soft_anno_strand.csv',index=False)

    
    #################
    #   Figure S1E  #
    #################
    # Pie plot of different subtypes in interloc from CIRI2
    data = pd.merge(total_dataset_anno[1],inter_loc,how='inner',on=['chrom','start','end'])
    a=np.sum(data['type'] == 'circRNA')/len(data['type'])
    b=np.sum(data['type'] == 'intron')/len(data['type'])
    c=np.sum(data['type'] == 'intergenic_region')/len(data['type'])
    fracs = [a,b,c]
    labels = ['exon_circRNA','ciRNA','intergenic_region']
    fig=plt.figure(figsize=(12,8),dpi = 300)
    plt.pie(fracs, labels=labels, autopct='%.1f%%', shadow=False,textprops={'fontsize':25})
    plt.show()
