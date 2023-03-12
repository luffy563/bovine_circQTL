import sys
import getopt


def usage():
    print('''
          Usage: python3 script.py [option] [patameter]
          -f/--fasta          input a fasta file
          -i/--idlist         input id list
          -o/--outfile        input outfile name
          -h/--help           show possible options
          ''')

opts,args = getopt.getopt(sys.argv[1:], 'hf:i:o:', ['help', 'fasta=', 'idlist=', 'outfile='])
for opt,val in opts:
    if opt == '-f' or opt == '--fasta':
        fasta = val 
    elif opt == '-i' or opt == '--idlist':
        idlist = val 
    elif opt == '-o' or opt == '--outfile':
        outfile = val 
    elif opt == '-h' or opt == '--help':
        usage()
        sys.exit(1)

outf = open(outfile, 'w')

dict = {}
with open(fasta, 'r') as fastaf:
    for line in fastaf:
        if line.startswith('>'):
            name = line.strip().split()[0][1:]
            dict[name] = ''
        else:
            dict[name] += line.replace('\n','')

with open(idlist, 'r') as listf:
    for row in listf:
        row = row.strip()
        for key in dict.keys():
            if row in key:
                outf.write('>' + key + '\n')
                outf.write(dict[key] + '\n')

outf.close()