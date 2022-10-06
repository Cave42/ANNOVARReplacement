import pandas as pd
import io
import os
from Bio import SeqIO
import re

import warnings

warnings.simplefilter(action='ignore', category=FutureWarning)

def read_vcf(path):
    with open(path, 'r') as f:
        lines = [l for l in f if not l.startswith('##')]
    return pd.read_csv(
        io.StringIO(''.join(lines)),
        dtype={'#CHROM': str, 'POS': int, 'ID': str, 'REF': str, 'ALT': str,
               'QUAL': str, 'FILTER': str, 'INFO': str},
        sep='\t'
    ).rename(columns={'#CHROM': 'CHROM'})

def translate(seq):
    table = {
        'ATA': 'I', 'ATC': 'I', 'ATT': 'I', 'ATG': 'M',
        'ACA': 'T', 'ACC': 'T', 'ACG': 'T', 'ACT': 'T',
        'AAC': 'N', 'AAT': 'N', 'AAA': 'K', 'AAG': 'K',
        'AGC': 'S', 'AGT': 'S', 'AGA': 'R', 'AGG': 'R',
        'CTA': 'L', 'CTC': 'L', 'CTG': 'L', 'CTT': 'L',
        'CCA': 'P', 'CCC': 'P', 'CCG': 'P', 'CCT': 'P',
        'CAC': 'H', 'CAT': 'H', 'CAA': 'Q', 'CAG': 'Q',
        'CGA': 'R', 'CGC': 'R', 'CGG': 'R', 'CGT': 'R',
        'GTA': 'V', 'GTC': 'V', 'GTG': 'V', 'GTT': 'V',
        'GCA': 'A', 'GCC': 'A', 'GCG': 'A', 'GCT': 'A',
        'GAC': 'D', 'GAT': 'D', 'GAA': 'E', 'GAG': 'E',
        'GGA': 'G', 'GGC': 'G', 'GGG': 'G', 'GGT': 'G',
        'TCA': 'S', 'TCC': 'S', 'TCG': 'S', 'TCT': 'S',
        'TTC': 'F', 'TTT': 'F', 'TTA': 'L', 'TTG': 'L',
        'TAC': 'Y', 'TAT': 'Y', 'TAA': 'STOP', 'TAG': 'STOP',
        'TGC': 'C', 'TGT': 'C', 'TGA': 'STOP', 'TGG': 'W',
    }
    protein = ""
    if len(seq) % 3 == 0:
        for i in range(0, len(seq), 3):
            codon = seq[i:i + 3]
            protein += table[codon]
    return protein

gff = pd.read_csv("C:/Users/theca/Downloads/lava_ref2.gff", sep='\t',header = None)

#variantFunction = pd.read_csv("/Users/administrator/Downloads/277D15_R1.trim.fastq.gz.fastq.avinput", sep="\t", header=None)
variantFunction = pd.read_csv("C:/Users/theca/Downloads/277D15_R1.trim.fastq.gz.fastq.avinput", sep="\t", header=None)

variantFunctionName = "C:/Users/theca/Downloads/277D15_R1.trim.fastq.gz.fastq.avinput"

variantFunctionName = os.path.basename(variantFunctionName)

#variantFunctionName.replace(".fastq.avinput", "")
variantFunctionName = re.sub('.fastq.avinput', '', str(variantFunctionName))

#print(variantFunctionName)

#fasta = ("C:/Users/theca/Downloads/consensus.fasta")

for fasta in SeqIO.parse("C:/Users/theca/Downloads/consensus.fasta", "fasta"):
 print("")
variantFunction.insert(0, "", "")
variantFunction.insert(0, " ", "")

numRegions = 0

for i in range(len(gff.index)):
    if(pd.isnull(gff.iloc[i,3]) != True):
        numRegions = numRegions +1
    else:
        break

test = 0

for j in range(numRegions):
    if(gff.iloc[j,2] == "CDS"):
        for k in range(len(variantFunction)):

            proteinName = gff.iloc[j, 8]
            proteinName = proteinName.replace("ID=CDS:", "gene:")
            sep = ';'
            proteinName = proteinName.split(sep, 1)[0]

            if(variantFunction.iloc[k,8] >= int(gff.iloc[j,3]) and variantFunction.iloc[k,8] <= int(gff.iloc[j,4])):
                if(variantFunction.iloc[k, 0] == ""):
                    variantFunction.iloc[k, 0] = "exonic"
                    variantFunction.iloc[k, 1] = proteinName
                elif (variantFunction.iloc[k, 0] != ""):
                    if (test == 0):
                        df = pd.DataFrame(variantFunction.loc[[k]])
                        df.iloc[0, 1] = proteinName
                    df2 = pd.DataFrame(variantFunction.loc[[k]])
                    df2.iloc[0, 1] = proteinName

                    if(test != 0):
                        df = df.append(df2)

                    test = test + 1

df = df.iloc[1: , :]

variantFunction = pd.concat([variantFunction, df])

variantFunction.to_csv('varriantfunction.csv',index = False,header = False)
variantFunction.to_csv(variantFunctionName + '.varriant_function', index = False, sep='\t', encoding='utf-8',header=False)

variantFunction.insert(0, "   ", "")
variantFunction.insert(0, "    ", "")

for l in range(len(variantFunction)):
    if((variantFunction.iloc[l,2]) == "exonic"):
        variantFunction.iloc[l, 0] = "line" + str(l+1)
        for m in range(len(gff.index)):
            if(gff.iloc[m, 2] == "gene"):
                if(variantFunction.iloc[l,5] >= gff.iloc[m,3] and variantFunction.iloc[l,5] <= gff.iloc[m,4]):
                    if((variantFunction.iloc[l,3] in gff.iloc[m,8])):

                        # adds deletions and insertions, detects if frameshift or not
                        if (len(variantFunction.iloc[l, 12]) > len(variantFunction.iloc[l, 13])):
                            if (len(variantFunction.iloc[l, 7]) % 3 == 0):
                                variantFunction.iloc[l, 1] = "nonframeshift deletion"
                            else:
                                variantFunction.iloc[l, 1] = "frameshift deletion"
                            #if ((int(variantFunction.iloc[l, 10]) - int(gff.iloc[m, 3])) % 3 == 0):
                                #if(translate(variantFunction.iloc[l,8] + fasta[int(variantFunction.iloc[l,10])] + fasta[int(variantFunction.iloc[l,10]+1)]) == "STOP"):

                                #if (translate(fasta[int(variantFunction.iloc[l, 10]) - 1] + fasta[int(variantFunction.iloc[l, 10])] + fasta[int(variantFunction.iloc[l, 10] + 1)]) == "STOP" and translate(variantFunction.iloc[l, 8] + fasta[int(variantFunction.iloc[l, 10])] + fasta[int(variantFunction.iloc[l, 10] + 1)]) != "STOP"):

                        if (len(variantFunction.iloc[l, 12]) < len(variantFunction.iloc[l, 13])):
                            if (len(variantFunction.iloc[l, 8]) % 3 == 0):
                                variantFunction.iloc[l, 1] = "nonframeshift insertion"
                            else:
                                variantFunction.iloc[l, 1] = "frameshift insertion"

                    #check if stopgain, not checking if stopgain in frameshift
                        if((int(variantFunction.iloc[l,10]) - int(gff.iloc[m,3]))%3 == 0):
                            if(variantFunction.iloc[l,8] !=  "-" and fasta[int(variantFunction.iloc[l,10])] !=  "-" and fasta[int(variantFunction.iloc[l,10]+1)] !=  "-"):
                                if(translate(variantFunction.iloc[l,8] + fasta[int(variantFunction.iloc[l,10])] + fasta[int(variantFunction.iloc[l,10]+1)]) == "STOP"):
                                    variantFunction.iloc[l, 1] = "stopgain"
                                if(translate(fasta[int(variantFunction.iloc[l,10])-1] + fasta[int(variantFunction.iloc[l,10])] + fasta[int(variantFunction.iloc[l,10]+1)]) == "STOP" and translate(variantFunction.iloc[l,8] + fasta[int(variantFunction.iloc[l,10])] + fasta[int(variantFunction.iloc[l,10]+1)]) != "STOP"):
                                    variantFunction.iloc[l, 1] = "stoploss"
                                if(translate(fasta[int(variantFunction.iloc[l,10])-1] + fasta[int(variantFunction.iloc[l,10])] + fasta[int(variantFunction.iloc[l,10]+1)]) == translate(variantFunction.iloc[l,8] + fasta[int(variantFunction.iloc[l,10])] + fasta[int(variantFunction.iloc[l,10]+1)])):
                                    variantFunction.iloc[l, 1] = "synonymous mutation"
                                if (translate(fasta[int(variantFunction.iloc[l, 10]) - 1] + fasta[int(variantFunction.iloc[l, 10])] + fasta[int(variantFunction.iloc[l, 10] + 1)]) != translate(variantFunction.iloc[l, 8] + fasta[int(variantFunction.iloc[l, 10])] + fasta[ int(variantFunction.iloc[l, 10] + 1)])):
                                    variantFunction.iloc[l, 1] = "nonsynonymous mutation"

                        if((int(variantFunction.iloc[l,10]) - int(gff.iloc[m,3]))%3 == 1):
                            if(variantFunction.iloc[l,8] !=  "-" and fasta[int(variantFunction.iloc[l,10]-2)] !=  "-" and fasta[int(variantFunction.iloc[l,10])] !=  "-"):
                                if(translate(fasta[int(variantFunction.iloc[l,10])-2] + variantFunction.iloc[l,8] + fasta[int(variantFunction.iloc[l,10])]) == "STOP"):
                                    variantFunction.iloc[l, 1] = "stopgain"
                                if(translate(fasta[int(variantFunction.iloc[l,10])-2] + fasta[int(variantFunction.iloc[l,10])-1] + fasta[int(variantFunction.iloc[l,10])]) == "STOP" and translate(fasta[int(variantFunction.iloc[l,10])-2] + variantFunction.iloc[l,8] + fasta[int(variantFunction.iloc[l,10])]) != "STOP"):
                                    variantFunction.iloc[l, 1] = "stoploss"
                                if (translate(fasta[int(variantFunction.iloc[l, 10]) - 2] + fasta[int(variantFunction.iloc[l, 10]) - 1] + fasta[int(variantFunction.iloc[l, 10])]) ==  translate(fasta[int(variantFunction.iloc[l, 10]) - 2] + variantFunction.iloc[l, 8] + fasta[int(variantFunction.iloc[l, 10])])):
                                    variantFunction.iloc[l, 1] =  "synonymous mutation"
                                if (translate(fasta[int(variantFunction.iloc[l, 10]) - 2] + fasta[int(variantFunction.iloc[l, 10]) - 1] + fasta[int(variantFunction.iloc[l, 10])]) !=  translate(fasta[int(variantFunction.iloc[l, 10]) - 2] + variantFunction.iloc[l, 8] + fasta[int(variantFunction.iloc[l, 10])])):
                                    variantFunction.iloc[l, 1] =  "nonsynonymous mutation"


                        if ((int(variantFunction.iloc[l, 10]) - int(gff.iloc[m, 3])) % 3 == 2):
                            if (variantFunction.iloc[l, 8] != "-" and fasta[int(variantFunction.iloc[l, 10] - 2)] != "-" and fasta[int(variantFunction.iloc[l, 10] - 3)] != "-"):
                                if (translate(fasta[int(variantFunction.iloc[l, 10]) - 3] +  fasta[int(variantFunction.iloc[l, 10])- 2] + variantFunction.iloc[l, 8]) == "STOP"):
                                    variantFunction.iloc[l, 1] = "stopgain"
                                if (translate(fasta[int(variantFunction.iloc[l, 10]) - 3] +  fasta[int(variantFunction.iloc[l, 10])- 2] + fasta[int(variantFunction.iloc[l, 10])- 1]) == "STOP" and translate(fasta[int(variantFunction.iloc[l, 10]) - 3] +  fasta[int(variantFunction.iloc[l, 10])- 2] + variantFunction.iloc[l, 8]) != "STOP"):
                                    variantFunction.iloc[l, 1] = "stoploss"
                                if (translate(fasta[int(variantFunction.iloc[l, 10]) - 3] +  fasta[int(variantFunction.iloc[l, 10])- 2] + fasta[int(variantFunction.iloc[l, 10])- 1]) == translate(fasta[int(variantFunction.iloc[l, 10]) - 3] +  fasta[int(variantFunction.iloc[l, 10])- 2] + variantFunction.iloc[l, 8])):
                                    variantFunction.iloc[l, 1] ="synonymous mutation"
                                if (translate(fasta[int(variantFunction.iloc[l, 10]) - 3] +  fasta[int(variantFunction.iloc[l, 10])- 2] + fasta[int(variantFunction.iloc[l, 10])- 1]) != translate(fasta[int(variantFunction.iloc[l, 10]) - 3] +  fasta[int(variantFunction.iloc[l, 10])- 2] + variantFunction.iloc[l, 8])):
                                    variantFunction.iloc[l, 1] = "nonsynonymous mutation"

                            if(len(variantFunction.iloc[l, 7]) == 2 and len(variantFunction.iloc[l, 8]) ==2):
                                variantFunction.iloc[l, 1] = "double AA mutation"

#indexNames = variantFunction[(variantFunction.columns[2] == 'exonic')].index
#variantFunction.drop(indexNames , inplace=True)



#variantFunction = variantFunction[variantFunction.columns[2].notnull()]

#variantFunction = variantFunction[variantFunction.iloc[:,2].notnull()]

variantFunction = variantFunction[variantFunction.iloc[:,2] == "exonic"]

variantFunction.drop(variantFunction.columns[2], axis=1,inplace=True)

variantFunction.to_csv('varriantfunction_exonic.csv',index = False,header = False)
variantFunction.to_csv(variantFunctionName + '.exonic_varriant_function', index = False, sep='\t', encoding='utf-8',header=False)

#print(gff.df.iloc[0,3])
#print(gff.df.iloc[0,4])


