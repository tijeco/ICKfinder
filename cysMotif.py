import re
import json
import pandas as pd
import sys

#test_file = "test/SPIDER.knottin.fa"
def fa2dict(fa):
    fasta_dict = {}
    with open(fa) as f:
        header = ""
        for line in f:
            if line.strip() != '':
                if line[0] == ">":
                    header = line.strip().strip(">")
                    fasta_dict[header] = ""
                else:
                    fasta_dict[header] += line.strip()
    return fasta_dict

def findICK(fa_dict,json=True):
    cys_dict = {}
    pattern_num = 1
    pattern_dict ={}
    cys_pd_dict = {"count":[],"pattern_num":[],"pattern":[],"header":[]}
    for header in fa_dict:
        seq = fa_dict[header]
        cys_count = seq.count("C")
        if cys_count != 0:
            cys_positions = [pos for pos, char in enumerate(seq) if char == "C"]
            cys_pattern = [x - (1+cys_positions[i - 1]) for i, x in enumerate(cys_positions)][1:]
            cys_seqs = [seq[cys_positions[i - 1]+1:x] for i,x in enumerate(cys_positions)][1:]
            cys_motif = ''.join(["-C-"+str(x - (1+cys_positions[i - 1])) for i, x in enumerate(cys_positions)][1:]+["-C"]).strip('-')
            if cys_count not in cys_dict:
                cys_dict[cys_count] = {}
            if cys_motif not in cys_dict[cys_count]:
                cys_dict[cys_count][cys_motif] = {}
                if cys_motif not in pattern_dict:
                    pattern_dict[cys_motif] = pattern_num
                    pattern_num += 1
            cys_dict[cys_count][cys_motif][header] = cys_seqs
            #print(cys_count,cys_motif,header)
            cys_pd_dict["count"].append(cys_count)
            cys_pd_dict["pattern"].append(cys_motif)
            cys_pd_dict["pattern_num"].append(pattern_dict[cys_motif])
            cys_pd_dict["header"].append(header)
    if json:
        return cys_dict
    else:
        return pd.DataFrame.from_dict(cys_pd_dict)


#spider_knottin_fa = fa2dict("test/SPIDER.knottin.fa")
# spider_knottin_pd = findICK(spider_knottin_fa,json=False)

# try:
#     input_file = "test/P.nigriventer.ick.fa" # replace with sys input
# except:
#     input_file = sys.argv[1]

json_out = True # replace with sys input

try:
    input_fa = fa2dict("test/P.nigriventer.ick.fa")
except:
    input_fa = fa2dict(sys.argv[1])
if json_out:
    input_dict = findICK(input_fa)
    for cys_num in input_dict:
        for pattern in input_dict[cys_num]:
            general_pattern = pattern.replace("-0-","")
            general_pattern = general_pattern.replace("-1-","-X-")
            general_pattern = re.sub(r'[0-9]+', '', general_pattern)
            general_pattern = general_pattern.replace("--","-")

            for header in input_dict[cys_num][pattern]:
                print(header+","+general_pattern+","+pattern)

else:
    input_pd = findICK(input_fa,json=False)
# print(spider_knottin_pd)
# print(P_nigriventer_dict)


        # print(pattern,len(P_nigriventer_dict[cys_num][pattern]))
