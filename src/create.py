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

def faDict2json(fa_dict):
    cys_dict = {}
    pattern_num = 1
    pattern_dict ={}
    for header in fa_dict:
        seq = fa_dict[header]
        cys_count = seq.count("C")
        if cys_count > 2: #NOTE make customizable
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

    return cys_dict
