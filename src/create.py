import re



def faDict2json(fa_dict):
    cys_dict = {}
    # pattern_num = 1
    # pattern_dict ={}
    for header in fa_dict:
        seq = fa_dict[header]
        cys_count = seq.count("C")
        if cys_count > 2 and not cys_count%2 : #NOTE make customizable
            cys_positions = [pos for pos, char in enumerate(seq) if char == "C"]
            cys_pattern = [x - (1+cys_positions[i - 1]) for i, x in enumerate(cys_positions)][1:]
            cys_seqs = [seq[cys_positions[i - 1]+1:x] for i,x in enumerate(cys_positions)][1:]
            cys_motif = ''.join(["-C-"+str(x - (1+cys_positions[i - 1])) for i, x in enumerate(cys_positions)][1:]+["-C"]).strip('-')

            general_pattern = cys_motif.replace("-0-","")
            general_pattern = general_pattern.replace("-1-","-X-")
            general_pattern = re.sub(r'[0-9]+', '', general_pattern)
            general_pattern = general_pattern.replace("--","-")

            if cys_count not in cys_dict:
                cys_dict[cys_count] = {}
            if general_pattern not in cys_dict[cys_count]:
                cys_dict[cys_count][general_pattern] = {}
            if cys_motif not in cys_dict[cys_count][general_pattern]:
                cys_dict[cys_count][general_pattern][cys_motif] = {}
                # if cys_motif not in pattern_dict:
                    # pattern_dict[cys_motif] = pattern_num
                    # pattern_num += 1
            cys_dict[cys_count][general_pattern][cys_motif][header] = cys_seqs
            #print(cys_count,cys_motif,header)

    return cys_dict
