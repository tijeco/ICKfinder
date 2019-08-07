import re
import json
import pandas as pd

def annotateFastaDict(fa_dict, json_file,motif=True):
    out_dict = {"header":[],"pattern":[],"Cys_num":[],"hits":[]}
    # print(fa_dict,json_file)
    with open(json_file) as f:
        json_dict = json.load(f)
    # print(json_dict)
    for header in fa_dict:
        seq = fa_dict[header]
        cys_count = seq.count("C")
        if cys_count > 2: #NOTE make customizable
            cys_positions = [pos for pos, char in enumerate(seq) if char == "C"]
            cys_pattern = [x - (1+cys_positions[i - 1]) for i, x in enumerate(cys_positions)][1:]
            cys_seqs = [seq[cys_positions[i - 1]+1:x] for i,x in enumerate(cys_positions)][1:]
            cys_motif = ''.join(["-C-"+str(x - (1+cys_positions[i - 1])) for i, x in enumerate(cys_positions)][1:]+["-C"]).strip('-')

            general_pattern = cys_motif.replace("-0-","")
            general_pattern = general_pattern.replace("-1-","-X-")
            general_pattern = re.sub(r'[0-9]+', '', general_pattern)
            general_pattern = general_pattern.replace("--","-")
            cys_count_str = str(cys_count)
            # print(cys_count_str in json_dict)
            if cys_count_str in json_dict:
                if motif:
                    if general_pattern in json_dict[cys_count_str]:
                        hits = []

                        for key in json_dict[cys_count_str][general_pattern]:
                            for sp in json_dict[cys_count_str][general_pattern][key]:
                                hits.append(sp)
                        hits = ';'.join(hits)
                        # print(header,general_pattern,hits)
                        out_dict["header"].append(header)
                        out_dict["pattern"].append(general_pattern)
                        out_dict["Cys_num"].append(cys_count_str)
                        out_dict["hits"].append(hits)
                        # keys = [key for key in json_dict[cys_count_str][general_pattern]]]
                        # hits = ';'.join([key for key in json_dict[cys_count_str][general_pattern]])
                else:
                    if general_pattern in json_dict[cys_count_str]:
                        if cys_motif in json_dict[cys_count_str][general_pattern]:
                            hits = ';'.join(json_dict[cys_count_str][general_pattern][cys_motif].keys())
                            # print(header,cys_motif,hits)
                            out_dict["header"].append(header)
                            out_dict["pattern"].append(cys_motif)
                            out_dict["Cys_num"].append(cys_count_str)
                            out_dict["hits"].append(hits)

    out_pd = pd.DataFrame.from_dict(out_dict)
    return out_pd
                # else:
                #     print(general_pattern, "not in", json_dict[cys_count_str].keys())
            # print(header,general_pattern)
