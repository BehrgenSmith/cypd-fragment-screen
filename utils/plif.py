import numpy as np 

PLIF_TYPES = {
    1: "DON SIDE",
    2: "ACC SIDE",
    3: "DON BACK",
    4: "ACC BACK",
    5: "DON SOL",
    6: "ACC SOL",
    7: "ION ATTR.",
    8: "SRF CON.",
    9: "MET CON.",
    10: "ARN CON.",
    11: "COV CON.",
    12: "HYD CON.",
    13: "CHG CON.",
    14: "PART CON.",
    15: "OTH CON.",
}

# PLIF Methods
def build_interaction_mapping(all_ruids, all_types):
    """
    all_ruids : iterable of all residue UIDs seen in the dataset
    all_types: iterable of all interaction type ints seen in the dataset

    Returns:
        mapping: dict[(ruid, type) -> index]
        dim: int, length of fingerprint vectors
    """
    unique_ruids = sorted(set(all_ruids))
    unique_types = sorted(set(all_types))

    mapping = {}
    idx = 0
    for r in unique_ruids:
        for t in unique_types:
            mapping[(r, t)] = idx
            idx += 1

    return mapping, idx


def make_interaction_fingerprint(ruids, types, mapping, dim):
    """
    ruids   : iterable of residue UIDs for this ligand
    types   : iterable of interaction type ints for this ligand
    mapping : dict[(ruid, type) -> index]
    dim     : length of fingerprint vector (from build_interaction_mapping)

    Returns:
        fp : 1D numpy array of shape (dim,), binary fingerprint
    """
    fp = np.zeros(dim, dtype=int)

    for r, t in zip(ruids, types):
        key = (r, t)
        if key in mapping:
            fp[mapping[key]] = 1
        else:
            raise KeyError

    return fp

def make_all_ruid_type(plifs):
    all_ruids = set()
    all_types = set()


    for plif in plifs: 
        if type(plif[0]) is list:
            if len(plif[0]) == 0:
                continue
            # print(plif[2])
            for i,ruid in enumerate(plif[2]):
                if ruid > 1000:
                    tmp = str(ruid)
                    ruid = int(tmp[2:])
                    plif[2][i] = ruid
                    # print(tmp,ruid)
                    
                all_ruids.add(ruid)
            for typ in plif[3]:
                all_types.add(typ)
        else:
            if plif[2] == '':
                continue
            if plif[2] > 1000:
                    tmp = str(plif[2])
                    plif[2] = int(tmp[2:])
            all_ruids.add(plif[2])
            all_types.add(plif[3])

    return all_ruids, all_types

def make_fps(plifs, mapping, length):
    fps = []
    skip_mask = []
    for plif in plifs:
        if type(plif[0]) is list:
            if len(plif[0]) == 0:
                skip_mask.append(True)
                continue
            fps.append(make_interaction_fingerprint(plif[2],plif[3],mapping,    length))
        else:
            fps.append(make_interaction_fingerprint([plif[2]],[plif[3]],    mapping, length))
        skip_mask.append(False)
    fps = np.array(fps)
    skip_mask = np.array(skip_mask,dtype=bool)
    print(f"Skipped : {len(np.where(skip_mask==1)[0])}")
    return fps, skip_mask
