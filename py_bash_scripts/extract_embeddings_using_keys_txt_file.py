import subprocess
import h5py
import numpy as np

def build_inmemory_index(keys_txt_path: str):
    """
    One‐time build: reads the 2‐column TSV (key \t index)
    Returns: { key_str: (int(index) - 1) }

    Avoid per‐line list allocations with partition (returns a 3‐tuple of (before, sep, after))
    Takes around 80 seconds to build.
    """
    idx_map = {}
    with open(keys_txt_path, 'r', encoding='utf-8') as f:
        for line in f:
            raw = line[:-1] # cheaper than rstrip('\n') but assumes every line ends in newline
            # partition at the first '\t'
            key, _, idx_str = raw.partition('\t')
            idx_map[key] = int(idx_str) - 1
    return idx_map

def _grep_index(keys_txt_path: str, query_id: str):
    """
    Returns the zero-based line index of query_id in keys_txt_path,
    or raises KeyError if not found.
    grep -n returns format like 123456:ID_QUERY
    """
    # -n : prefix each line with its 1-based line number
    proc = subprocess.run(
        ['grep', '-n', query_id, keys_txt_path],
        capture_output=True, text=True
    )
    if proc.returncode != 0 or not proc.stdout:
        # !r tells Python to call repr() on query_id rather than str()
        # e.g. with quotes in case there is newlines, or non-printable characters
        raise KeyError(f"{query_id!r} not found in keys file")
    # stdout like "12345:KEY\t6789\n"
    lineno = int(proc.stdout.split('\t', 1)[1].rstrip())
    return lineno - 1  # convert to 0-based

def get_embeddings_multi(
    h5_path: str,
    keys_txt_path: str,
    query_ids,
    threshold: int = 20
):
    """
    Returns a dict { key: embedding_array } for all found keys.
    
    - If 1 <= n < threshold: uses grep-per-key.
    - If n >= threshold: builds an in-memory index once.
    """
    # Normalise to a list
    if isinstance(query_ids, str):
        query_ids = [query_ids]
    total = len(query_ids)
    result = {}
    
    with h5py.File(h5_path, 'r') as f:
        emb_ds = f['embeddings']
        
        if total < threshold:
            # use grep for each
            for q in query_ids:
                try:
                    idx = _grep_index(keys_txt_path, q)
                    result[q] = emb_ds[idx, :]
                except KeyError:
                    # skip missing
                    pass
        else:
            # build in-memory index
            idx_map = build_inmemory_index(keys_txt_path)
            for q in query_ids:
                idx = idx_map.get(q)
                if idx is not None:
                    result[q] = emb_ds[idx, :]
                # else skip missing

    # Summary
    found = len(result)
    missing = total - found
    print(f"Total queried IDs : {total}")
    print(f"Found embeddings  : {found}")
    print(f"Missing embeddings: {missing}")
    if missing:
        missing_keys = set(query_ids) - set(result.keys())
        print("Missing keys:")
        for k in missing_keys:
            print(" ", k)
    
    return result

# example usages
if __name__ == "__main__":
    # Load embeddings and the corresponding text file of keys
    h5_path        = "/lisc/fileshare/globdb/globdb_r226/pLM_embeddings/globdb_r226_linclust_2plus_ProtT5_embeddings.h5"
    keys_txt_path  = "globdb_r226_linclust_2plus_ID_index.tsv"

    # For a single ID or list of IDs:
    query_ids      = ["BCRBG_15609___2917", "GCA_013288945___541", "MOTU40_023145___1271"]  
    embeddings_dict = get_embeddings_multi(
        h5_path,
        keys_txt_path,
        query_ids,
        threshold=20
    )
    # Or from a txt file of IDs (one per line, exactly as in GlobDB):
    with open('/lisc/scratch/dome/pullen/GlobDB/Testing/clustered_IDs_head10000.txt', 'r') as f:
         query_ids = [line.strip() for line in f if line.strip()]
    print(len(query_ids))
    embeddings_dict = get_embeddings_multi(
        h5_path,
        keys_txt_path,
        query_ids,
        threshold=20
    )
