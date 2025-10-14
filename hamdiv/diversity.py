import numpy as np
from rdkit import Chem, DataStructs
from rdkit.Chem import AllChem
from rdkit.Chem.Scaffolds import MurckoScaffold
import six
import sys
sys.modules['sklearn.externals.six'] = six
import random
from tqdm import tqdm

import networkx as nx
from python_tsp.exact import solve_tsp_dynamic_programming
from python_tsp.heuristics import solve_tsp_local_search
from .utils import identify_functional_groups, GetRingSystems


def dist_array(smiles=None, mols=None):
    """
    Calculate distance matrix for molecules using Tanimoto distance of ECFPs.
    
    Args:
        smiles: List of SMILES strings
        mols: List of RDKit molecule objects
        
    Returns:
        numpy.ndarray: Distance matrix
    """
    if mols == None:
        l = len(smiles)
        mols = [Chem.MolFromSmiles(s) for s in smiles]
    else:
        l = len(mols)
    '''
    You can replace the Tanimoto distances of ECFPs with other molecular distance metrics!
    '''
    sims = np.zeros((l, l))
    # Use modern MorganGenerator to avoid deprecation warnings
    try:
        from rdkit.Chem.rdMolDescriptors import GetMorganFingerprintAsBitVect
        fps = [GetMorganFingerprintAsBitVect(x, radius=2, nBits=2048) for x in mols]
    except ImportError:
        # Fallback to older method if needed
        fps = [Chem.AllChem.GetMorganFingerprintAsBitVect(x, radius=2, nBits=2048) for x in mols]
    disable = (l <= 2000)
    for i in tqdm(range(l), disable=(l < 2000)):
        sims[i, i] = 1
        for j in range(i + 1, l):
            sims[i, j] = DataStructs.FingerprintSimilarity(fps[i], fps[j])
            sims[j, i] = sims[i, j]
    dists = 1 - sims
    return dists

    
def diversity_all(smiles=None, mols=None, dists=None, mode="HamDiv", args=None, disable=False):
    """
    Calculate various molecular diversity metrics.
    
    Args:
        smiles: List of SMILES strings
        mols: List of RDKit molecule objects
        dists: Precomputed distance matrix
        mode: Diversity metric to calculate
        args: Additional arguments
        disable: Disable progress bar
        
    Returns:
        float: Diversity value
    """
    if mode == "Richness":
        if smiles != None:
            return len(set(smiles))
        else:
            smiles = set()
            for mol in mols:
                smiles.add(Chem.MolToSmiles(mol))
            return len(smiles)
    elif mode == "FG":
        func_groups = set()
        for i in range(len(mols) if mols is not None else len(smiles)):
            smi = Chem.MolToSmiles(mols[i]) if smiles is None else smiles[i]
            grps = identify_functional_groups(smi)
            func_groups.update(grps)

        return(len(func_groups))

    elif mode == "RS":
        ring_sys = set()
        for i in range(len(mols) if mols is not None else len(smiles)):
            smi = Chem.MolToSmiles(mols[i]) if smiles is None else smiles[i]
            grps = GetRingSystems(smi)
            ring_sys.update(grps)

        return(len(ring_sys))

    elif mode == "BM":
        scaffolds = set()
        for i in range(len(mols) if mols is not None else len(smiles)):
            if mols is not None:
                scaf = MurckoScaffold.GetScaffoldForMol(mols[i])
            else:
                mol = Chem.MolFromSmiles(smiles[i])
                scaf = MurckoScaffold.GetScaffoldForMol(mol)
            scaf_smi = Chem.MolToSmiles(scaf)
            scaffolds.update([scaf_smi])

        return(len(scaffolds))


    if type(dists) is np.ndarray:
        l = len(dists)
    elif mols == None:
        l = len(smiles)
        assert l >= 2
        dists = dist_array(smiles)
    else:
        l = len(mols)
        assert l >= 2
        dists = dist_array(smiles, mols)
    
    if mode == "IntDiv":
        if l == 1:
            return 0
        return np.sum(dists) / l / (l - 1)
    elif mode == "SumDiv":
        if l == 1:
            return 0
        return np.sum(dists)/ (l - 1)
    elif mode == "Diam":
        if l == 1:
            return 0
        d_max = 0
        for i in range(l):
            for j in range(i + 1, l):
                if d_max < dists[i, j]:
                    d_max = dists[i, j]
        return d_max
    elif mode == "SumDiam":
        if l == 1:
            return 0
        sum_d_max = 0
        for i in range(l):
            d_max_i = 0
            for j in range(l):
                if j != i and d_max_i < dists[i, j]:
                    d_max_i = dists[i, j]
            sum_d_max += d_max_i
        return sum_d_max
    elif mode == "Bot":
        if l == 1:
            return 0
        d_min = 1
        for i in range(l):
            for j in range(i + 1, l):
                if d_min > dists[i, j]:
                    d_min = dists[i, j]
        return d_min
    elif mode == "SumBot":
        if l == 1:
            return 0
        sum_d_min = 0
        for i in range(l):
            d_min_i = 1
            for j in range(l):
                if j != i and d_min_i > dists[i, j]:
                    d_min_i = dists[i, j]
            sum_d_min += d_min_i
        return sum_d_min
    elif mode == "DPP":
        return np.linalg.det(1 - dists)
    elif mode.split('-')[0] == 'NCircles':
        threshold = float(mode.split('-')[1])
        circs_sum = []
        for k in tqdm(range(1), disable=disable):
            circs = np.zeros(l)
            rs = np.arange(l)
            # random.shuffle(rs)
            for i in rs:
                circs_i = 1
                for j in range(l):
                    if j != i and circs[j] == 1 and dists[i, j] <= threshold:
                        circs_i = 0
                        break
                circs[i] = circs_i
            circs_sum.append(np.sum(circs))
        return np.max(np.array(circs_sum))            
    elif mode == "HamDiv":
        total = HamDiv(dists=dists)
        return total
    
    else:
        raise Exception('Mode Undefined!')


def HamDiv(smiles=None, mols=None, dists=None, method="greedy_tsp"):
    """
    Calculate Hamiltonian diversity using traveling salesman problem solution.
    
    Args:
        smiles: List of SMILES strings
        mols: List of RDKit molecule objects
        dists: Precomputed distance matrix
        method: TSP solving method
        
    Returns:
        float: Hamiltonian diversity value
    """
    l = dists.shape[0] if dists is not None else len(mols) if mols is not None else len(smiles)
    if l == 1:
        return 0
    dists = dist_array(smiles) if dists is None else dists
    
    remove = np.zeros(l)
    for i in range(l):
        for j in range(i + 1, l):
            if dists[i, j] == 0:
                remove[i] = 1
    remove = np.argwhere(remove == 1)
    dists = np.delete(dists, remove, axis=0)
    dists = np.delete(dists, remove, axis=1)
    
    G = nx.from_numpy_array(dists)
    
    if method == "exact_dp":
        tsp, total = solve_tsp_dynamic_programming(dists)
    elif method == "christofides":
        tsp = nx.approximation.christofides(G, weight='weight')
    elif method == "greedy_tsp":
        tsp = nx.approximation.greedy_tsp(G, weight='weight')
    elif method == "simulated_annealing_tsp":
        tsp = nx.approximation.simulated_annealing_tsp(G, init_cycle="greedy", weight='weight')
    elif method == "threshold_accepting_tsp":
        tsp = nx.approximation.threshold_accepting_tsp(G, init_cycle="greedy", weight='weight')
    elif method == "local_search":
        tsp, total = solve_tsp_local_search(dists, max_processing_time=300)
    else:
        raise Exception("Undefined method")
    
    if method not in ["exact_dp", "local_search"]:
        total = 0
        for i in range(1, len(tsp)):
            total += dists[tsp[i - 1], tsp[i]]
    
    return total


def diversity_vector(smiles=None, mols=None, dists=None, ncircles_threshold=0.7, disable=False):
    """
    Calculate all molecular diversity metrics and return as a vector.
    
    Args:
        smiles: List of SMILES strings
        mols: List of RDKit molecule objects
        dists: Precomputed distance matrix
        ncircles_threshold: Threshold for NCircles metric
        disable: Disable progress bar
        
    Returns:
        dict: Dictionary containing all diversity metrics
    """
    # Count-based metrics (don't need distance matrix)
    richness = diversity_all(smiles=smiles, mols=mols, mode="Richness")
    fg_diversity = diversity_all(smiles=smiles, mols=mols, mode="FG")
    rs_diversity = diversity_all(smiles=smiles, mols=mols, mode="RS")
    bm_diversity = diversity_all(smiles=smiles, mols=mols, mode="BM")
    
    # Calculate distance matrix once for all distance-based metrics
    if dists is None:
        if smiles is not None:
            l = len(smiles)
            if l >= 2:
                dists = dist_array(smiles)
        elif mols is not None:
            l = len(mols)
            if l >= 2:
                dists = dist_array(mols=mols)
        else:
            raise ValueError("Either smiles or mols must be provided")
    
    # Distance-based metrics
    if dists is not None and len(dists) >= 2:
        int_div = diversity_all(dists=dists, mode="IntDiv", disable=disable)
        sum_div = diversity_all(dists=dists, mode="SumDiv", disable=disable)
        diam = diversity_all(dists=dists, mode="Diam", disable=disable)
        sum_diam = diversity_all(dists=dists, mode="SumDiam", disable=disable)
        bot = diversity_all(dists=dists, mode="Bot", disable=disable)
        sum_bot = diversity_all(dists=dists, mode="SumBot", disable=disable)
        dpp = diversity_all(dists=dists, mode="DPP", disable=disable)
        ncircles = diversity_all(dists=dists, mode=f"NCircles-{ncircles_threshold}", disable=disable)
        ham_div = diversity_all(dists=dists, mode="HamDiv", disable=disable)
    else:
        # If only one molecule or no distance matrix, set distance-based metrics to 0
        int_div = sum_div = diam = sum_diam = bot = sum_bot = dpp = ncircles = ham_div = 0.0
    
    return {
        "richness": richness,
        "functional_groups": fg_diversity,
        "ring_systems": rs_diversity,
        "bemis_murcko": bm_diversity,
        "internal_diversity": int_div,
        "sum_diversity": sum_div,
        "diameter": diam,
        "sum_diameter": sum_diam,
        "bottleneck": bot,
        "sum_bottleneck": sum_bot,
        "determinantal_point_process": dpp,
        "n_circles": ncircles,
        "hamiltonian_diversity": ham_div
    }


def diversity_vector_array(smiles=None, mols=None, dists=None, ncircles_threshold=0.7, disable=False):
    """
    Calculate all molecular diversity metrics and return as a numpy array.
    
    Args:
        smiles: List of SMILES strings
        mols: List of RDKit molecule objects
        dists: Precomputed distance matrix
        ncircles_threshold: Threshold for NCircles metric
        disable: Disable progress bar
        
    Returns:
        tuple: (metric_names, values_array) where metric_names is a list of metric names
               and values_array is a numpy array of corresponding values
    """
    diversity_dict = diversity_vector(smiles=smiles, mols=mols, dists=dists, 
                                    ncircles_threshold=ncircles_threshold, disable=disable)
    
    metric_names = list(diversity_dict.keys())
    values_array = np.array(list(diversity_dict.values()))
    
    return metric_names, values_array