"""Additional term implementation for molecular memory in generation."""

import pandas as pd
import numpy as np


def memory_update(smiles, scores, memory=None, sigma1=0.1, threshold=0.5, calc_fingerprints=None, calc_dists=None):
    """
    Update molecular memory with diversity enhancement.
    
    Args:
        smiles: List of SMILES strings
        scores: List of scores for molecules
        memory: DataFrame with molecular memory
        sigma1: Diversity enhancement factor
        threshold: Score threshold for memory update
        calc_fingerprints: Function to calculate fingerprints
        calc_dists: Function to calculate distances
        
    Returns:
        tuple: Updated smiles and scores
    """
    if memory is None:
        memory = pd.DataFrame(columns=["smiles", "scores", "fps"])
    
    for i in range(len(smiles)):
        if calc_fingerprints is not None:
            fp = calc_fingerprints([smiles[i]])
            
            if calc_dists is not None and not memory.empty:
                dists = calc_dists(fp, memory["fps"])
                # Adding the diversity enhancement term
                scores[i] += sigma1 * np.min(dists)

            # Update the memory
            if scores[i] > threshold:
                new_data = pd.DataFrame({"smiles": [smiles[i]], "scores": [scores[i]], "fps": [fp]})
                memory = pd.concat([memory, new_data], ignore_index=True, sort=False)

    return smiles, scores, memory