import networkx as nx
import numpy as np
from math import sqrt
from tqdm import tqdm

NUMBER_OF_ORBITS = 73

class GDDA:
    """
    This class computes the GDD-agreement.
    """
    
    def agreement(self, G, H, method="arith", verbose=False):
        """
        This method computes the GDD-agreement between graphs G and H.

        Arguments
        ----------
        G, H : networkx.Graph
            Networkx graphs.

        method : "arith" or "geo" or None
            Method to calculate the average.
            If None, return a vector before calculating the GDD-agreement.
        
        verbose : bool
            Show progress bar if True.
        """
        assert (method in ["arith", "geo"]) or (method is None)
        
        # Computing GDDs
        if verbose: print("Computing the GDD of G.")
        gdd_G = self.distribution(G.to_undirected(), verbose=verbose)
        if verbose: print("Computing the GDD of H.")
        gdd_H = self.distribution(H.to_undirected(), verbose=verbose)
        
        # Computing the vector before the GDD-agreement.
        gdda_vector = GDD_agreement_vector(gdd_G, gdd_H)
        
        if method is None:
            return gdda_vector
        elif method == "arith":
            return gdda_vector.mean()
        else:
            return np.exp(np.log(gdda_vector).mean())
    
    
    def distribution(self, G, verbose=False, orbit_mode=False):
        """
        This method computes th GDD of graph G.
        
        Arguments
        ----------
        G : networkx.Graph
            A networkx graph.
        
        verbose : bool
            Show progress bar if True.
        """
        adjacency_matrix = nx.to_scipy_sparse_matrix(G, format="lil")
        return graphlet_degree_distribution(adjacency_matrix, verbose=verbose, orbit_mode=orbit_mode)


def graphlet_degree_distribution(adjacency_matrix, verbose=False, orbit_mode=False):
    """
    This function computes the GDD of a graph with adjacency_matrix.
    
    Arguments
    ----------
    adjacency_matrix : scipy.lil_matrix
        An adjacency matrix.
    
    verbose : bool
        Show progress bar if True.
    
    orbit_mode : bool
        Return orbit features.
    """
    G = nx.from_scipy_sparse_matrix(adjacency_matrix)
    
    V = G.number_of_nodes()
    N = [ set(G.neighbors(v)) for v in G.nodes() ]
    orbits = np.zeros((V, NUMBER_OF_ORBITS), dtype=np.int64)
    
    # Orbit 0 is a degree distribution
    orbits[:, 0] = [ d for v, d in G.degree() ]
    
    
    def comb(n):
        return n * (n - 1) // 2
    
    if verbose:
        edges = tqdm(G.edges(), total=G.number_of_edges())
    else:
        edges = G.edges()
    
    # Orbits 1-72
    for u, v in edges:
        N_u = N[u] - {v}
        N_v = N[v] - {u}
        U_e = N[u] | N[v]
        
        Star_u = N_u - N_v
        Star_v = N_v - N_u
        Tri_e  = N_u & N_v
        
        n_u   = len(Star_u)
        n_v   = len(Star_v)
        n_e   = len(Tri_e)
        n_uu  = comb(n_u)
        n_vv  = comb(n_v)
        n_ee  = comb(n_e)
        n_uv  = n_u * n_v
        n_ue  = n_u * n_e
        n_ve  = n_v * n_e
        n_uuu = n_uu * (n_u - 2) // 3
        n_vvv = n_vv * (n_v - 2) // 3
        n_eee = n_ee * (n_e - 2) // 3
        n_uuv = n_uu * n_v
        n_uvv = n_vv * n_u
        n_uue = n_uu * n_e
        n_vve = n_vv * n_e
        n_uee = n_u  * n_ee
        n_vee = n_v  * n_ee
        n_uve = n_uv * n_e
        
        # 4-node graphlets
        n_3_1_u = 0
        n_3_1_v = 0
        #n_3_2   = 0
        #n_4_1_u = 0
        #n_4_1_v = 0
        n_5_1   = 0
        n_6_1_u = 0
        n_6_1_v = 0
        #n_6_2_u = 0
        #n_6_2_v = 0
        n_6_3   = 0
        n_7_1_u = 0
        n_7_1_v = 0
        #n_7_2   = 0
        n_8_1   = 0
        
        # 5-node graphlets
        n_9_1_u  = 0
        n_9_1_v  = 0
        n_9_2_u  = 0
        n_9_2_v  = 0
        n_10_1_u = 0
        n_10_1_v = 0
        #n_10_2_u = 0
        #n_10_2_v = 0
        n_10_3_u = 0
        n_10_3_v = 0
        #n_11_1_u = 0
        #n_11_1_v = 0
        n_12_1_u = 0
        n_12_1_v = 0
        n_12_2_u = 0
        n_12_2_v = 0
        n_12_3   = 0
        n_13_1_u = 0
        n_13_1_v = 0
        n_13_2_u = 0
        n_13_2_v = 0
        n_13_3_u = 0
        n_13_3_v = 0
        n_13_4   = 0
        n_14_1_u = 0
        n_14_1_v = 0
        #n_14_2_u = 0
        #n_14_2_v = 0
        n_14_3   = 0
        n_15_1   = 0
        n_16_1_u = 0
        n_16_1_v = 0
        n_16_2_u = 0
        n_16_2_v = 0
        n_16_3_u = 0
        n_16_3_v = 0
        n_17_1_u = 0
        n_17_1_v = 0
        n_17_2_u = 0
        n_17_2_v = 0
        n_17_3_u = 0
        n_17_3_v = 0
        n_17_4_u = 0
        n_17_4_v = 0
        n_18_1   = 0
        n_18_2_u = 0
        n_18_2_v = 0
        n_19_1_u = 0
        n_19_1_v = 0
        n_19_2_u = 0
        n_19_2_v = 0
        n_19_3_u = 0
        n_19_3_v = 0
        n_19_4   = 0
        n_20_1_u = 0
        n_20_1_v = 0
        n_21_1   = 0
        n_21_2_u = 0
        n_21_2_v = 0
        n_21_3_u = 0
        n_21_3_v = 0
        n_21_4   = 0
        n_22_1_u = 0
        n_22_1_v = 0
        #n_22_2   = 0
        n_23_1_u = 0
        n_23_1_v = 0
        n_23_2_u = 0
        n_23_2_v = 0
        n_23_3   = 0
        n_24_1_u = 0
        n_24_1_v = 0
        n_24_2_u = 0
        n_24_2_v = 0
        n_24_3   = 0
        n_24_4_u = 0
        n_24_4_v = 0
        n_25_1_u = 0
        n_25_1_v = 0
        n_25_2_u = 0
        n_25_2_v = 0
        n_25_3   = 0
        n_26_1_u = 0
        n_26_1_v = 0
        n_26_2_u = 0
        n_26_2_v = 0
        n_26_3   = 0
        n_26_4   = 0
        n_27_1   = 0
        n_27_2_u = 0
        n_27_2_v = 0
        n_28_1_u = 0
        n_28_1_v = 0
        n_28_2   = 0
        n_29_1   = 0
        
        for w in Star_u:
            N_w  = N[w] - {u}
            NwUe = N_w - U_e
            NwSu = N_w & Star_u
            NwSv = N_w & Star_v
            NwTe = N_w & Tri_e
            UeNw = U_e | N_w
            SuNw = Star_u - N_w - {w}
            SvNw = Star_v - N_w
            TeNw = Tri_e  - N_w
            
            n_NwUe = len(NwUe)
            n_NwSu = len(NwSu)
            n_NwSv = len(NwSv)
            n_NwTe = len(NwTe)
            n_SuNw = len(SuNw)
            n_SvNw = len(SvNw)
            n_TeNw = len(TeNw)
            
            # 4-node graphlets
            n_3_1_u += n_NwUe
            n_6_1_u += n_NwSu
            n_5_1   += n_NwSv
            n_7_1_u += n_NwTe
            
            # 5-node graphlets
            for r in N_w:
                N_r = N[r] - {w}
                if r in Star_u:
                    n_23_1_u += len(N_r & NwSu)
                    n_25_1_u += len(N_r & NwSv)
                    n_26_1_u += len(N_r & NwTe)
                elif r in Star_v:
                    n_27_1   += len(N_r & NwTe)
                elif r in Tri_e:
                    n_28_1_u += len(N_r & NwTe)
                else:
                    n_9_1_u  += len(N_r - UeNw)
                    n_13_1_u += len(N_r & NwUe)
                    n_15_1   += len(N_r & SvNw)
                    n_16_1_u += len(N_r & SuNw)
                    n_19_1_u += len(N_r & NwSu)
                    n_21_1   += len(N_r & NwSv)
                    n_21_3_u += len(N_r & TeNw)
                    n_24_2_u += len(N_r & NwTe)
            
            n_9_2_u  += n_NwUe * n_SvNw
            n_10_1_u += comb(n_NwUe)
            n_10_3_u += n_NwUe * n_SuNw
            n_12_1_u += n_NwUe * n_NwSu
            n_13_2_u += n_NwSu * n_SvNw
            n_13_3_u += n_NwUe * n_TeNw
            n_14_1_u += n_NwSu * n_SuNw
            n_16_2_u += n_NwSv * n_SvNw
            n_16_3_u += n_NwUe * n_NwSv
            n_17_1_u += comb(n_NwSu)
            n_17_2_v += n_NwTe * n_SuNw
            n_18_2_u += n_NwSu * n_TeNw
            n_19_2_u += n_NwTe * n_SvNw
            n_19_3_u += n_NwUe * n_NwTe
            n_20_1_u += comb(n_NwSv)
            n_21_2_u += n_NwSu * n_NwSv
            n_21_4   += n_NwSv * n_TeNw
            n_24_1_u += n_NwSu * n_NwTe
            n_24_4_u += n_NwTe * n_TeNw
            n_25_2_u += n_NwSv * n_NwTe
            n_27_2_u += comb(n_NwTe)
        
        
        for w in Star_v:
            N_w  = N[w] - {v}
            NwUe = N_w - U_e
            NwSu = N_w & Star_u
            NwSv = N_w & Star_v
            NwTe = N_w & Tri_e
            UeNw = U_e | N_w
            SuNw = Star_u - N_w
            SvNw = Star_v - N_w - {w}
            TeNw = Tri_e  - N_w
            
            n_NwUe = len(NwUe)
            n_NwSu = len(NwSu)
            n_NwSv = len(NwSv)
            n_NwTe = len(NwTe)
            n_SuNw = len(SuNw)
            n_SvNw = len(SvNw)
            n_TeNw = len(TeNw)
            
            # 4-node graphlets
            n_3_1_v += n_NwUe
            n_6_1_v += n_NwSv
            n_7_1_v += n_NwTe
            
            # 5-node graphlets
            for r in N_w:
                N_r = N[r] - {w}
                if r in Star_v:
                    n_23_1_v += len(N_r & NwSv)
                    n_25_1_v += len(N_r & NwSu)
                    n_26_1_v += len(N_r & NwTe)
                elif r in Star_u:
                    continue
                elif r in Tri_e:
                    n_28_1_v += len(N_r & NwTe)
                else:
                    n_9_1_v  += len(N_r - UeNw)
                    n_13_1_v += len(N_r & NwUe)
                    n_16_1_v += len(N_r & SvNw)
                    n_19_1_v += len(N_r & NwSv)
                    n_21_3_v += len(N_r & TeNw)
                    n_24_2_v += len(N_r & NwTe)
            
            n_9_2_v  += n_NwUe * n_SuNw
            n_10_1_v += comb(n_NwUe)
            n_10_3_v += n_NwUe * n_SvNw
            n_12_1_v += n_NwUe * n_NwSv
            n_13_2_v += n_NwSv * n_SuNw
            n_13_3_v += n_NwUe * n_TeNw
            n_14_1_v += n_NwSv * n_SvNw
            n_16_2_v += n_NwSu * n_SuNw
            n_16_3_v += n_NwUe * n_NwSu
            n_17_1_v += comb(n_NwSv)
            n_17_2_u += n_NwTe * n_SvNw
            n_18_2_v += n_NwSv * n_TeNw
            n_19_2_v += n_NwTe * n_SuNw
            n_19_3_v += n_NwUe * n_NwTe
            n_20_1_v += comb(n_NwSu)
            n_21_2_v += n_NwSv * n_NwSu
            n_24_1_v += n_NwSv * n_NwTe
            n_24_4_v += n_NwTe * n_TeNw
            n_25_2_v += n_NwSu * n_NwTe
            n_27_2_v += comb(n_NwTe)
        
        
        for w in Tri_e:
            N_w  = N[w] - {u, v}
            NwUe = N_w - U_e
            NwSu = N_w & Star_u
            NwSv = N_w & Star_v
            NwTe = N_w & Tri_e
            UeNw = U_e | N_w
            SuNw = Star_u - N_w
            SvNw = Star_v - N_w
            TeNw = Tri_e  - N_w - {w}
            
            n_NwUe = len(NwUe)
            n_NwSu = len(NwSu)
            n_NwSv = len(NwSv)
            n_NwTe = len(NwTe)
            n_SuNw = len(SuNw)
            n_SvNw = len(SvNw)
            n_TeNw = len(TeNw)
            
            # 4-node graphlets
            n_6_3 += n_NwUe
            n_8_1 += n_NwTe
            
            # 5-node graphlets
            for r in N_w:
                N_r = N[r] - {w}
                if r in Star_u:
                    continue
                elif r in Star_v:
                    continue
                elif r in Tri_e:
                    n_29_1 += len(N_r & NwTe)
                else:
                    n_13_4 += len(N_r - UeNw)
                    n_18_1 += len(N_r & NwUe)
                    n_25_3 += len(N_r & TeNw)
                    n_26_3 += len(N_r & NwTe)
            
            n_12_2_u += n_NwUe * n_SvNw
            n_12_2_v += n_NwUe * n_SuNw
            n_14_3   += comb(n_NwUe)
            n_17_3_u += n_NwUe * n_NwSv
            n_17_3_v += n_NwUe * n_NwSu
            n_19_4   += n_NwUe * n_TeNw
            n_22_1_u += comb(n_NwSv)
            n_22_1_v += comb(n_NwSu)
            n_23_2_u += n_NwTe * n_SvNw
            n_23_2_v += n_NwTe * n_SuNw
            n_23_3   += n_NwUe * n_NwTe
            n_24_3   += n_NwSu * n_NwSv
            n_26_2_u += n_NwSv * n_NwTe
            n_26_2_v += n_NwSu * n_NwTe
            n_26_4   += n_NwTe * n_TeNw
            n_28_2   += comb(n_NwTe)
            
        
        # 4-node graphlets
        n_6_1_u //= 2
        n_6_1_v //= 2
        n_8_1   //= 2
        
        n_3_2   = n_uv - n_5_1
        n_4_1_u = n_uu - n_6_1_u
        n_4_1_v = n_vv - n_6_1_v
        n_6_2_u = n_ue - n_7_1_u
        n_6_2_v = n_ve - n_7_1_v
        n_7_2   = n_ee - n_8_1
        
        # 5-node graphlets
        n_13_1_u //= 2
        n_13_1_v //= 2
        n_16_1_u //= 2
        n_16_1_v //= 2
        n_18_1   //= 2
        n_19_1_u //= 2
        n_19_1_v //= 2
        n_23_1_u //= 6
        n_23_1_v //= 6
        n_25_1_u //= 2
        n_25_1_v //= 2
        n_25_3   //= 2
        n_26_1_u //= 2
        n_26_1_v //= 2
        n_26_3   //= 2
        n_28_1_u //= 2
        n_28_1_v //= 2
        n_29_1   //= 6
        
        n_9_2_u  -= n_15_1
        n_9_2_v  -= n_15_1
        n_10_1_u -= n_13_1_u
        n_10_1_v -= n_13_1_v
        n_10_3_u -= 2 * n_16_1_u
        n_10_3_v -= 2 * n_16_1_v
        n_12_1_u -= 2 * n_19_1_u
        n_12_1_v -= 2 * n_19_1_v
        n_12_2_u -= n_21_3_v
        n_12_2_v -= n_21_3_u
        n_13_3_u -= n_21_3_u
        n_13_3_v -= n_21_3_v
        n_14_3   -= n_18_1
        n_16_3_u -= n_21_1
        n_16_3_v -= n_21_1
        n_17_1_u -= 3 * n_23_1_u
        n_17_1_v -= 3 * n_23_1_v
        n_17_3_u -= n_24_2_v
        n_17_3_v -= n_24_2_u
        n_19_3_u -= n_24_2_u
        n_19_3_v -= n_24_2_v
        n_19_4   -= 2 * n_25_3
        n_20_1_u -= n_25_1_v
        n_20_1_v -= n_25_1_u
        n_21_2_u -= 2 * n_25_1_u
        n_21_2_v -= 2 * n_25_1_v
        n_22_1_u -= n_26_1_v
        n_22_1_v -= n_26_1_u
        n_23_3   -= 2 * n_26_3
        n_24_1_u -= 2 * n_26_1_u
        n_24_1_v -= 2 * n_26_1_v
        n_24_3   -= n_27_1
        n_25_2_u -= n_27_1
        n_25_2_v -= n_27_1
        n_26_2_u -= 2 * n_28_1_v
        n_26_2_v -= 2 * n_28_1_u
        n_27_2_u -= n_28_1_u
        n_27_2_v -= n_28_1_v
        n_28_2    = (n_28_2 - 3 * n_29_1) // 2
        
        n_13_2_u  = (n_13_2_u - n_21_2_u) // 2
        n_13_2_v  = (n_13_2_v - n_21_2_v) // 2
        n_14_1_u  = (n_14_1_u - 2 * n_17_1_u) // 2
        n_14_1_v  = (n_14_1_v - 2 * n_17_1_v) // 2
        n_16_2_u -= n_21_2_v
        n_16_2_v -= n_21_2_u
        n_17_2_u -= 2 * n_22_1_u
        n_17_2_v -= 2 * n_22_1_v
        n_18_2_u  = (n_18_2_u - n_24_1_u) // 2
        n_18_2_v  = (n_18_2_v - n_24_1_v) // 2
        n_19_2_u -= n_24_3
        n_19_2_v -= n_24_3
        n_21_4   -= n_25_2_v
        n_23_2_u  = (n_23_2_u - n_26_2_u) // 2
        n_23_2_v  = (n_23_2_v - n_26_2_v) // 2
        n_24_4_u -= n_26_2_v
        n_24_4_v -= n_26_2_u
        n_26_4    = (n_26_4 - 2 * n_28_2) // 2
        
        n_11_1_u = n_uuu - (n_14_1_u + n_17_1_u + n_23_1_u)
        n_11_1_v = n_vvv - (n_14_1_v + n_17_1_v + n_23_1_v)
        n_10_2_u = n_uuv - (n_13_2_u + n_16_2_v + n_21_2_u + n_20_1_v + n_25_1_u)
        n_10_2_v = n_uvv - (n_13_2_v + n_16_2_u + n_21_2_v + n_20_1_u + n_25_1_v)
        n_14_2_u = n_uue - (n_17_2_v + n_18_2_u + n_22_1_v + n_24_1_u + n_26_1_u)
        n_14_2_v = n_vve - (n_17_2_u + n_18_2_v + n_22_1_u + n_24_1_v + n_26_1_v)
        n_12_3   = n_uve - (n_19_2_u + n_19_2_v + n_21_4   + n_24_3   + n_25_2_u + n_25_2_v + n_27_1)
        n_17_4_u = n_uee - (n_23_2_v + n_24_4_u + n_26_2_v + n_27_2_u + n_28_1_u)
        n_17_4_v = n_vee - (n_23_2_u + n_24_4_v + n_26_2_u + n_27_2_v + n_28_1_v)
        n_22_2   = n_eee - (n_26_4   + n_28_2   + n_29_1)
        
        # Counting degrees
        orbits[u, 1] += n_v
        orbits[v, 1] += n_u
        orbits[u, 2] += n_u
        orbits[v, 2] += n_v
        orbits[u, 3] += n_e
        orbits[v, 3] += n_e
        
        orbits[u,  4] += n_3_1_v
        orbits[v,  4] += n_3_1_u
        orbits[u,  5] += n_3_1_u + n_3_2
        orbits[v,  5] += n_3_1_v + n_3_2
        orbits[u,  6] += n_4_1_v
        orbits[v,  6] += n_4_1_u
        orbits[u,  7] += n_4_1_u
        orbits[v,  7] += n_4_1_v
        orbits[u,  8] += n_5_1
        orbits[v,  8] += n_5_1
        orbits[u,  9] += n_6_1_v
        orbits[v,  9] += n_6_1_u
        orbits[u, 10] += n_6_2_v + n_6_3
        orbits[v, 10] += n_6_2_u + n_6_3
        orbits[u, 11] += n_6_1_u + n_6_2_u
        orbits[v, 11] += n_6_1_v + n_6_2_v
        orbits[u, 12] += n_7_1_v
        orbits[v, 12] += n_7_1_u
        orbits[u, 13] += n_7_1_u + n_7_2
        orbits[v, 13] += n_7_1_v + n_7_2
        orbits[u, 14] += n_8_1
        orbits[v, 14] += n_8_1
        
        orbits[u, 15] += n_9_1_v
        orbits[v, 15] += n_9_1_u
        orbits[u, 16] += n_9_1_u + n_9_2_v
        orbits[v, 16] += n_9_1_v + n_9_2_u
        orbits[u, 17] += n_9_2_u
        orbits[v, 17] += n_9_2_v
        orbits[u, 18] += n_10_1_v
        orbits[v, 18] += n_10_1_u
        orbits[u, 19] += n_10_3_v
        orbits[v, 19] += n_10_3_u
        orbits[u, 20] += n_10_1_u + n_10_2_v
        orbits[v, 20] += n_10_1_v + n_10_2_u
        orbits[u, 21] += n_10_2_u + n_10_3_u
        orbits[v, 21] += n_10_2_v + n_10_3_v
        orbits[u, 22] += n_11_1_v
        orbits[v, 22] += n_11_1_u
        orbits[u, 23] += n_11_1_u
        orbits[v, 23] += n_11_1_v
        orbits[u, 24] += n_12_1_v
        orbits[v, 24] += n_12_1_u
        orbits[u, 25] += n_12_2_u
        orbits[v, 25] += n_12_2_v
        orbits[u, 26] += n_12_1_u + n_12_2_v + n_12_3
        orbits[v, 26] += n_12_1_v + n_12_2_u + n_12_3
        orbits[u, 27] += n_13_1_v
        orbits[v, 27] += n_13_1_u
        orbits[u, 28] += n_13_1_u + n_13_2_v
        orbits[v, 28] += n_13_1_v + n_13_2_u
        orbits[u, 29] += n_13_3_v + n_13_4
        orbits[v, 29] += n_13_3_u + n_13_4
        orbits[u, 30] += n_13_2_u + n_13_3_u
        orbits[v, 30] += n_13_2_v + n_13_3_v
        orbits[u, 31] += n_14_1_v
        orbits[v, 31] += n_14_1_u
        orbits[u, 32] += n_14_2_v + n_14_3
        orbits[v, 32] += n_14_2_u + n_14_3
        orbits[u, 33] += n_14_1_u + n_14_2_u
        orbits[v, 33] += n_14_1_v + n_14_2_v
        orbits[u, 34] += n_15_1
        orbits[v, 34] += n_15_1
        orbits[u, 35] += n_16_1_v
        orbits[v, 35] += n_16_1_u
        orbits[u, 36] += n_16_3_v
        orbits[v, 36] += n_16_3_u
        orbits[u, 37] += n_16_2_u + n_16_3_u
        orbits[v, 37] += n_16_2_v + n_16_3_v
        orbits[u, 38] += n_16_1_u + n_16_2_v
        orbits[v, 38] += n_16_1_v + n_16_2_u
        orbits[u, 39] += n_17_1_v
        orbits[v, 39] += n_17_1_u
        orbits[u, 40] += n_17_2_u + n_17_3_u
        orbits[v, 40] += n_17_2_v + n_17_3_v
        orbits[u, 41] += n_17_3_v + n_17_4_v
        orbits[v, 41] += n_17_3_u + n_17_4_u
        orbits[u, 42] += n_17_1_u + n_17_2_v + n_17_4_u
        orbits[v, 42] += n_17_1_v + n_17_2_u + n_17_4_v
        orbits[u, 43] += n_18_1   + n_18_2_v
        orbits[v, 43] += n_18_1   + n_18_2_u
        orbits[u, 44] += n_18_2_u
        orbits[v, 44] += n_18_2_v
        orbits[u, 45] += n_19_1_v
        orbits[v, 45] += n_19_1_u
        orbits[u, 46] += n_19_3_v
        orbits[v, 46] += n_19_3_u
        orbits[u, 47] += n_19_1_u + n_19_2_v
        orbits[v, 47] += n_19_1_v + n_19_2_u
        orbits[u, 48] += n_19_2_u + n_19_3_u + n_19_4
        orbits[v, 48] += n_19_2_v + n_19_3_v + n_19_4
        orbits[u, 49] += n_20_1_u
        orbits[v, 49] += n_20_1_v
        orbits[u, 50] += n_20_1_v
        orbits[v, 50] += n_20_1_u
        orbits[u, 51] += n_21_1   + n_21_2_v
        orbits[v, 51] += n_21_1   + n_21_2_u
        orbits[u, 52] += n_21_3_v
        orbits[v, 52] += n_21_3_u
        orbits[u, 53] += n_21_2_u + n_21_3_u + n_21_4
        orbits[v, 53] += n_21_2_v + n_21_3_v + n_21_4
        orbits[u, 54] += n_22_1_u
        orbits[v, 54] += n_22_1_v
        orbits[u, 55] += n_22_1_v + n_22_2
        orbits[v, 55] += n_22_1_u + n_22_2
        orbits[u, 56] += n_23_1_v
        orbits[v, 56] += n_23_1_u
        orbits[u, 57] += n_23_2_u + n_23_3
        orbits[v, 57] += n_23_2_v + n_23_3
        orbits[u, 58] += n_23_1_u + n_23_2_v
        orbits[v, 58] += n_23_1_v + n_23_2_u
        orbits[u, 59] += n_24_1_v + n_24_2_v
        orbits[v, 59] += n_24_1_u + n_24_2_u
        orbits[u, 60] += n_24_2_u + n_24_3   + n_24_4_v
        orbits[v, 60] += n_24_2_v + n_24_3   + n_24_4_u
        orbits[u, 61] += n_24_1_u + n_24_4_u
        orbits[v, 61] += n_24_1_v + n_24_4_v
        orbits[u, 62] += n_25_1_v
        orbits[v, 62] += n_25_1_u
        orbits[u, 63] += n_25_1_u + n_25_2_v
        orbits[v, 63] += n_25_1_v + n_25_2_u
        orbits[u, 64] += n_25_2_u + n_25_3
        orbits[v, 64] += n_25_2_v + n_25_3
        orbits[u, 65] += n_26_1_v
        orbits[v, 65] += n_26_1_u
        orbits[u, 66] += n_26_2_u + n_26_3
        orbits[v, 66] += n_26_2_v + n_26_3
        orbits[u, 67] += n_26_1_u + n_26_2_v + n_26_4
        orbits[v, 67] += n_26_1_v + n_26_2_u + n_26_4
        orbits[u, 68] += n_27_1   + n_27_2_v
        orbits[v, 68] += n_27_1   + n_27_2_u
        orbits[u, 69] += n_27_2_u
        orbits[v, 69] += n_27_2_v
        orbits[u, 70] += n_28_1_v
        orbits[v, 70] += n_28_1_u
        orbits[u, 71] += n_28_1_u + n_28_2
        orbits[v, 71] += n_28_1_v + n_28_2
        orbits[u, 72] += n_29_1
        orbits[v, 72] += n_29_1
    
    if orbit_mode:
        return orbits
    
    # Converting to GDD
    def to_DD(x):
        return np.bincount(x)
    
    gdd = [ None for _ in range(NUMBER_OF_ORBITS) ]
    
    # 2-node graphlet
    gdd[0]  = to_DD(orbits[:,  0])
    
    # 3-node graphlets
    gdd[1]  = to_DD(orbits[:,  1])
    gdd[2]  = to_DD(orbits[:,  2] // 2)
    gdd[3]  = to_DD(orbits[:,  3] // 2)
    
    # 4-node graphlets
    gdd[4]  = to_DD(orbits[:,  4])
    gdd[5]  = to_DD(orbits[:,  5] // 2)
    gdd[6]  = to_DD(orbits[:,  6])
    gdd[7]  = to_DD(orbits[:,  7] // 3)
    gdd[8]  = to_DD(orbits[:,  8] // 2)
    gdd[9]  = to_DD(orbits[:,  9])
    gdd[10] = to_DD(orbits[:, 10] // 2)
    gdd[11] = to_DD(orbits[:, 11] // 3)
    gdd[12] = to_DD(orbits[:, 12] // 2)
    gdd[13] = to_DD(orbits[:, 13] // 3)
    gdd[14] = to_DD(orbits[:, 14] // 3)
    
    # 5-node graphlets
    gdd[15] = to_DD(orbits[:, 15])
    gdd[16] = to_DD(orbits[:, 16] // 2)
    gdd[17] = to_DD(orbits[:, 17] // 2)
    gdd[18] = to_DD(orbits[:, 18])
    gdd[19] = to_DD(orbits[:, 19])
    gdd[20] = to_DD(orbits[:, 20] // 2)
    gdd[21] = to_DD(orbits[:, 21] // 3)
    gdd[22] = to_DD(orbits[:, 22])
    gdd[23] = to_DD(orbits[:, 23] // 4)
    gdd[24] = to_DD(orbits[:, 24])
    gdd[25] = to_DD(orbits[:, 25] // 2)
    gdd[26] = to_DD(orbits[:, 26] // 3)
    gdd[27] = to_DD(orbits[:, 27])
    gdd[28] = to_DD(orbits[:, 28] // 2)
    gdd[29] = to_DD(orbits[:, 29] // 2)
    gdd[30] = to_DD(orbits[:, 30] // 3)
    gdd[31] = to_DD(orbits[:, 31])
    gdd[32] = to_DD(orbits[:, 32] // 2)
    gdd[33] = to_DD(orbits[:, 33] // 4)
    gdd[34] = to_DD(orbits[:, 34] // 2)
    gdd[35] = to_DD(orbits[:, 35])
    gdd[36] = to_DD(orbits[:, 36] // 2)
    gdd[37] = to_DD(orbits[:, 37] // 2)
    gdd[38] = to_DD(orbits[:, 38] // 3)
    gdd[39] = to_DD(orbits[:, 39])
    gdd[40] = to_DD(orbits[:, 40] // 2)
    gdd[41] = to_DD(orbits[:, 41] // 3)
    gdd[42] = to_DD(orbits[:, 42] // 4)
    gdd[43] = to_DD(orbits[:, 43] // 2)
    gdd[44] = to_DD(orbits[:, 44] // 4)
    gdd[45] = to_DD(orbits[:, 45])
    gdd[46] = to_DD(orbits[:, 46] // 2)
    gdd[47] = to_DD(orbits[:, 47] // 3)
    gdd[48] = to_DD(orbits[:, 48] // 3)
    gdd[49] = to_DD(orbits[:, 49] // 2)
    gdd[50] = to_DD(orbits[:, 50] // 3)
    gdd[51] = to_DD(orbits[:, 51] // 2)
    gdd[52] = to_DD(orbits[:, 52] // 2)
    gdd[53] = to_DD(orbits[:, 53] // 3)
    gdd[54] = to_DD(orbits[:, 54] // 2)
    gdd[55] = to_DD(orbits[:, 55] // 4)
    gdd[56] = to_DD(orbits[:, 56])
    gdd[57] = to_DD(orbits[:, 57] // 3)
    gdd[58] = to_DD(orbits[:, 58] // 4)
    gdd[59] = to_DD(orbits[:, 59] // 2)
    gdd[60] = to_DD(orbits[:, 60] // 3)
    gdd[61] = to_DD(orbits[:, 61] // 4)
    gdd[62] = to_DD(orbits[:, 62] // 2)
    gdd[63] = to_DD(orbits[:, 63] // 3)
    gdd[64] = to_DD(orbits[:, 64] // 3)
    gdd[65] = to_DD(orbits[:, 65] // 2)
    gdd[66] = to_DD(orbits[:, 66] // 3)
    gdd[67] = to_DD(orbits[:, 67] // 4)
    gdd[68] = to_DD(orbits[:, 68] // 3)
    gdd[69] = to_DD(orbits[:, 69] // 4)
    gdd[70] = to_DD(orbits[:, 70] // 3)
    gdd[71] = to_DD(orbits[:, 71] // 4)
    gdd[72] = to_DD(orbits[:, 72] // 4)
    
    return gdd


def GDD_agreement_vector(gdd_G, gdd_H):
    """
    This function computes the GDD-agreement between gdd_G and gdd_H.
    
    Arguments
    ----------
    gdd_G, gdd_H : list of numpy.ndarray
        GDDs.
    """
    # Normalization
    norm_gdd_G = normalized_distribution(gdd_G)
    norm_gdd_H = normalized_distribution(gdd_H)
    sq2 = sqrt(2)
    
    gdda_vector = np.array([
        1 - padded_norm(norm_gdd_G[j], norm_gdd_H[j]) / sq2 for j in range(NUMBER_OF_ORBITS)
    ])
    
    return gdda_vector


def normalized_distribution(gdd):
    # Scaling
    S = [ np.array([ d / k if k != 0 else 0 for k, d in enumerate(dd_j) ])  for dd_j in gdd ]
    
    # Total
    T = [ s_j.sum() for s_j in S ]
    
    # Normalization
    N = [ s_j / T[j] if T[j] != 0 else np.array([0.0]) for j, s_j in enumerate(S) ]
    
    return N


def padded_norm(x, y):
    nx = len(x)
    ny = len(y)
    
    if nx > ny:
        y_pad = np.pad(y, (0, nx - ny))
        return np.linalg.norm(x - y_pad)
    elif nx < ny:
        x_pad = np.pad(x, (0, ny - nx))
        return np.linalg.norm(x_pad - y)
    else:
        return np.linalg.norm(x - y)
