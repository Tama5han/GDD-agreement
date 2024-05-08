"""
Reference
----------
1. N. Przulj. (2007) Biological network comparison using graphlet degree distribution.
   Bioinformatics. 23(2), e177-e183. <https://doi.org/10.1093/bioinformatics/btl301>
"""

import networkx as nx
import numpy as np
from math import sqrt
from tqdm import tqdm

# for 2~5-node graphlets
NUMBER_OF_ORBITS = 73



def agreement(G, H, method="arith", verbose=False):
    """
    This function computes the GDD-agreement between two graphs G and H.

    Arguments
    ----------
    G, H : networkx.Graph
        Two networkx graphs.

    method : str
        * "arith" : arithmetic mean.
        * "geo" : geometric mean.
        * "vec" : vector before computing the GDD-agreement.

    kwargs : optional
        Parameters of function `_orbit_count()`.
    """

    assert (method=="arith") or (method=="geo") or (method=="vec")


    if verbose: print("Computing the GDD of the 1st graph.")
    gdd_G = distribution(G, verbose=verbose)

    if verbose: print("Computing the GDD of the 2nd graph.")
    gdd_H = distribution(H, verbose=verbose)


    gdd_diff = gdd_difference(gdd_G, gdd_H)


    if method == "arith":
        return gdd_diff.mean()

    elif method == "geo":
        return np.exp(np.log(gdd_diff).mean())

    else:
        return gdd_diff





def distribution(G, **kwargs):
    """
    This function computes the GDD of graph G.

    Arguments
    ----------
    G : networkx.Graph
        A networkx graph.

    kwargs : optional
        Parameters of function `_orbit_count()`.
    """

    counts = orbit_counter(G, **kwargs)

    return [ np.bincount(counts[:, j]) for j in range(NUMBER_OF_ORBITS) ]





def orbital_features(G, **kwargs):
    """
    This function computes the orbital features of nodes.

    Arguments
    ----------
    G : networkx.Graph
        A networkx graph.

    kwargs : optional
        Parameters of function `_orbit_count()`.
    """

    counts = orbit_counter(G, **kwargs)

    return dict(zip(G.nodes(), counts))





def orbit_counter(G, **kwargs):
    """
    This function counts the number of orbits for each node.

    Arguments
    ----------
    A : scipy.sparse.lil_matrix
        An adjacency matrix.

    kwargs : optional
        Parameters of function `_orbit_count()`.
    """

    assert (not G.is_directed())


    A = nx.to_scipy_sparse_matrix(G, format="lil")

    return _orbit_counter(A, **kwargs)





def _orbit_counter(A, verbose=False):
    """
    This function counts the number of orbits for each node.

    Arguments
    ----------
    A : scipy.sparse.lil_matrix
        An adjacency matrix.

    verbose : bool
        Show a progress bar if True.
    """

    def comb(n):
        return n * (n - 1) // 2


    G = nx.from_scipy_sparse_matrix(A)

    V = G.number_of_nodes()
    N = [ set(G.neighbors(v)) for v in G.nodes() ]


    counts = np.zeros((V, NUMBER_OF_ORBITS), dtype=np.int64)

    # Orbit 0 (degree distribution)
    counts[:, 0] = [ d for v, d in G.degree() ]


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

                if (r in Star_u) or (r in Star_v):
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



        counts[u, 1] += n_v
        counts[v, 1] += n_u
        counts[u, 2] += n_u
        counts[v, 2] += n_v
        counts[u, 3] += n_e
        counts[v, 3] += n_e

        counts[u,  4] += n_3_1_v
        counts[v,  4] += n_3_1_u
        counts[u,  5] += n_3_1_u + n_3_2
        counts[v,  5] += n_3_1_v + n_3_2
        counts[u,  6] += n_4_1_v
        counts[v,  6] += n_4_1_u
        counts[u,  7] += n_4_1_u
        counts[v,  7] += n_4_1_v
        counts[u,  8] += n_5_1
        counts[v,  8] += n_5_1
        counts[u,  9] += n_6_1_v
        counts[v,  9] += n_6_1_u
        counts[u, 10] += n_6_2_v + n_6_3
        counts[v, 10] += n_6_2_u + n_6_3
        counts[u, 11] += n_6_1_u + n_6_2_u
        counts[v, 11] += n_6_1_v + n_6_2_v
        counts[u, 12] += n_7_1_v
        counts[v, 12] += n_7_1_u
        counts[u, 13] += n_7_1_u + n_7_2
        counts[v, 13] += n_7_1_v + n_7_2
        counts[u, 14] += n_8_1
        counts[v, 14] += n_8_1

        counts[u, 15] += n_9_1_v
        counts[v, 15] += n_9_1_u
        counts[u, 16] += n_9_1_u + n_9_2_v
        counts[v, 16] += n_9_1_v + n_9_2_u
        counts[u, 17] += n_9_2_u
        counts[v, 17] += n_9_2_v
        counts[u, 18] += n_10_1_v
        counts[v, 18] += n_10_1_u
        counts[u, 19] += n_10_3_v
        counts[v, 19] += n_10_3_u
        counts[u, 20] += n_10_1_u + n_10_2_v
        counts[v, 20] += n_10_1_v + n_10_2_u
        counts[u, 21] += n_10_2_u + n_10_3_u
        counts[v, 21] += n_10_2_v + n_10_3_v
        counts[u, 22] += n_11_1_v
        counts[v, 22] += n_11_1_u
        counts[u, 23] += n_11_1_u
        counts[v, 23] += n_11_1_v
        counts[u, 24] += n_12_1_v
        counts[v, 24] += n_12_1_u
        counts[u, 25] += n_12_2_u
        counts[v, 25] += n_12_2_v
        counts[u, 26] += n_12_1_u + n_12_2_v + n_12_3
        counts[v, 26] += n_12_1_v + n_12_2_u + n_12_3
        counts[u, 27] += n_13_1_v
        counts[v, 27] += n_13_1_u
        counts[u, 28] += n_13_1_u + n_13_2_v
        counts[v, 28] += n_13_1_v + n_13_2_u
        counts[u, 29] += n_13_3_v + n_13_4
        counts[v, 29] += n_13_3_u + n_13_4
        counts[u, 30] += n_13_2_u + n_13_3_u
        counts[v, 30] += n_13_2_v + n_13_3_v
        counts[u, 31] += n_14_1_v
        counts[v, 31] += n_14_1_u
        counts[u, 32] += n_14_2_v + n_14_3
        counts[v, 32] += n_14_2_u + n_14_3
        counts[u, 33] += n_14_1_u + n_14_2_u
        counts[v, 33] += n_14_1_v + n_14_2_v
        counts[u, 34] += n_15_1
        counts[v, 34] += n_15_1
        counts[u, 35] += n_16_1_v
        counts[v, 35] += n_16_1_u
        counts[u, 36] += n_16_3_v
        counts[v, 36] += n_16_3_u
        counts[u, 37] += n_16_2_u + n_16_3_u
        counts[v, 37] += n_16_2_v + n_16_3_v
        counts[u, 38] += n_16_1_u + n_16_2_v
        counts[v, 38] += n_16_1_v + n_16_2_u
        counts[u, 39] += n_17_1_v
        counts[v, 39] += n_17_1_u
        counts[u, 40] += n_17_2_u + n_17_3_u
        counts[v, 40] += n_17_2_v + n_17_3_v
        counts[u, 41] += n_17_3_v + n_17_4_v
        counts[v, 41] += n_17_3_u + n_17_4_u
        counts[u, 42] += n_17_1_u + n_17_2_v + n_17_4_u
        counts[v, 42] += n_17_1_v + n_17_2_u + n_17_4_v
        counts[u, 43] += n_18_1   + n_18_2_v
        counts[v, 43] += n_18_1   + n_18_2_u
        counts[u, 44] += n_18_2_u
        counts[v, 44] += n_18_2_v
        counts[u, 45] += n_19_1_v
        counts[v, 45] += n_19_1_u
        counts[u, 46] += n_19_3_v
        counts[v, 46] += n_19_3_u
        counts[u, 47] += n_19_1_u + n_19_2_v
        counts[v, 47] += n_19_1_v + n_19_2_u
        counts[u, 48] += n_19_2_u + n_19_3_u + n_19_4
        counts[v, 48] += n_19_2_v + n_19_3_v + n_19_4
        counts[u, 49] += n_20_1_u
        counts[v, 49] += n_20_1_v
        counts[u, 50] += n_20_1_v
        counts[v, 50] += n_20_1_u
        counts[u, 51] += n_21_1   + n_21_2_v
        counts[v, 51] += n_21_1   + n_21_2_u
        counts[u, 52] += n_21_3_v
        counts[v, 52] += n_21_3_u
        counts[u, 53] += n_21_2_u + n_21_3_u + n_21_4
        counts[v, 53] += n_21_2_v + n_21_3_v + n_21_4
        counts[u, 54] += n_22_1_u
        counts[v, 54] += n_22_1_v
        counts[u, 55] += n_22_1_v + n_22_2
        counts[v, 55] += n_22_1_u + n_22_2
        counts[u, 56] += n_23_1_v
        counts[v, 56] += n_23_1_u
        counts[u, 57] += n_23_2_u + n_23_3
        counts[v, 57] += n_23_2_v + n_23_3
        counts[u, 58] += n_23_1_u + n_23_2_v
        counts[v, 58] += n_23_1_v + n_23_2_u
        counts[u, 59] += n_24_1_v + n_24_2_v
        counts[v, 59] += n_24_1_u + n_24_2_u
        counts[u, 60] += n_24_2_u + n_24_3   + n_24_4_v
        counts[v, 60] += n_24_2_v + n_24_3   + n_24_4_u
        counts[u, 61] += n_24_1_u + n_24_4_u
        counts[v, 61] += n_24_1_v + n_24_4_v
        counts[u, 62] += n_25_1_v
        counts[v, 62] += n_25_1_u
        counts[u, 63] += n_25_1_u + n_25_2_v
        counts[v, 63] += n_25_1_v + n_25_2_u
        counts[u, 64] += n_25_2_u + n_25_3
        counts[v, 64] += n_25_2_v + n_25_3
        counts[u, 65] += n_26_1_v
        counts[v, 65] += n_26_1_u
        counts[u, 66] += n_26_2_u + n_26_3
        counts[v, 66] += n_26_2_v + n_26_3
        counts[u, 67] += n_26_1_u + n_26_2_v + n_26_4
        counts[v, 67] += n_26_1_v + n_26_2_u + n_26_4
        counts[u, 68] += n_27_1   + n_27_2_v
        counts[v, 68] += n_27_1   + n_27_2_u
        counts[u, 69] += n_27_2_u
        counts[v, 69] += n_27_2_v
        counts[u, 70] += n_28_1_v
        counts[v, 70] += n_28_1_u
        counts[u, 71] += n_28_1_u + n_28_2
        counts[v, 71] += n_28_1_v + n_28_2
        counts[u, 72] += n_29_1
        counts[v, 72] += n_29_1


    # Deleting duplicates

    targets = [
         2,  3,  5,  8, 10, 12, 16, 17, 20, 25,
        28, 29, 32, 34, 36, 37, 40, 43, 46, 49,
        51, 52, 54, 59, 62, 65
    ]

    for j in targets:
        counts[:, j] //= 2

    targets = [
         7, 11, 13, 14, 21, 26, 30, 38, 41, 47,
        48, 50, 53, 57, 60, 63, 64, 66, 68, 70
    ]

    for j in targets:
        counts[:, j] //= 3

    targets = [23, 33, 42, 44, 55, 58, 61, 67, 69, 71, 72]

    for j in targets:
        counts[:, j] //= 4


    return counts





def gdd_difference(gdd_G, gdd_H):
    """
    This function computes the difference between two GDDs gdd_G and gdd_H.

    Arguments
    ----------
    gdd_G, gdd_H : list of numpy.ndarray
        The GDDs.
    """

    sq2 = sqrt(2)

    norm_G = normalized_gdd(gdd_G)
    norm_H = normalized_gdd(gdd_H)

    return np.array([ 1 - padded_norm(norm_G[j], norm_H[j]) / sq2 for j in range(NUMBER_OF_ORBITS) ])



def normalized_gdd(gdd):
    # Scaling
    S = [ np.array([ d / k if k != 0 else 0 for k, d in enumerate(dd_j) ])  for dd_j in gdd ]

    # Total
    T = [ s_j.sum() for s_j in S ]

    return [ s_j / T[j] if T[j] != 0 else np.array([0.0]) for j, s_j in enumerate(S) ]



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
