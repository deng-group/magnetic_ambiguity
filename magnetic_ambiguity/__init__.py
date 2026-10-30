"""Geometry-magnetic ambiguity of DFT datasets.

Site ambiguity  a_i = median_j |m_i - m_j|  over the adaptive SOAP twins of site i
(exact K nearest neighbours in L2-normalised SOAP space, twins d_ij <= alpha * d_i,m).
Dataset ambiguity T_p = mean of the upper-p mass of the a_i distribution (default p = 8 %).

The defaults reproduce the definition in force since 2026-10-28 ("ambiguity = T8"):
SOAP r_cut 6 A, n_max 8, l_max 3; K = 100; alpha = 1.25; m = 5; same-frame sites are
candidates and only the query site is excluded; signed moments; frame-force filter
|F_DFT| > 20 eV/A on query sites; p = 0.08.  Every one of these is an option.
"""

from .settings import (AmbiguitySettings, DatasetSettings, NeighborSettings,
                       PopulationSettings, Settings, SoapSettings)
from .dataset_ambiguity import dataset_statistics, tail, tail_mean
from .neighbors import adaptive_k, exact_knn
from .site_ambiguity import prefix_median, site_ambiguity

__version__ = "0.1.0"
__all__ = [
    "AmbiguitySettings", "DatasetSettings", "NeighborSettings", "PopulationSettings",
    "Settings", "SoapSettings", "dataset_statistics", "tail", "tail_mean",
    "adaptive_k", "exact_knn", "prefix_median", "site_ambiguity", "__version__",
]
