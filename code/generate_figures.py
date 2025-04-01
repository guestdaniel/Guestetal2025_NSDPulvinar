# generate_figures.py
#
# This script generates all Python-derived figure elements for the manuscript.
# This script does fairly little heavy lifting, and as such can be expected
# to run on a normal workstation and should finish within 30 minutes.

# Import OS and navigate to correct directory to run this file; we assume that the script
# is initiated from the main repository folder!
import os
os.chdir("code")

# Import functions from various figure function files
from figs_gen.colorbars import *
from figs_gen.fig_prf_contrast import *
from figs_gen.fig_prf_body import *
from figs_gen.fig_prf_wta import *
# from figs_gen.fig_corr_methods import *
# from figs_gen.fig_corr_cor_to_sub import *
# from figs_gen.fig_corr_sub_to_cor import *
# from figs_gen.supp_fig_prf_wta import *

# Function to generate all colorbars
def run_colorbars():
    colorbar_feature()
    colorbar_angle()
    colorbar_size()
    colorbar_variance_explained()
    colorbar_roi()
    colorbar_laterality()
    colorbar_pearsons()
    colorbar_consistency()

# Function to generate components of Figure 2 (winner-take-all pRF comparison)
def run_fig2():
    fig_prf_wta_contextual_anatomy()
    fig_prf_wta_mip_maps()
    fig_prf_wta_maps()

# Function to generate components of Figure 3 (pRF spatial coding)
def run_fig3():
    fig_prf_contrast_anatomy()
    fig_prf_contrast_maps()
    fig_prf_contrast_rf_coverage()
    fig_prf_body_maps()
    fig_prf_body_rf_coverage()

# def run_fig6():
#     plot_fig_corr_cor_to_sub_ventral_stream_avg_group_contours()
#     plot_fig_corr_cor_to_sub_individual_subject_consistency_maps()
#     plot_fig_corr_cor_to_sub_group_consistency_maps()

# def run_figsupp():
#     supp_fig_prf_wta.supp_fig_prf_wta_rainbow_maps()
#     supp_fig_prf_wta.supp_fig_prf_wta_correlation_matrix()


if __name__ == '__main__':
    # Generate colorbars
    run_colorbars()

    # Run figures in order
    run_fig2()
    run_fig3()
    # run_fig6()
    # run_figsupp()