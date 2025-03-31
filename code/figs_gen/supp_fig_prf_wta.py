import os
import numpy as np
import nibabel as nib
import matplotlib
import matplotlib.pyplot as plt
matplotlib.use('Agg')
from figure_funcs import ff


def supp_fig_prf_wta_rainbow_maps():
    # Load in all data and store it in a list in order
    data = list()
    maps = ['contrastNEW', 'foregroundauto', 'faceauto', 'bodyauto', 'backgroundauto', 'wordauto', 'salience']
    for map in maps:
        data.append(np.mean(np.array(ff.load_volume(volume=map + '_R2')), axis=0))
    T1 = np.mean(np.array(ff.load_volume(volume='T1')), axis=0)
    THA = nib.load(os.path.join(ff.dir_data, 'group', 'mni', 'postthalamus.nii.gz')).get_fdata()

    # Pick maps and colors
    selected_maps = ['foregroundauto', 'contrastNEW', 'bodyauto']
    colors = np.array([[255, 89, 89],
                       [122, 253, 121],
                       [56, 143, 255]])/255

    # Stack all selected maps into one array
    map = np.stack([data[maps.index(selected_maps[0])][:, :, :],
                    data[maps.index(selected_maps[1])][:, :, :],
                    data[maps.index(selected_maps[2])][:, :, :]], axis=3)

    # Clean up the maps and construct the colormap
    map[map < 0] = 0
    map = map/np.nanmax(map[THA == 1])
    map = 1-(1-np.tile(np.expand_dims(map[:, :, :, 0], axis=3), 3)*colors[0])*(1-np.tile(np.expand_dims(map[:, :, :, 1], axis=3), 3)*colors[1])*(1-np.tile(np.expand_dims(map[:, :, :, 2], axis=3), 3)*colors[2])
    map = map*255*2
    map[map > 255] = 255
    map = map.astype('uint8')

    # construct the figure
    slices = [100, 97, 95, 93]
    plt.figure(figsize=(12, 2.5))
    n_slice = len(slices)
    cutoff = 0.2
    for idx_slice, slice in enumerate(slices):
        # R2
        plt.subplot(1, n_slice, idx_slice+1)
        ff.plot_rainbow_slice(T1, map, slice, np.zeros(THA.shape), (0 + 50, 182 - 50), (0 + 55, 182 - 80))

    plt.tight_layout()
    plt.savefig(os.path.join('../figures', 'supp_fig_prf_wta_rainbow_maps.png'), dpi=300)


def supp_fig_prf_wta_correlation_matrix():
    # Load in all data and store it in a list in order
    data = list()
    maps = ['contrastNEW', 'salience', 'foregroundauto', 'bodyauto', 'faceauto']  # discard word and background
    for map in maps:
        data.append(ff.load_volume(volume=map + '_R2'))
    THA = nib.load(os.path.join(ff.dir_data, 'group', 'mni', 'postthalamus.nii.gz')).get_fdata()

    # Compute WTA map
    R2_max = np.max(np.mean(np.array(data), axis=1), axis=0)

    # Plot function
    def plot_ND_histogram(data, mask, bins_1=np.linspace(0, 8, 16), bins_2=np.linspace(0, 8, 16)):
        f, axs = plt.subplots(len(data)+1, len(data)+1, gridspec_kw={'width_ratios': [2] + [6]*len(data), 'height_ratios': [6]*len(data) + [2.5]}, figsize=(8, 8))
        # Vertical histograms histogram
        for idx_subdata, subdata in enumerate(data):
            axs[idx_subdata][0].hist(np.concatenate([dat[mask] for dat in subdata]), orientation='horizontal', bins=bins_1)
            axs[idx_subdata][0].get_xaxis().set_visible(False)
            axs[idx_subdata][0].set_ylim((np.min(bins_1), np.max(bins_1)))
        # Horizontal histograms histogram
        for idx_subdata, subdata in enumerate(data):
            axs[-1][idx_subdata+1].hist(np.concatenate([dat[mask] for dat in subdata]), bins=bins_2)
            axs[-1][idx_subdata+1].get_yaxis().set_visible(False)
            axs[-1][idx_subdata+1].set_xlim((np.min(bins_2), np.max(bins_2)))
        # 2D histograms
        for idx_subdata1, subdata1 in enumerate(data):
            for idx_subdata2, subdata2 in enumerate(data):
                if idx_subdata1 <= idx_subdata2:
                    axs[idx_subdata1][idx_subdata2+1].hist2d(np.concatenate([dat[mask] for dat in subdata2]), np.concatenate([dat[mask] for dat in subdata1]),
                              bins=(bins_1, bins_2))
                    axs[idx_subdata1][idx_subdata2+1].get_xaxis().set_visible(False)
                    axs[idx_subdata1][idx_subdata2+1].get_yaxis().set_visible(False)
                else:
                    x1 = np.concatenate([dat[mask] for dat in subdata2])
                    x2 = np.concatenate([dat[mask] for dat in subdata1])
                    include = np.logical_and(x1 > 0, x2 > 0)
                    x1 = x1[include]
                    x2 = x2[include]
                    corrcoef = np.corrcoef(x1, x2)[1, 0]
                    axs[idx_subdata1][idx_subdata2+1].get_xaxis().set_visible(False)
                    axs[idx_subdata1][idx_subdata2+1].get_yaxis().set_visible(False)
                    axs[idx_subdata1][idx_subdata2+1].text(0.75/2 + 0.1, 0.75/2, str(round(corrcoef, 2)), horizontalalignment='center', fontdict={'weight': 'bold', 'size': 16}, color=matplotlib.cm.get_cmap('RdBu')(0.5+corrcoef*2))
                    axs[idx_subdata1][idx_subdata2+1].spines['right'].set_visible(False)
                    axs[idx_subdata1][idx_subdata2+1].spines['top'].set_visible(False)
                    axs[idx_subdata1][idx_subdata2 + 1].spines['bottom'].set_visible(False)
                    axs[idx_subdata1][idx_subdata2 + 1].spines['left'].set_visible(False)
        # Hide spare axes
        axs[-1][0].set_visible(False)

    # Save
    plot_ND_histogram(data, np.logical_and(R2_max > 0.2, THA == 1), bins_1=np.linspace(0.0, 0.75, 25), bins_2=np.linspace(0.0, 0.75, 25))
    plt.tight_layout(h_pad=0.0, w_pad=0.0)
    plt.savefig(os.path.join('../figures', 'supp_fig_prf_wta_correlation_matrix.png'), dpi=300)