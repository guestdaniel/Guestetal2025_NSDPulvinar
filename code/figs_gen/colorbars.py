import os
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import matplotlib as mpl
matplotlib.use('Agg')
from figure_funcs import ff


def colorbar_feature():
    # Save dummy colorbars
    a = np.array([[0, 1]])
    plt.figure(figsize=(2, 2))
    img = plt.imshow(a, cmap=ff.cmap_feat.reversed())
    plt.gca().set_visible(False)
    cax = plt.axes([0.05, 0.05, 0.15, 0.9])
    cbar = plt.colorbar(cax=cax)
    cbar.set_ticks(np.linspace(0.05, 0.95, 5))
    cbar.set_ticklabels(['Faces', 'Bodies', 'Foreground', 'Salience', 'Contrast'])
    plt.savefig(os.path.join('../figures', 'fig_colorbar_features.png'), dpi=300)


def colorbar_variance_explained():
    # Variance explained colorbar
    a = np.array([[0, 1]])
    plt.figure(figsize=(1, 2.25))
    img = plt.imshow(a, cmap=ff.cmap_rsqr)
    plt.gca().set_visible(False)
    cax = plt.axes([0.2, 0.1, 0.2, 0.8])
    plt.colorbar(cax=cax)
    plt.savefig(os.path.join('../figures', 'fig_colorbar_variance_explained.png'))


def colorbar_angle():
    # Angle colorbar
    fig = plt.figure()
    display_axes = fig.add_axes([0.1, 0.1, 0.8, 0.8], projection='polar')
    quant_steps = 512
    norm = mpl.colors.Normalize(0, 2 * np.pi)
    hsv = mpl.cm.get_cmap('hsv', quant_steps)
    #cmap = mpl.colors.ListedColormap(ff.cmap_angle(np.roll(np.tile(np.linspace(0, 1, quant_steps), 2), int(quant_steps/2))))
    cmap = mpl.colors.ListedColormap(ff.cmap_angle(np.roll(np.concatenate([np.linspace(1, 0, quant_steps), np.linspace(0, 1, quant_steps)]), int(quant_steps/2))))

    cb = mpl.colorbar.ColorbarBase(display_axes,
                                   cmap=cmap,
                                   norm=norm,
                                   orientation='horizontal')
    cb.outline.set_visible(False)
    display_axes.set_axis_off()
    plt.plot([0, np.pi/2], [0, 12], color='white')
    plt.plot([0, -np.pi/2], [0, 12], color='white')
    plt.savefig(os.path.join('../figures', 'fig_colorbar_angle.png'))


def colorbar_size():
    # Size/eccentricity colorbar
    a = np.array([[0, 6]])
    plt.figure(figsize=(1, 2.25))
    img = plt.imshow(a, cmap=ff.cmap_size)
    plt.gca().set_visible(False)
    cax = plt.axes([0.2, 0.1, 0.2, 0.8])
    plt.colorbar(cax=cax)
    plt.savefig(os.path.join('../figures', 'fig_colorbar_size.png'))
    plt.close('all')


def colorbar_laterality():
    a = np.array([[-1, 1]])
    plt.figure(figsize=(1.25, 2.25))
    img = plt.imshow(a, cmap=ff.cmap_lat)
    plt.gca().set_visible(False)
    cax = plt.axes([0.2, 0.12, 0.2, 0.78])
    cbar = plt.colorbar(cax=cax)
    cbar.set_ticks([-1, 0, 1])
    cbar.set_ticklabels(['', '', ''])
    plt.savefig(os.path.join('../figures', 'fig_colorbar_lat.png'))
    plt.close('all')


def colorbar_roi():
    # Save dummy colorbars
    a = np.array([[0, 1]])
    plt.figure(figsize=(2, 2))
    img = plt.imshow(a, cmap=ff.cmap_thal.reversed())
    plt.gca().set_visible(False)
    cax = plt.axes([0.05, 0.05, 0.15, 0.9])
    cbar = plt.colorbar(cax=cax)
    cbar.set_ticks(np.linspace(0.10, 0.90, 5))
    cbar.set_ticklabels(['SC', 'dmPul', 'dlPul', 'vPul', 'LGN'])
    plt.savefig(os.path.join('../figures', 'fig_colorbar_thalamus.png'), dpi=300)


def colorbar_pearsons():
    a = np.array([[0, 1]])
    plt.figure(figsize=(1, 2.25))
    img = plt.imshow(a, cmap=ff.cmap_rsqr)
    plt.gca().set_visible(False)
    cax = plt.axes([0.15, 0.1, 0.2, 0.85])
    cbar = plt.colorbar(cax=cax)
    cbar.set_ticks([0, 0.5, 1])
    cbar.set_ticklabels(['0.0', '0.05', '0.1'])
    plt.savefig(os.path.join('../figures', 'fig_colorbar_pearsons.png'))


def colorbar_consistency():
    a = np.array([[0, 1]])
    plt.figure(figsize=(1, 2.25))
    img = plt.imshow(a, cmap=matplotlib.colors.ListedColormap(ff.cmap_cons(np.linspace(0, 1, 9))))
    plt.gca().set_visible(False)
    cax = plt.axes([0.15, 0.1, 0.2, 0.85])
    cbar = plt.colorbar(cax=cax)
    cbar.set_ticks(np.linspace(0+0.05, 1-0.05, 9))
    cbar.set_ticklabels(['0', '', '2', '', '4', '', '6', '', '8'])
    plt.savefig(os.path.join('../figures', 'fig_colorbar_consistency.png'))
