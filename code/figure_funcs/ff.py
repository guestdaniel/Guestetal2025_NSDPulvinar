import numpy as np
import matplotlib.pyplot as plt
import matplotlib
from skimage import measure
import nibabel as nib
import os
import scipy
import random
import math

# Set location for data (we expect these functions to be executed in the `code` directory,
# and here we define the other paths relative to that directory)
dir_data = "../data/prepared"
dir_plots = "../figures"

# Handle some general sizing across figures
SMALL_SIZE = 18
MEDIUM_SIZE = 22
BIGGER_SIZE = 20
plt.rc('font', size=SMALL_SIZE)          # controls default text sizes
plt.rc('axes', titlesize=SMALL_SIZE)     # fontsize of the axes title
plt.rc('axes', labelsize=MEDIUM_SIZE)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('ytick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('legend', fontsize=SMALL_SIZE)    # legend fontsize
plt.rc('figure', titlesize=BIGGER_SIZE)  # fontsize of the figure title

# Set up several colormaps
def create_cmap_arcaro():
    """ Creates a colormap similar to that from Arcaro et al. (2015)

    Returns:
        cmap_arcaro (ListedColormap): ListedColormap object containing the colormap
    """
    cmap_arcaro = np.array([[255, 0, 0],
                            [255, 68, 0],
                            [255, 153, 0],
                            [255, 255, 0],
                            [0, 255, 0],
                            [50, 205, 50],
                            [0, 255, 255],
                            [0, 153, 255],
                            [0, 105, 255],
                            [0, 68, 255]])
    cmap_arcaro = cmap_arcaro / 255
    cmap_arcaro = np.append(cmap_arcaro, np.ones((10, 1)), 1)
    cmap_arcaro = np.flip(cmap_arcaro, axis=0)
    cmap_arcaro = matplotlib.colors.ListedColormap(cmap_arcaro)
    return cmap_arcaro

# colormap for all angle maps
cmap_angle = create_cmap_arcaro()

# colormap for variance explained / R^2 maps
cmap_rsqr = matplotlib.cm.get_cmap('hot')
cmap_r = matplotlib.cm.get_cmap('RdBu_r')

# colormap for size and eccentricity maps
cmap_size = matplotlib.cm.get_cmap('plasma_r')

# colormap for lateralization maps
cmap_lat  = matplotlib.cm.get_cmap('RdBu_r')

# colormap for consistency maps
cmap_cons = matplotlib.cm.get_cmap('winter')

# colormap for thalamus ROIs (LGN, pul, and SC)
cmap_thal = matplotlib.colors.ListedColormap([[0.0, 0.5, 0.4, 1.0],
                                             [0.3, 0.8, 1.0, 1.0],
                                             [0.2, 0.3, 0.7, 1.0],
                                             [0.4, 0.2, 0.7, 1.0],
                                             [0.7, 0.1, 0.6, 1.0]])                                                                     

# colormap for features (contrast, bodies, etc.)
cmap_feat = matplotlib.colors.ListedColormap([[0.8, 0.2, 0.6, 1.0],
                                              [0.8, 0.1, 0.2, 1.0],
                                              [0.9, 0.4, 0.0, 1.0],
                                              [0.9, 0.8, 0.0, 1.0],
                                              [0.1, 0.6, 0.2, 1.0]])
cmap_feat2 = matplotlib.colors.ListedColormap([[27/255, 158/255, 119/255, 1.0],
                                              [217/255, 95/255, 2/255, 1.0],
                                              [117/255, 112/255, 179/255, 1.0],
                                              [231/255, 41/255, 138/255, 1.0],
                                              [230/255, 171/255, 2/255, 1.0]])

cmap_corr = matplotlib.colors.ListedColormap([[0.0, 0.6, 0.5, 1.0],
                                              [0.0, 0.9, 0.5, 1.0],
                                              [0.4, 0.8, 0.3, 1.0],
                                              [0.9, 0.8, 0.0, 1.0],
                                              [0.0, 0.6, 0.8, 1.0],
                                              [0.1, 0.4, 0.8, 1.0]])


def fibonacci_sphere(samples=1000):
    """ Quick and dirty function to evenly sample points along sphere surface

    Directly copied from https://stackoverflow.com/questions/9600801/evenly-distributing-n-points-on-a-sphere
    """
    points = []
    phi = math.pi * (math.sqrt(5.) - 1.)  # golden angle in radians

    for i in range(samples):
        y = 1 - (i / float(samples - 1)) * 2  # y goes from 1 to -1
        radius = math.sqrt(1 - y * y)  # radius at y

        theta = phi * i  # golden angle increment

        x = math.cos(theta) * radius
        z = math.sin(theta) * radius

        points.append((x, y, z))

    return points


def load_volume(subjs=None, space='mni', volume='T1'):
    """ Load NIFTI data of multiple subjects from Guestetal2021_data

    Args:
        subjs (list): List of integers indicating which subjects' data to load. If None, defaults to all eight subjects.
            The contents of the list are not explicitly checked, and strange inputs will raised equally strange errors.
        space (string): a space to load from, either 'mni' or 'func1mm'.
        volume (string): name of datafile to load (e.g., 'T1')
    """
    # Set default subjs
    if subjs is None:
        subjs = [1, 2, 3, 4, 5, 6, 7, 8]
    elif not isinstance(subjs, list):
        raise TypeError('subject should be a list of integers!')
    # Handle space
    if space not in ['mni', 'func1mm']:
        raise ValueError('space must be either mni or func1mm')
    # Loop over subjs and load data
    storage = []
    for subj in subjs:
        storage.append(nib.load(os.path.join(dir_data, 'subj0' + str(subj), space, volume + '.nii.gz')).get_fdata())
    # Return
    return storage


def calc_areas(angle, size, ecc):
    # Convert angle to radians
    angle = np.pi/180 * angle
    # Calculate x, y, and radius
    x = np.cos(angle) * ecc
    y = np.sin(angle) * ecc
    radius = size*2
    # Calculate r, the distance from the center to the midline
    r = np.abs(x)
    # Calculate h, the distance from the edge of the circle in towards the midline
    h = radius - np.abs(r)
    # Calculate total area
    A = 2*np.pi*radius**2
    # If h is less than or equal to zero, that means the circle is entirely on one side... return now
    if h <= 0:
        if x >= 0:
            return 1
        else:
            return -1
    # Calculate R, the length of the line going from the midpoint to the intersection of the midline and circle
    R = h + r
    # Calculate a, the length of the segment
    a = 2 * np.sqrt(R**2 - r**2)
    # Calculate theta, the angle of the triangle formed by R and a
    theta = 2 * np.arctan2(a, 2*r)
    # Calculate the area of the sector and the area of the triangle
    A_sector = A * theta / (2 * np.pi)
    A_triangle = a * r
    # Calculate final area of the segment
    A_segment = A_sector - A_triangle
    if x >= 0:
        return ((A-A_segment)/A - 1/2)*2
    else:
        return (A_segment/A - 1/2)*2


def calc_lat(angle, size, ecc):
    """ Calculates a laterality metric for a set of angle, size, and eccentricity RF measurements

    Args:
        angle (ndarray): angle measurements
        size (ndarray): size measurements
        ecc (ndarray): eccentricity measurements

    Returns:
        lat (ndarray): laterality metric
    """
    lat = np.zeros_like(angle)
    for ii in range(angle.shape[0]):
        for jj in range(angle.shape[1]):
            for kk in range(angle.shape[2]):
                lat[ii, jj, kk] = calc_areas(angle[ii, jj, kk], size[ii, jj, kk], ecc[ii, jj, kk])
    return lat


def quick_bounding_box_viz(T1, THA):
    """ Some quick and dirty code to plot T1 and THA
    """
    plt.figure(figsize=(20, 20))
    slices = np.arange(106, 89, -1)
    vmax = 0.8 * np.max(T1)
    for idx_slice, slice in enumerate(slices):
        plt.subplot(int(np.ceil(np.sqrt(len(slices)))), int(np.ceil(np.sqrt(len(slices)))), idx_slice+1)
        plt.imshow(T1[:, slice, :].transpose(), cmap='gray', vmin=0, vmax=vmax, origin='lower')
        plt.imshow(np.ma.masked_where(THA[:, slice, :].transpose() <= 0, THA[:, slice, :].transpose()), origin='lower')
        if np.mod(idx_slice, int(np.ceil(np.sqrt(len(slices))))) == 0:
            plt.xticks(np.arange(0, 182, 10))
            plt.yticks(np.arange(0, 182, 10))
        else:
            plt.xticks(np.arange(0, 182, 10))
            plt.yticks([])
        plt.xlim(0+55, 182-55)
        plt.ylim(0+55, 182-85)
        plt.title(str(slice))
    plt.gcf().tight_layout(pad=0, h_pad=0.1, w_pad=0.1)
    plt.show()


def char_prf_angle(an, r2, thal, slices=None, n_rep_shuffle=5, side='left', orientation='coronal'):
    """ Characterizes the organization of retinotopic maps in terms of their angle progression

    Analyzes angle data in `AN` to characterize how well it obeys an oriented linear
    gradient. We assume the angle data are from 0 to 180 degrees, where 0 degrees
    [[corresponds to the downward vertical meridian and 180 degrees corresponds to the
    upeward vertical meridian.]] For several coronal slices that cover our areas of
    interest, we correlate the angles against their voxel coordinates' projections onto
    vectors at different orientations. We compare the result for the data to corresponding
    simulated datasets generated by shuffling the data within each coronal slice.
    """
    # Note:
        # indices [0:93, :, :] -> LH
        # indices [94:, :, :] -> RH

    # Create empty figure
    fig = plt.figure(figsize=(12, 2.5), dpi=300)

    # Handle slices
    if slices is None:
        slices = np.arange(102, 93, -1)
    
    # Determine matrix indiceds that are valid (needs to be left hemisphere and variance
    # explained greater than 0.1% to be considered for testing)
    hemi = np.ones(an.shape, dtype=bool)
    if side == 'left':
        hemi[94:, :, :] = False
    else:
        hemi[1:93, :, :] = False
    idxs_valid = np.logical_and(r2 > 0.1, thal == 1)
    idxs_valid = np.logical_and(idxs_valid, hemi)

    # Determine vector orientations we want to test and represent as an (n_orientation, 2) 
    # matrix of x and y coordinates for each oriented vector
    angles = np.linspace(0.0, 2*np.pi, 100)
    x = np.cos(angles)
    y = np.sin(angles)
    vs = np.stack([x, y], axis=1)

    # Create matrix to store result
    results = np.zeros([len(angles), len(slices)])
    results_shuffled = np.zeros([len(angles), len(slices), n_rep_shuffle])

    # Loop through slices, for each, we want to project all angles onto unit vectors at
    # different angles 
    for i, slice in enumerate(slices):
        # First, for this slice, determine which of its flattened indices are valid
        if orientation == 'coronal':
            idxs_valid_flat = np.squeeze(idxs_valid[:, slice, :]).flatten()
        elif orientation == 'sagittal':
            idxs_valid_flat = np.squeeze(idxs_valid[slice, :, :]).flatten()
        else:
            idxs_valid_flat = np.squeeze(idxs_valid[:, :, slice]).flatten()

        # Next, squeeze the data in the slice and subset to include only valid elements. We
        # also determine and subset the cartesian coordintes of each element
        if orientation == 'coronal':
            ss = np.squeeze(an[:, slice, :])  # now, dim1 -> left-right, dim 2 -> inferior-superior
        elif orientation == 'sagittal':
            ss = np.squeeze(an[slice, :, :])  # now, dim1 -> posterior-anterior, dim 2 -> inferior-superior
        else:
            ss = np.squeeze(an[:, :, slice])  # now, dim1 -> left-right, dim 2 -> posterior-anterior

        I, J = np.indices(ss.shape)    
        idxs = np.stack([I.flatten(), J.flatten()], axis=1)
        idxs = idxs[idxs_valid_flat]
        ss = ss.flatten()[idxs_valid_flat]

        # Check to make sure ss is not empty
        if len(ss) == 0:
            continue

        # Now, for each possible orientation , we determine how well angle pRF params can be
        # explained by linear regression between angle and the projection of the voxel
        # coordinates onto that oriented vector
        curve = np.zeros(len(angles))
        for j, v in enumerate(vs):
            # Determine projections of coordinates
            v_mag = np.linalg.norm(v)
            scalar_projections = np.dot(idxs, v) / v_mag

            # Correlate projections with true angles
            _, _, r_value, _, _= scipy.stats.linregress(scalar_projections, ss)
            # curve[j] = np.arctan(r_value)/(1/np.sqrt(len(ss) - 3))
            curve[j] = r_value

        results[:, i] = curve

        # Now, we repeat the above process on shuffled copies of the data
        for idx_rep in range(n_rep_shuffle):
            idxs_shuffle = random.sample(range(len(ss)), len(ss))
            ss_shuffled = ss[idxs_shuffle]
            curve = np.zeros(len(angles))
            for j, v in enumerate(vs):
                # Determine projections of coordinates
                v_mag = np.linalg.norm(v)
                scalar_projections = np.dot(idxs, v) / v_mag

                # Correlate projections with true angles
                _, _, r_value, _, _= scipy.stats.linregress(scalar_projections, ss_shuffled)
                # curve[j] = np.arctan(r_value)/(1/np.sqrt(len(ss) - 3))
                curve[j] = r_value
            
            results_shuffled[:, i, idx_rep] = curve

    # Now, we loop through slices and plot each for one half of the unit circle
    # For each set of observed data in a slice, we plot a thick black line showing its 
    # Fisher r-to-Z transform versus map orientation; we plot all shuffled replicates in 
    # thin gray lines. The significance test pits the observed peak Z value versus the
    # distribution of shuffled peak Z values and tests one-sided at 95%.
    for idx_slice in range(len(slices)):
        # Create subplot for this plot
        plt.subplot(1, len(slices), idx_slice+1, projection='polar')

        # Test for signif
        Z_real = np.max(np.squeeze(results[:, idx_slice]))
        Z_shuffle = np.max(np.squeeze(results_shuffled[:, idx_slice, :]), 0)
        if Z_real > np.percentile(Z_shuffle, 95.0):
            signif = True
        else:
            signif = False

        # Plot data
        if signif:
            color = 'red'
        else:
            color = 'black'
        plt.plot(angles, results[:, idx_slice] ** 2, color=color)

        # Plot shuffled replicates
        plt.plot(angles, np.squeeze(results_shuffled[:, idx_slice, :]) ** 2, color='gray', alpha=0.05, linewidth=0.2)
        
        # Set limits and ticks and labels
        plt.xlim(0.0, 2*np.pi)
        plt.yticks(np.arange(0.0, 0.5, 0.1), [])
        plt.xticks([0, np.pi/2, np.pi, 2*np.pi/3], [])
        plt.ylim(0.0, 0.5)
        plt.gca().set_thetamin(0)
        plt.gca().set_thetamax(360.0)
        plt.gca().set_rmax(0.4)
        plt.gca().set_rmin(0) 
        plt.gca().set_rorigin(0) 
        plt.gca().set_rticks(np.arange(0.0, 0.5, 0.1))
        plt.gca().set_thetagrids(np.linspace(0.0, 270.0, 4))

        # Add title
        plt.title(str(slices[idx_slice]))

    fig.tight_layout()
    plt.subplots_adjust(wspace=0, hspace=0)

    # Last step: print R2 max and corresponding angle for each slice
    for idx_slice in range(len(slices)):
        r2_max = np.max(results[:, idx_slice]) ** 2
        angle_max = angles[np.argmax(results[:, idx_slice])]
        print("Slice " + str(slices[idx_slice]) + ": r^2 = " + str(round(r2_max, 3)) + " at " + str(round(np.rad2deg(angle_max))))


def plot_sphere():
    azimuths = np.linspace(0.0, 2*np.pi, 200)
    elevations = np.linspace(0.0, 1*np.pi, 200)
    theta, phi = np.meshgrid(azimuths, elevations)
    R = np.abs(np.sin(2*np.pi * 3*theta))
    x = R * np.sin(phi) * np.cos(theta)
    y = np.sin(phi) * np.sin(theta)
    z = np.cos(phi)
    vs = np.stack([x, y, z], axis=1)
    n_vec = len(azimuths) * len(elevations)

    fig = plt.figure()
    ax = fig.add_subplot(projection='3d')
    ax.plot_surface(x, y, z, linewidth = 0.5, edgecolors = 'k')
    plt.show()


def char_prf_angle_3D_plot_one_angle(an, r2, thal, v, side='left', r2_cutoff=0.1):
    """ Characterizes the organization of pRF maps angle in 3D space
    """
    # Note:
        # indices [0:93, :, :] -> LH
        # indices [94:, :, :] -> RH

    # Create empty figure
    fig = plt.figure(figsize=(5.0, 5.0), dpi=300)
    
    # Determine matrix indiceds that are valid (needs to be left hemisphere and variance
    # explained greater than 0.1% to be considered for testing)
    hemi = np.ones(an.shape, dtype=bool)
    if side == 'left':
        hemi[94:, :, :] = False
    else:
        hemi[1:93, :, :] = False
    idxs_valid = np.logical_and(np.logical_and(r2 > r2_cutoff, thal == 1), hemi)

    # First, determine all linear indices of valid coordinates
    idxs_valid_flat = idxs_valid.flatten()

    # Next, squeeze the data in the slice and subset to include only valid elements. We
    # also determine and subset the cartesian coordintes corresponding to each element
    ss = an.flatten()[idxs_valid_flat]  # flattened subset of data we want 
    I, J, K = np.indices(r2.shape)    
    idxs = np.stack([I.flatten(), J.flatten(), K.flatten()], axis=1)
    idxs = idxs[idxs_valid_flat]        # flattened cartesian coordinates 

    # Determine projections of coordinates
    v_mag = np.linalg.norm(v)
    scalar_projections = np.dot(idxs, v) / v_mag

    # Correlate projections with true angles
    _, _, r_value, _, _= scipy.stats.linregress(scalar_projections, ss)
    plt.scatter(scalar_projections, ss)
    plt.show()


def char_prf_angle_3D(an, r2, thal, n_rep_shuffle=5, side='left', view='coronal', n_orientation=100, r2_cutoff=0.2):
    """ Characterizes the organization of pRF maps angle in 3D space
    """
    # Note:
        # indices [0:93, :, :] -> LH
        # indices [94:, :, :] -> RH

    # Determine matrix indiceds that are valid (needs to be left hemisphere and variance
    # explained greater than 0.1% to be considered for testing)
    hemi = np.ones(an.shape, dtype=bool)
    if side == 'left':
        hemi[94:, :, :] = False
    else:
        hemi[1:93, :, :] = False
    idxs_valid = np.logical_and(np.logical_and(r2 > r2_cutoff, thal == 1), hemi)

    # From azimuths and elevations, we first need to calculate the vectors that represent
    # all of the candidate pRF map orientations we want to test
    azimuths = np.linspace(0.0, 2*np.pi, n_orientation)
    elevations = np.linspace(0.0, np.pi, n_orientation)
    theta, phi = np.meshgrid(azimuths, elevations)
    v_X = np.sin(phi) * np.cos(theta)
    v_Y = np.sin(phi) * np.sin(theta)
    v_Z = np.cos(phi)

    # First, determine all linear indices of valid coordinates
    idxs_valid_flat = idxs_valid.flatten()

    # Next, squeeze the data in the slice and subset to include only valid elements. We
    # also determine and subset the cartesian coordintes corresponding to each element
    ss = an.flatten()[idxs_valid_flat]  # flattened subset of data we want 
    I, J, K = np.indices(r2.shape)    
    idxs = np.stack([I.flatten(), J.flatten(), K.flatten()], axis=1)
    idxs = idxs[idxs_valid_flat]        # flattened cartesian coordinates 

    # For each possible orientation, we determine how well angle pRF
    # params can be explained by linear regression between angle (in 0 to 180 degrees)
    # and the projection of the voxel coordinates onto that oriented vector
    rsqr = np.zeros([len(azimuths), len(elevations)])
    for i in range(len(azimuths)):
        for j in range(len(elevations)):
            # Create vector from cartesian coordinates
            v = np.array([v_X[i, j], v_Y[i, j], v_Z[i, j]])

            # Determine projections of coordinates
            v_mag = np.linalg.norm(v)
            scalar_projections = np.dot(idxs, v) / v_mag

            # Correlate projections with true angles
            _, _, r_value, _, _= scipy.stats.linregress(scalar_projections, ss)
            rsqr[i, j] = r_value

    # Check if we need to do additional work 
    if n_rep_shuffle > 1:
        # For n_rep_shuffle > 1, we are creating a control map by shuffling the data,
        # performing the analysis again, and then rotating the shuffled map result to
        # peak-match the true result

        # First, we need to identify the best fitting vector for reference
        idx_best = np.argmax(rsqr)
        v_best = np.array([v_X.flatten()[idx_best], v_Y.flatten()[idx_best], v_Z.flatten()[idx_best]])

        # Now, we need to loop through each shuffle
        rsqr_shuffled = np.zeros([len(azimuths), len(elevations), n_rep_shuffle])
        for k in range(n_rep_shuffle):
            # Shuffle the correspondence between the data and their coordinates
            idxs_shuffle = np.random.permutation(idxs)

            # Do analysis
            for i in range(len(azimuths)):
                for j in range(len(elevations)):
                    # Create vector from cartesian coordinates
                    v = np.array([v_X[i, j], v_Y[i, j], v_Z[i, j]])

                    # Determine projections of coordinates
                    v_mag = np.linalg.norm(v)
                    scalar_projections = np.dot(idxs_shuffle, v) / v_mag

                    # Correlate projections with true angles
                    _, _, r_value, _, _= scipy.stats.linregress(scalar_projections, ss)
                    rsqr_shuffled[i, j, k] = r_value
    
    # Select colormap norm
    norm = matplotlib.colors.Normalize(vmin=-0.45, vmax=0.45)

    # Select coordiantes
    if view == 'coronal':
        camera = np.array([0.0, -90.0, 0.0])
    elif view == 'axial':
        camera = np.array([90.0, -90.0, 0.0])
    else:
        camera = np.array([0.0, 0.0, 0.0])

    # Plot colored surface based on which plot we want
    if n_rep_shuffle >= 1 and n_rep_shuffle <= 20:
        # Create empty figure
        fig = None
        if n_rep_shuffle == 1:
            fig = plt.figure(figsize=(5.0, 5.0), dpi=300)
        else:
            fig = plt.figure(figsize=(5.0, 10.0), dpi=300)
        for i in range(n_rep_shuffle):
            # Create axis and plot into
            ax = None
            if n_rep_shuffle == 1:
                ax = fig.add_subplot(1, 1, 1, projection='3d')
            else:
                ax = fig.add_subplot(5, 4, i+1, projection='3d')
            if n_rep_shuffle == 1:
                ax.plot_surface(v_X, v_Y, v_Z, facecolors=cmap_r(norm(rsqr)), rstride=1, cstride=1, shade=False)
            else:
                ax.plot_surface(v_X, v_Y, v_Z, facecolors=cmap_r(norm(np.squeeze(rsqr_shuffled[:, :, i]))), rstride=1, cstride=1, shade=False)

            # Set aspect ratio and view
            ax.set_aspect('equal')
            ax.view_init(camera[0], camera[1], camera[2])

            # Set ticks, lims, axes, etc.
            if n_rep_shuffle == 1:
                ax.set_title(view, fontsize=4.0)
            ax.set_xticks([])
            ax.set_yticks([])
            ax.set_zticks([])
            ax.set_axis_off()
            ax.set_xlim([-0.9, 0.9])
            ax.set_ylim([-0.9, 0.9])
            ax.set_zlim([-0.9, 0.9])
            ax.set_box_aspect([1,1,1])
            if view == 'axial':
                if side == 'left':
                    ax.set_xlabel('L-M', fontsize=20.0)
                else:
                    ax.set_xlabel('M-L', fontsize=20.0)
                ax.set_ylabel('P-A', fontsize=4.0)
                ax.plot([-1.0, 1.0], [0.0, 0.0], [1.0, 1.0], color='k', zorder=10)
                ax.plot([0.0, 0.0], [-1.0, 1.0], [1.0, 1.0], color='k', zorder=10)
            elif view == 'coronal':
                if side == 'left':
                    ax.set_xlabel('L-M', fontsize=4.0)
                else:
                    ax.set_xlabel('M-L', fontsize=4.0)
                ax.set_zlabel('I-S', fontsize=4.0)
                ax.plot([-1.0, 1.0], [-1.0, -1.0], [0.0, 0.0], color='k', zorder=10)
                ax.plot([0.0, 0.0], [-1.0, -1.0], [-1.0, 1.0], color='k', zorder=10)
            else:
                ax.set_ylabel('P-A', fontsize=4.0)
                ax.set_zlabel('I-S', fontsize=4.0)
                ax.plot([1.0, 1.0], [-1.0, 1.0], [0.0, 0.0], color='k', zorder=10)
                ax.plot([1.0, 1.0], [0.0, 0.0], [-1.0, 1.0], color='k', zorder=10)
            fig.tight_layout(pad=0.0)
    else:
        # Create empty figure
        fig = plt.figure(figsize=(10.0, 5.0), dpi=300)
        ax = fig.add_subplot()
        ax.hist(np.max(np.max(rsqr_shuffled, 0), 0))
        ax.vlines([np.max(rsqr.flatten())], 0.0, 5.0, color='red')
        plt.xlim([0.0, np.max(rsqr.flatten()) * 1.2])
        plt.ylabel('Bin frequency', fontsize=45.0)
        plt.xlabel('Corr. coef.', fontsize=45.0)
        plt.xticks(np.arange(0.0, 0.5, 0.1), fontsize=35.0)
        ax.spines[['right', 'top']].set_visible(False)
        fig.tight_layout(pad=0.0)

    # Print map properties
    idx_best = np.argmax(rsqr)
    print("Best fit: r = " + str(round(rsqr.flatten()[idx_best], ndigits=2)))
    print("Best fit: v = " + str([v_X.flatten()[idx_best], v_Y.flatten()[idx_best], v_Z.flatten()[idx_best]]))
    print("Best fit: N = " + str(len(ss)))


def plot_empty_polar_axis():
    # Create figure
    fig = plt.figure(figsize=(3, 3), dpi=300)
    plt.subplot(1, 1, 1, projection='polar')

    # Set limits and ticks and labels
    plt.xlim(0.0, 2*np.pi)
    plt.yticks(np.arange(0.0, 0.5, 0.1), [])
    plt.xticks([0, np.pi/2, np.pi, 2*np.pi/3], [])
    plt.ylim(0.0, 0.5)
    plt.gca().set_thetamin(0)
    plt.gca().set_thetamax(360.0)
    plt.gca().set_rmax(0.4)
    plt.gca().set_rmin(0) 
    plt.gca().set_rorigin(0) 
    plt.gca().set_rticks(np.arange(0.0, 0.5, 0.1))
    plt.gca().set_thetagrids(np.linspace(0.0, 270.0, 4))
    fig.tight_layout()
    plt.savefig(os.path.join('../figures', 'polar_axis_empty.png'), bbox_inches='tight', pad_inches=0.1)


def char_prf_eccentricity(ec, r2, thal, slices=None, n_rep_shuffle=5):
    """ Characterizes the organization of retinotopic maps in terms of their ecc progression

    Analyzes size data in `ec` to characterize how well it obeys a Gaussian gradient w.r.t
    the center of mass of the (inverse) pRF size data within the mask.
    """
    # Note:
        # indices [0:93, :, :] -> LH
        # indices [94:, :, :] -> RH

    # Create empty figure
    fig = plt.figure(figsize=(12, 2.5), dpi=300)

    # Handle slices
    if slices is None:
        slices = np.arange(102, 93, -1)
    
    # Determine matrix indiceds that are valid (needs to be left hemisphere and variance
    # explained greater than 0.1% to be considered for testing)
    hemi = np.ones(ec.shape, dtype=bool)
    hemi[94:, :, :] = False
    idxs_valid = np.logical_and(r2 > 0.1, thal == 1)
    idxs_valid = np.logical_and(idxs_valid, hemi)

    # Determine sigma values we want to test for size variation with isometric Gaussian in voxel size units
    sigmas = np.exp(np.linspace(np.log(1.0), np.log(20.0), 45))

    # Create matrix to store result
    results = np.zeros([len(sigmas), len(slices)])
    results_shuffled = np.zeros([len(sigmas), len(slices), n_rep_shuffle])

    # Loop through slices, for each, we want to correlate size with each voxel's distance
    # from the center of mass expressed in Guassian sigma units
    for i, slice in enumerate(slices):
        # First, for this slice, determine which of its flattened indices are valid
        idxs_valid_flat = np.squeeze(idxs_valid[:, slice, :]).flatten()

        # Next, squeeze the data in the slice and subset to include only valid elements. We
        # also determine and subset the cartesian coordintes of each resulting element
        ss = np.squeeze(ec[:, slice, :])  # now, dim1 -> left-right, dim 2 -> inferior-superior
        ss_r2 = np.squeeze(r2[:, slice, :])
        I, J = np.indices(ss.shape)    
        idxs = np.stack([I.flatten(), J.flatten()], axis=1)
        idxs = idxs[idxs_valid_flat]
        ss = ss.flatten()[idxs_valid_flat]
        ss_r2 = ss_r2.flatten()[idxs_valid_flat]
    
        # Determine the center of mass of the voxels under investigation?
        com = np.sum(np.stack([ss_r2, ss_r2], axis=1) * idxs, 0)/sum(ss_r2)

        # Now, for each possible sigma, we determine how well size pRF parameters can be 
        # predicted by a linear regression betwen size pRF parameter and the voxel's location
        # expressed in terms of distance from the center of the Gaussian distribution.
        curve = np.zeros(len(sigmas))
        for j, sigma in enumerate(sigmas):
            # Express all coordinate values relative to their position w.r.t. Guassian centered
            # at center of mass
            norm = 1/(2*np.pi*sigma**2)
            inner = ((idxs[:, 0] - com[0])/sigma)**2  + ((idxs[:, 1] - com[1])/sigma)**2
            distances = norm * np.exp(-1/2 * inner)

            # Correlate projections with true angles
            _, _, r_value, _, _= scipy.stats.linregress(distances, ss)
            # curve[j] = np.arctan(r_value)/(1/np.sqrt(len(ss) - 3))
            curve[j] = r_value

        results[:, i] = curve

        # Now, we repeat the above process on shuffled copies of the data
        for idx_rep in range(n_rep_shuffle):
            idxs_shuffle = random.sample(range(len(ss)), len(ss))
            ss_shuffled = ss[idxs_shuffle]
            curve = np.zeros(len(sigmas))
            for j, sigma in enumerate(sigmas):
                # Express all coordinate values relative to their position w.r.t. Guassian centered
                # at center of mass
                norm = 1/(2*np.pi*sigma**2)
                inner = ((idxs[:, 0] - com[0])/sigma)**2  + ((idxs[:, 1] - com[1])/sigma)**2
                distances = norm * np.exp(-1/2 * inner)

                # Correlate projections with true angles
                _, _, r_value, _, _= scipy.stats.linregress(distances, ss_shuffled)
                # curve[j] = np.arctan(r_value)/(1/np.sqrt(len(ss) - 3))
                curve[j] = r_value
            
            results_shuffled[:, i, idx_rep] = curve

    # Flip results polarity; if Gaussian value goes down we expect eccentricity to increase
    # (anti-correlated) but we want to test for positive correlation
    results = -1 * results
    results_shuffled = -1 * results_shuffled

    # Now, we loop through slices and plot 
    for idx_slice in range(len(slices)):
        # Create subplot for this plot
        plt.subplot(1, len(slices), idx_slice+1)

        # Plot horizontal line at zero
        plt.hlines(0.0, 0.0, np.max(sigmas), color='black', linewidth=0.5)

        # Test for signif
        Z_real = np.max(np.squeeze(results[:, idx_slice]))
        Z_shuffle = np.max(np.squeeze(results_shuffled[:, idx_slice, :]), 0)
        if Z_real > np.percentile(Z_shuffle, 95.0):
            signif = True
        else:
            signif = False

        # Plot data
        if signif:
            color = 'red'
        else:
            color = 'black'
        plt.semilogx(sigmas, results[:, idx_slice], color=color)

        # Plot shuffled replicates
        plt.plot(sigmas, np.squeeze(results_shuffled[:, idx_slice, :]), color='gray', alpha=0.05, linewidth=0.2)
        
        # Set limits and ticks and labels
        plt.xlim(np.min(sigmas), np.max(sigmas))
        if idx_slice > 0:
            plt.yticks(np.arange(-1.0, 1.5, 0.5), [])
        else:
            plt.yticks(np.arange(-1.0, 1.5, 0.5))
            
        plt.xticks([1, 10], [])
        # plt.ylim(-10, 10)
        plt.ylim(-1.0, 1.0)

        # Add title
        plt.title(str(slices[idx_slice]))

    fig.tight_layout()
    plt.subplots_adjust(wspace=0, hspace=0)

    # Last step: print R2 max and corresponding eccentricity for each slice
    for idx_slice in range(len(slices)):
        r2_max = np.max(results[:, idx_slice])
        sigma_max = sigmas[np.argmax(results[:, idx_slice])]
        print("Slice " + str(slices[idx_slice]) + ": r = " + str(round(r2_max, 3)) + " at " + str(round(sigma_max, 2)))


def plot_prf_parameter_sequence(t1, r2, an, an_orig, ec, sz, thal, roi, slices=None):
    """ Plots variance explained, angle maps, eccentricity maps, and laterality maps on T1-weighted anatomy

    Args:
        t1 (ndarray): T1-weighted anatomy, one subject or average across subjects, full volume
        r2 (ndarray): pRF variance explained, one subject or average across subjects, full volume
        an (ndarray): pRF angle maps, in degrees from 0 to 180, one subject or average across subjects, full volume
        an (ndarray): pRF angle maps, in degrees from 0 to 360, ...
        ec (ndarray): pRF eccentricity maps, one subject or average across subjects, full volume
        sz (ndarray): pRF size maps, one subject or average across subjects, full volume
        thal (ndarray): ROI label indicating the extent of the thalamus
        roi (ndarray): ROI labels indicating different subcortical structures
        slices (list): list of slices to plot
    """
    # Construct figure
    if slices is None:
        slices = [100, 97, 95]
    n_slice = len(slices)
    cutoff = 0.1
    plt.figure(figsize=(10.5, 12))

    lat = calc_lat(an_orig, sz, ec)
    for idx_slice, slice in enumerate(slices):
        # R2
        plt.subplot(4, n_slice, idx_slice + 1)
        plot_slice_with_overlay(t1, r2, slice, 'coronal', np.logical_or(r2 < cutoff, thal != 1), (0 + 50, 182 - 50),
                                (0 + 55, 182 - 80), cmap2=cmap_rsqr, clim2=(0, 1))
        plot_roi_overlay(roi, 'coronal', slice)
        # Angle
        plt.subplot(4, n_slice, idx_slice + 1 + n_slice * 1)
        plot_slice_with_overlay(t1, an, slice, 'coronal', np.logical_or(r2 < cutoff, thal != 1), (0 + 50, 182 - 50),
                                (0 + 55, 182 - 80), cmap2=cmap_angle, clim2=(0, 180))
        # Eccentricity
        plt.subplot(4, n_slice, idx_slice + 1 + n_slice * 2)
        plot_slice_with_overlay(t1, ec, slice, 'coronal', np.logical_or(r2 < cutoff, thal != 1), (0 + 50, 182 - 50),
                                (0 + 55, 182 - 80), cmap2=cmap_size, clim2=(0, 6))
        # Laterality
        plt.subplot(4, n_slice, idx_slice + 1 + n_slice * 3)
        plot_slice_with_overlay(t1, lat, slice, 'coronal', np.logical_or(r2 < cutoff, thal != 1), (0 + 50, 182 - 50),
                                (0 + 55, 182 - 80), cmap2=cmap_lat, clim2=(-1, 1))
        # RF size
        #plt.subplot(4, n_slice, idx_slice + 1 + n_slice * 3)
        #plot_slice_with_overlay(t1, sz, slice, 'coronal', np.logical_or(r2 < cutoff, thal != 1), (0 + 50, 182 - 50),
        #                        (0 + 55, 182 - 80), cmap2=cmap_size, clim2=(0, 6))
    plt.tight_layout(pad=0.0, h_pad=0.1, w_pad=0.1)


def plotsave_prf_parameter_sequence_reduced(t1, r2, an, an_orig, ec, sz, thal, roi, savename, slices=None, outline_width=1.0):
    """ Plots variance explained, angle maps, and laterality maps on T1-weighted anatomy

    Args:
        t1 (ndarray): T1-weighted anatomy, one subject or average across subjects, full volume
        r2 (ndarray): pRF variance explained, one subject or average across subjects, full volume
        an (ndarray): pRF angle maps, in degrees from 0 to 180, one subject or average across subjects, full volume
        an (ndarray): pRF angle maps, in degrees from 0 to 360, ...
        ec (ndarray): pRF eccentricity maps, one subject or average across subjects, full volume
        sz (ndarray): pRF size maps, one subject or average across subjects, full volume
        thal (ndarray): ROI label indicating the extent of the thalamus
        roi (ndarray): ROI labels indicating different subcortical structures
        slices (list): list of slices to plot
    """
    # Construct figure
    if slices is None:
        slices = [100, 97, 95]
    n_slice = len(slices)
    cutoff = 0.1
    lat = calc_lat(an_orig, sz, ec)
    # Plot R2
    plt.figure(figsize=(7, 12))
    for idx_slice, slice in enumerate(slices):
        # R2
        plt.subplot(n_slice, 1, idx_slice + 1)
        plot_slice_with_overlay(t1, r2, slice, 'coronal', np.logical_or(r2 < cutoff, thal != 1), (0 + 50, 182 - 50),
                                (0 + 55, 182 - 80), cmap2=cmap_rsqr, clim2=(0, 1))
        plot_roi_overlay(roi, 'coronal', slice, outline_width)
    plt.tight_layout()
    plt.savefig(savename + 'R2.png')
    # Plot Angle
    plt.figure(figsize=(7, 12))
    for idx_slice, slice in enumerate(slices):
        plt.subplot(n_slice, 1, idx_slice + 1)
        plot_slice_with_overlay(t1, an, slice, 'coronal', np.logical_or(r2 < cutoff, thal != 1), (0 + 50, 182 - 50),
                                (0 + 55, 182 - 80), cmap2=cmap_angle, clim2=(0, 180))
    plt.tight_layout()
    plt.savefig(savename + 'angle.png')
    # Plot laterality
    plt.figure(figsize=(7, 12))
    for idx_slice, slice in enumerate(slices):
        plt.subplot(n_slice, 1, idx_slice + 1)
        plot_slice_with_overlay(t1, lat, slice, 'coronal', np.logical_or(r2 < cutoff, thal != 1), (0 + 50, 182 - 50),
                                (0 + 55, 182 - 80), cmap2=cmap_lat, clim2=(-1, 1))
    plt.tight_layout()
    plt.savefig(savename + 'laterality.png')
    # Plot eccentricity
    plt.figure(figsize=(7, 12))
    for idx_slice, slice in enumerate(slices):
        plt.subplot(n_slice, 1, idx_slice + 1)
        plot_slice_with_overlay(t1, ec, slice, 'coronal', np.logical_or(r2 < cutoff, thal != 1), (0 + 50, 182 - 50),
                                (0 + 55, 182 - 80), cmap2=cmap_size, clim2=(0, 6))
    plt.tight_layout()
    plt.savefig(savename + 'eccentricity.png')



def plotsave_prf_parameter_sequence_reduced_arcaro_custom(t1, r2, an, an_orig, ec, sz, thal, roi, savename, slices=None, outline_width=1.0, cutoff=0.1):
    """ Plots variance explained, angle maps, and laterality maps on T1-weighted anatomy

    Args:
        t1 (ndarray): T1-weighted anatomy, one subject or average across subjects, full volume
        r2 (ndarray): pRF variance explained, one subject or average across subjects, full volume
        an (ndarray): pRF angle maps, in degrees from 0 to 180, one subject or average across subjects, full volume
        an (ndarray): pRF angle maps, in degrees from 0 to 360, ...
        ec (ndarray): pRF eccentricity maps, one subject or average across subjects, full volume
        sz (ndarray): pRF size maps, one subject or average across subjects, full volume
        thal (ndarray): ROI label indicating the extent of the thalamus
        roi (ndarray): ROI labels indicating different subcortical structures
        slices (list): list of slices to plot
    """
    # Construct figure
    if slices is None:
        slices = [100, 97, 95]
    n_slice = len(slices)
    lat = calc_lat(an_orig, sz, ec)
    # Plot R2
    plt.figure(figsize=(7, 12))
    for idx_slice, slice in enumerate(slices):
        # R2
        plt.subplot(n_slice, 1, idx_slice + 1)
        plot_slice_with_overlay(t1, r2, slice, 'coronal', np.logical_or(r2 < cutoff, thal != 1), (0 + 50, 182 - 50),
                                (0 + 55, 182 - 80), cmap2=cmap_rsqr, clim2=(0, 1))
        plot_roi_overlay(roi, 'coronal', slice, outline_width)
    plt.tight_layout()
    plt.savefig(savename + 'R2.png')


def plot_mip_roi_overlay(ax, ROI, view, color=None, linestyle='solid', linewidth=2, selected_rois=None):
    # Handle which ROIs to plot
    if selected_rois is None:
        selected_rois = [1, 2, 3, 4, 5]
    # Extract and overlay ROIs
    for roi_val in selected_rois:
        # Look only for one ROI val
        roi_present = np.array(ROI) == roi_val
        # Average over subjects
        roi_present = np.mean(roi_present, axis=0) > 2/8
        # Branch based on requested view
        if view == 'coronal':
            roi_present = np.max(roi_present, axis=1)
        elif view == 'sagittal':
            roi_present = np.max(roi_present, axis=0)
        elif view == 'axial':
            roi_present = np.max(roi_present, axis=2)
        contour = measure.find_contours(roi_present, 0.8)
        if len(contour) > 0:
            for cont in contour:
                if color is None:
                    ax.plot(cont[:, 0], cont[:, 1], linewidth=linewidth, linestyle=linestyle,
                             color=cmap_thal(np.linspace(0, 1, 5))[roi_val-1])
                elif isinstance(color, list):
                    ax.plot(cont[:, 0], cont[:, 1], linewidth=linewidth, linestyle=linestyle, color=color[roi_val-1])
                else:
                    ax.plot(cont[:, 0], cont[:, 1], linewidth=linewidth, linestyle=linestyle, color=color)


def plot_roi_overlay(ROI, view, idx, outline_width=1.0):
    # Extract and overlay ROIs
    for roi_val in [1, 2, 3, 4, 5]:
        # Look only for one ROI val
        roi_present = np.array(ROI) == roi_val
        # Average over subjects
        roi_present = np.mean(roi_present, axis=0) > 2/8
        # Branch based on requested view
        if view == 'coronal':
            roi_present = roi_present[:, idx, :]
        elif view == 'sagittal':
            roi_present = roi_present[idx, :, :]
        elif view == 'axial':
            roi_present = roi_present[:, :, idx]
        contour = measure.find_contours(roi_present, 0.8)
        if len(contour) > 0:
            for cont in contour:
                 plt.plot(cont[:, 0], cont[:, 1], linewidth=outline_width,
                          color=cmap_thal(np.linspace(0, 1, 5))[roi_val-1])


def plot_slice(img1, idx, view, lims_x=(0, 200), lims_y=(0, 200), cmap1='gray', clim1=None, marker=None,
               marker_color='black', marker_angle=0, marker_length=8):
    """ Plots slice of 3D brick of data.

    Args:
        img1 (ndarray): Volume of data, 3D
        idx (int): Index of the slice to plot
        view (str): either 'sagittal', 'coronal', or 'axial'
        lims_x (tuple): Limits in the plotted left-right axis
        lims_y (tuple): Limits in the plotted up-down axis
        cmap1 (str or np.array): Colormap to use for img1
        marker (tuple, list): a tuple or list of coordinates, an arrow is placed at these points
        marker_color (ndarray, string): some matplotlib-acceptable representation of a color to make the plotted arrow
        marker_angle (float): the angle of the marker in radians
        marker_length (float): the length of the marker from tip to tail
    """
    # Calculate marker start and stop positions
    x = marker_length * np.cos(marker_angle)
    y = marker_length * np.sin(marker_angle)
    # Plot marker
    if marker is not None:
        plt.arrow(marker[0]+x, marker[1]+y, -x, -y, color=marker_color, width=0.75, head_width=2, head_length=2,
                  length_includes_head=True)
    # Plot slice
    if view == 'sagittal':
        plt.imshow(np.squeeze(img1[idx, :, :]).T, cmap=cmap1)
    elif view == 'coronal':
        plt.imshow(np.squeeze(img1[:, idx, :]).T, cmap=cmap1)
    elif view == 'axial':
        plt.imshow(np.squeeze(img1[:, :, idx]).T, cmap=cmap1)
    if clim1 is not None:
        plt.clim(clim1)
    plt.xlim(lims_x)
    plt.ylim(lims_y)
    plt.gca().get_xaxis().set_visible(False)
    plt.gca().get_yaxis().set_visible(False)


def plot_slice_oo(ax, img1, idx, view, lims_x=(0, 200), lims_y=(0, 200), cmap1='gray', clim1=None, marker=None,
               marker_color='black', marker_angle=0, marker_length=8):
    """ Plots slice of 3D brick of data.

    Args:
        img1 (ndarray): Volume of data, 3D
        idx (int): Index of the slice to plot
        view (str): either 'sagittal', 'coronal', or 'axial'
        lims_x (tuple): Limits in the plotted left-right axis
        lims_y (tuple): Limits in the plotted up-down axis
        cmap1 (str or np.array): Colormap to use for img1
        marker (tuple, list): a tuple or list of coordinates, an arrow is placed at these points
        marker_color (ndarray, string): some matplotlib-acceptable representation of a color to make the plotted arrow
        marker_angle (float): the angle of the marker in radians
        marker_length (float): the length of the marker from tip to tail
    """
    # Calculate marker start and stop positions
    x = marker_length * np.cos(marker_angle)
    y = marker_length * np.sin(marker_angle)
    # Plot marker
    if marker is not None:
        ax.arrow(marker[0]+x, marker[1]+y, -x, -y, color=marker_color, width=0.75, head_width=2, head_length=2,
                  length_includes_head=True)
    # Plot slice
    if view == 'sagittal':
        ax.imshow(np.squeeze(img1[idx, :, :]).T, cmap=cmap1)
    elif view == 'coronal':
        ax.imshow(np.squeeze(img1[:, idx, :]).T, cmap=cmap1)
    elif view == 'axial':
        ax.imshow(np.squeeze(img1[:, :, idx]).T, cmap=cmap1)
    if clim1 is not None:
        ax.set_clim(clim1)
    ax.set_xlim(lims_x)
    ax.set_ylim(lims_y)
    ax.get_xaxis().set_visible(False)
    ax.get_yaxis().set_visible(False)


def plot_slice_token_markers(img1, idx, view, lims_x=(0, 200), lims_y=(0, 200), cmap1='gray', clim1=None, marker=None,
                             marker_color='black', marker_size=5):
    """ Plots slice of 3D brick of data.

    Args:
        img1 (ndarray): Volume of data, 3D
        idx (int): Index of the slice to plot
        view (str): either 'sagittal', 'coronal', or 'axial'
        lims_x (tuple): Limits in the plotted left-right axis
        lims_y (tuple): Limits in the plotted up-down axis
        cmap1 (str or np.array): Colormap to use for img1
        marker (tuple, list): a tuple or list of coordinates, an arrow is placed at these points
        marker_color (ndarray, string): some matplotlib-acceptable representation of a color to make the plotted arrow
        marker_angle (float): the angle of the marker in radians
        marker_length (float): the length of the marker from tip to tail
    """
    # Plot marker
    if marker is not None:
        if isinstance(marker, list):
            for _m, _c in zip(marker, marker_color):
                plt.plot(_m[0], _m[1], color=_c, marker='o', markersize=5)
        else:
            plt.plot(marker[0], marker[1], color=marker_color, marker='o', markersize=5)
    # Plot slice
    if view == 'sagittal':
        plt.imshow(np.squeeze(img1[idx, :, :]).T, cmap=cmap1)
    elif view == 'coronal':
        plt.imshow(np.squeeze(img1[:, idx, :]).T, cmap=cmap1)
    elif view == 'axial':
        plt.imshow(np.squeeze(img1[:, :, idx]).T, cmap=cmap1)
    if clim1 is not None:
        plt.clim(clim1)
    plt.xlim(lims_x)
    plt.ylim(lims_y)
    plt.gca().get_xaxis().set_visible(False)
    plt.gca().get_yaxis().set_visible(False)


def plot_slice_with_overlay(img1, img2, idx, view, mask, lims_x=(0, 200), lims_y=(0, 200), cmap1='gray', cmap2=cmap_rsqr,
                            clim1=None, clim2=None, marker=None, marker_color='black', marker_angle=0, marker_length=8):
    """ Plots slice of 3D brick of data with a corresponding slice of an overlay.

    Args:
        img1 (ndarray): Volume of data, 3D
        img2 (ndarray): Volume of data, 3D, plotted as overlay
        idx (int): Index of the slice to plot
        view (str): either 'sagittal', 'coronal', or 'axial'
        lims_x (tuple): Limits in the plotted left-right axis
        lims_y (tuple): Limits in the plotted up-down axis
        cmap1 (str or np.array): Colormap to use for img1
        cmap2 (str or np.array): Colormap to use for img2
        clim1 (tuple, list): limits of the first colormap. If None, default limits are applied.
        clim2 (tuple, list): limits of the second colormap. If None, default limits are applied.
        marker (tuple, list): a tuple or list of coordinates, an arrow is placed at these points
        marker_color (ndarray, string): some matplotlib-acceptable representation of a color to make the plotted arrow
        marker_angle (float): the angle of the marker in radians
        marker_length (float): the length of the marker from tip to tail
    """
    # Mask overlay image where mask indicates
    img2 = np.ma.masked_where(mask, img2)
    # Plot underlay
    plot_slice(img1, idx, view, lims_x, lims_y, cmap1, clim1, marker, marker_color, marker_angle, marker_length)
    # Plot overlay
    plot_slice(img2, idx, view, lims_x, lims_y, cmap2, clim2)


def plot_mip(img1, view, lims_x=(0, 200), lims_y=(0, 200), cmap1='gray', clim1=None, marker=None,
            marker_color='black', marker_angle=0, marker_length=8):
    """ Plots maximum intensity projection of 3D volume of data

    Args:
        img1 (ndarray): Volume of data, 3D
        view (str): either 'sagittal', 'coronal', or 'axial'
        lims_x (tuple): Limits in the plotted left-right axis
        lims_y (tuple): Limits in the plotted up-down axis
        cmap1 (str or np.array): Colormap to use for img1
        marker (tuple, list): a tuple or list of coordinates, an arrow is placed at these points
        marker_color (ndarray, string): some matplotlib-acceptable representation of a color to make the plotted arrow
        marker_angle (float): the angle of the marker in radians
        marker_length (float): the length of the marker from tip to tail
    """
    # Calculate marker start and stop positions
    x = marker_length * np.cos(marker_angle)
    y = marker_length * np.sin(marker_angle)
    # Plot marker
    if marker is not None:
        plt.arrow(marker[0]+x, marker[1]+y, -x, -y, color=marker_color, width=0.75, head_width=2, head_length=2,
                  length_includes_head=True)
    # Plot slice
    if view == 'sagittal':
        plt.imshow(np.squeeze(np.nanmax(img1, axis=0)).T, cmap=cmap1)
    elif view == 'coronal':
        plt.imshow(np.squeeze(np.nanmax(img1, axis=1)).T, cmap=cmap1)
    elif view == 'axial':
        plt.imshow(np.squeeze(np.nanmax(img1, axis=2)).T, cmap=cmap1)
    if clim1 is not None:
        plt.clim(clim1)
    plt.xlim(lims_x)
    plt.ylim(lims_y)
    plt.gca().get_xaxis().set_visible(False)
    plt.gca().get_yaxis().set_visible(False)


def plot_rainbow_slice(img1, img2, idx_slice, mask, lims_x=(50, 132), lims_y=(40, 102)):
    '''
    Plots base layer of data and overlay layer of data with 3D color
    Args:
        img1 (np.array): Volume of data, plotted as base layer
        img2 (np.array): Volume of data, plotted as overlay
        idx_slice (int): Index of the coronal slice to plot
        mask (np.array): Mask of which voxels to *exclude* in overlay
        lims_x (tuple): Limits in the left-right axis
        lims_y (tuple): Limits in the anterior-posterior axis
        cmap1 (str or np.array): Colormap to use for img1
        cmap2 (str or np.array): Colormap to use for img2
        transp_mask (np.array): Mask of voxels to visually suppress in img2 data via gray transparency overlay
    '''
    img2 = np.concatenate([img2, np.expand_dims(255*np.logical_not(mask), axis=3)], axis=3)
    #mask = np.stack([mask, mask, mask], axis=3)  # stack mask to match the fact that the data has color dimension
    #img2 = np.ma.masked_where(mask, img2)
    plt.imshow(np.squeeze(img1[:, idx_slice, :]).T, cmap='gray')
    plt.imshow(np.transpose(np.squeeze(img2[:, idx_slice, :, :]), [1, 0, 2]))
    plt.axis('off')
    plt.xlim(lims_x)
    plt.ylim(lims_y)


def plot_visual_field_coverage(R2, AN, EC, SZ, include, cutoff=0.25):
    '''
    Plots visual fields in stimulus space
    Args:
        R2 (list): list of 3D array of size (182, 218, 182) with variance explained metric
        AN (list): list of 3D array of size (182, 218, 182) with angle parameter estimates
        EC (list): list of 3D array of size (182, 218, 182) with eccentricity parameter estimates
        SZ (list): list of 3D array of size (182, 218, 182) with size parameter estimates
        include (ndarray): 3D array of size (182, 218, 182) with boolean values indicating which voxels to plot
        cutoff (float): cutoff of which voxels to include based on R2
    Returns:
        None

    '''
    fig, ax = plt.subplots(8, 2, figsize=(3.75, 16))
    for subj in range(8):
        # Prepare variables to calculate coverage maps
        maps = []
        x = np.linspace(-8, 8, 200)
        y = np.linspace(-8, 8, 200)
        # Design plots
        for row in range(2):
            # Lay out basic markings on plot
            ax[subj][row].plot([-10, 10], [0, 0], 'r--')
            ax[subj][row].plot([0, 0], [-10, 10], 'r--')
            rect = plt.Rectangle(xy=(-4.2, -4.2), width=8.4, height=8.4, color='r', fill=False)
            ax[subj][row].add_artist(rect)
            # Set limits
            ax[subj][row].set_xlim((-8.4, 8.4))
            ax[subj][row].set_ylim((-8.4, 8.4))
            # Set ticks
            ax[subj][row].set_xticks([-8.4, -4.2, 0, 4.2, 8.4])
            ax[subj][row].set_xticks([-6.3, -2.1, 2.1, 6.3], minor=True)
            ax[subj][row].set_yticks([-8.4, -4.2, 0, 4.2, 8.4])
            ax[subj][row].set_yticks([-6.3, -2.1, 2.1, 6.3], minor=True)
            # Force square axes
            ax[subj][row].set(adjustable='box', aspect='equal')
        # Remove unnecessary x/y-axes
        if subj != 7:
            ax[subj][0].get_xaxis().set_visible(False)
            ax[subj][1].get_xaxis().set_visible(False)
        ax[subj][1].get_yaxis().set_visible(False)
        for hemi in range(2):
            # Set all R2 to zero where voxels are not in include
            R2[subj][include != 1] = 0
            # Loop through hemispheres and plot RFs
            if hemi == 0:  # left hemisphere
                R2_subset = R2[subj][0:93, :, :]
                AN_subset = AN[subj][0:93, :, :]
                EC_subset = EC[subj][0:93, :, :]
                SZ_subset = SZ[subj][0:93, :, :]
            else:          # right hemisphere
                R2_subset = R2[subj][93:, :, :]
                AN_subset = AN[subj][93:, :, :]
                EC_subset = EC[subj][93:, :, :]
                SZ_subset = SZ[subj][93:, :, :]
            n_vox = len(R2_subset[R2_subset > cutoff])
            # Plot fields
            for vox in range(n_vox):
                # Extract parameters
                ang = AN_subset[R2_subset > cutoff][vox]
                ecc = EC_subset[R2_subset > cutoff][vox]
                size = SZ_subset[R2_subset > cutoff][vox]
                # Plot center
                if hemi == 0:
                    ax[subj][0].plot(ecc*np.cos(ang*np.pi/180), ecc*np.sin(ang*np.pi/180), color=[0, 1, 0], markersize=3, marker='o')
                else:
                    ax[subj][0].plot(ecc*np.cos(ang*np.pi/180), ecc*np.sin(ang*np.pi/180), color=[0, 0, 1], markersize=3, marker='o')
                # Plot surround
                circ = plt.Circle((ecc*np.cos(ang*np.pi/180), ecc*np.sin(ang*np.pi/180)), size, color='k', fill=False, linewidth=0.5, alpha=0.1)
                ax[subj][0].add_artist(circ)
            # Calculate coverage maps
            X, Y = np.meshgrid(x, y)
            z = np.zeros((200, 200))
            for vox in range(n_vox):
                # Extract parameters
                ang = AN_subset[R2_subset > cutoff][vox]
                ecc = EC_subset[R2_subset > cutoff][vox]
                size = SZ_subset[R2_subset > cutoff][vox]
                in_circle = np.sqrt(
                     (X - ecc * np.cos(ang * np.pi / 180)) ** 2 + (Y - ecc * np.sin(ang * np.pi / 180)) ** 2) < size
                z = z + in_circle
            maps.append(z)
        # Combine coverage maps
        coverage_map_oppositional = -1*maps[0] + maps[1]
        rgb_lh = np.zeros((200, 200, 3))
        rgb_lh[:, :, 1] = maps[0] / max(np.max(maps[0]), np.max(maps[1]))
        rgb_rh = np.zeros((200, 200, 3))
        rgb_rh[:, :, 2] = maps[1] / max(np.max(maps[0]), np.max(maps[1]))
        coverage_map = rgb_lh + rgb_rh
        # Plot
        #ax[subj][1].pcolormesh(x, y, coverage_map/np.max(np.abs(coverage_map)), cmap='GnBu', clim=(-1, 1))
        ax[subj][1].imshow(coverage_map, extent=(-8.4, 8.4, -8.4, 8.4), origin='lower')
        #ax[subj][1].set_ylim((-8.4, 8.4))  # todo: investigate
        # Calculate LH peak
        #matcoord_peak_LH = np.unravel_index(np.argmax(-1*coverage_map), coverage_map.shape)
        #matcoord_peak_RH = np.unravel_index(np.argmax(coverage_map), coverage_map.shape)
        #ax[subj][1].plot(X[matcoord_peak_LH], Y[matcoord_peak_LH], marker='x', color='r')
        #ax[subj][1].plot(X[matcoord_peak_RH], Y[matcoord_peak_RH], marker='x', color='r')
    plt.tight_layout(w_pad=0)


def summarize_centrality_rf_coverage(R2, EC, include, cutoff=0.1):
    '''
    Companion function to `plot_visual_field_coverage_horizontal` below. Prints to console
    the proportions of each set of eccentricity estimates that fall within inclusion 
    criteria that are within vs without 3 degrees of central visual field.

    Args:
        R2 (list): list of 3D array of size (182, 218, 182) with variance explained metric
        EC (list): list of 3D array of size (182, 218, 182) with eccentricity parameter estimates
        include (array): 3D array of size (182, 218, 182) with boolean values indicating which voxels to plot
        cutoff (float): cutoff of which voxels to include based on R2
    '''
    for subj in range(9):
        # Set all R2 to zero where voxels are not in include
        R2[subj][include != 1] = 0
        # Temporarily reassign R2[subj] and EC[subj]
        R2_temp = R2[subj]
        EC_temp = EC[subj]
        # Identify how many voxels exceed cutoff
        n_vox = len(R2_temp[R2_temp > cutoff])
        # Identify how many voxels exceed cutoff and also are within 3 degrees of center
        n_within = len(R2_temp[np.logical_and(EC_temp <= 3.0, R2_temp > cutoff)])
        n_without = len(R2_temp[np.logical_and(EC_temp > 3.0, R2_temp > cutoff)])
        assert(n_within+n_without == n_vox)
        print("Subj " + str(subj) + ": " + str(n_within/n_vox*100) + "%")


def plot_visual_field_coverage_horizontal(R2, AN, EC, SZ, include, cutoff=0.25, max_vox=100):
    '''
    Plots visual fields in stimulus space
    Args:
        AN (list): list of 3D array of size (182, 218, 182) with angle parameter estimates
        EC (list): list of 3D array of size (182, 218, 182) with eccentricity parameter estimates
        SZ (list): list of 3D array of size (182, 218, 182) with size parameter estimates
        R2 (list): list of 3D array of size (182, 218, 182) with variance explained metric
        include (array): 3D array of size (182, 218, 182) with boolean values indicating which voxels to plot
        cutoff (float): cutoff of which voxels to include based on R2
    Returns:
        None

    '''
    fig, ax = plt.subplots(2, 9, figsize=(18, 4.875))  # make a figure with two rows (different maps) and 9 columns (different subjects and mean subject)
    for subj in range(9):
        # Prepare variables to calculate coverage maps
        maps = []
        x = np.linspace(-8, 8, 200)
        y = np.linspace(-8, 8, 200)
        # Design plots
        for row in range(2):
            # Lay out basic markings on plot
            ax[row][subj].plot([-10, 10], [0, 0], linestyle='dashed', color='gray')
            ax[row][subj].plot([0, 0], [-10, 10], linestyle='dashed', color='gray')
            rect = plt.Rectangle(xy=(-4.2, -4.2), width=8.4, height=8.4, color='gray', fill=False)
            ax[row][subj].add_artist(rect)
            # Set limits
            ax[row][subj].set_xlim((-8.4, 8.4))
            ax[row][subj].set_ylim((-8.4, 8.4))
            # Set ticks
            #ax[row][subj].set_xticks([-8.4, -4.2, 0, 4.2, 8.4])
            #ax[row][subj].set_xticklabels(['', '-4.2', '0', '4.2', ''])
            #ax[row][subj].set_xticklabels(['', '', '', '', ''])
            ax[row][subj].set_xticks([])
            #ax[row][subj].set_xticks([-6.3, -2.1, 2.1, 6.3], minor=True)
            #ax[row][subj].set_yticks([-8.4, -4.2, 0, 4.2, 8.4])
            #ax[row][subj].set_yticks([-6.3, -2.1, 2.1, 6.3], minor=True)
            #ax[row][subj].set_yticklabels(['', '', '', '', ''])
            ax[row][subj].set_yticks([])
            # Force square axes
            ax[row][subj].set(adjustable='box', aspect='equal')
        # Remove unnecessary x/y-axes
        if subj > 0:
            ax[0][subj].get_yaxis().set_visible(False)
            ax[1][subj].get_yaxis().set_visible(False)
        ax[0][subj].get_xaxis().set_visible(False)
        for hemi in range(2):
            # Set all R2 to zero where voxels are not in include
            R2[subj][include != 1] = 0
            # Loop through and plot RFs
            if hemi == 0:  # left hemisphere
                R2_subset = R2[subj][0:93, :, :]
                AN_subset = AN[subj][0:93, :, :]
                EC_subset = EC[subj][0:93, :, :]
                SZ_subset = SZ[subj][0:93, ::, :]
            else:          # right hemisphere
                R2_subset = R2[subj][93:, :, :]
                AN_subset = AN[subj][93:, :, :]
                EC_subset = EC[subj][93:, :, :]
                SZ_subset = SZ[subj][93:, :, :]
            n_vox = len(R2_subset[R2_subset > cutoff])
            R2vals = R2_subset[R2_subset > cutoff]   # grab vector of R2 values that exceed cutoff
            idxs_sort = np.flip(np.argsort(R2vals))  # grab vector to sort R2 values in descending order
            # Plot fields
            for vox in range(min(n_vox, max_vox)):
                # Extract parameters (note that we index into voxels sorted in descending order of R2)
                ang = AN_subset[R2_subset > cutoff][idxs_sort][vox]
                ecc = EC_subset[R2_subset > cutoff][idxs_sort][vox]
                size = SZ_subset[R2_subset > cutoff][idxs_sort][vox]
                # Plot center & surround
                if hemi == 0:
                    ax[0][subj].plot(ecc*np.cos(ang*np.pi/180), ecc*np.sin(ang*np.pi/180), color=[1, 0, 0], markersize=3, marker='o')
                    circ = plt.Circle((ecc * np.cos(ang * np.pi / 180), ecc * np.sin(ang * np.pi / 180)), size,
                                      color=[1, 0, 0], fill=False, linewidth=0.5, alpha=0.1)
                    ax[0][subj].add_artist(circ)
                else:
                    ax[0][subj].plot(ecc*np.cos(ang*np.pi/180), ecc*np.sin(ang*np.pi/180), color=[0, 0, 1], markersize=3, marker='o')
                    circ = plt.Circle((ecc * np.cos(ang * np.pi / 180), ecc * np.sin(ang * np.pi / 180)), size,
                                      color=[0, 0, 1], fill=False, linewidth=0.5, alpha=0.1)
                    ax[0][subj].add_artist(circ)
            # Calculate coverage maps
            X, Y = np.meshgrid(x, y)
            z = np.zeros((200, 200))
            for vox in range(min(n_vox, max_vox)):
                # Extract parameters (note that we index into voxels sorted in descending order of R2)
                ang = AN_subset[R2_subset > cutoff][idxs_sort][vox]
                ecc = EC_subset[R2_subset > cutoff][idxs_sort][vox]
                size = SZ_subset[R2_subset > cutoff][idxs_sort][vox]
                in_circle = np.sqrt(
                     (X - ecc * np.cos(ang * np.pi / 180)) ** 2 + (Y - ecc * np.sin(ang * np.pi / 180)) ** 2) < size
                z = z + in_circle
            maps.append(z)
        # Combine coverage maps
        coverage_map_oppositional = -1*maps[0] + maps[1]
        rgb_lh = np.zeros((200, 200, 3))
        rgb_lh[:, :, 0] = maps[0] / max(np.max(maps[0]), np.max(maps[1]))
        rgb_rh = np.zeros((200, 200, 3))
        rgb_rh[:, :, 2] = maps[1] / max(np.max(maps[0]), np.max(maps[1]))
        coverage_map = rgb_lh + rgb_rh
        # Plot
        #ax[1][subj].pcolormesh(x, y, coverage_map/np.max(np.abs(coverage_map)), cmap='GnBu', clim=(-1, 1))
        ax[1][subj].imshow(coverage_map, extent=(-8.4, 8.4, -8.4, 8.4), origin='lower')
        #ax[subj][1].set_ylim((-8.4, 8.4))  # todo: investigate
        # Calculate LH peak
        #matcoord_peak_LH = np.unravel_index(np.argmax(-1*coverage_map), coverage_map.shape)
        #matcoord_peak_RH = np.unravel_index(np.argmax(coverage_map), coverage_map.shape)
        #ax[subj][1].plot(X[matcoord_peak_LH], Y[matcoord_peak_LH], marker='x', color='r')
        #ax[subj][1].plot(X[matcoord_peak_RH], Y[matcoord_peak_RH], marker='x', color='r')
    plt.tight_layout(w_pad=0)