% Configuration file
% Defines the locations of data and figure directories

% Branch based on whether we're on Surly (computer() == 'GLNXA64') or another machine
if strcmp(computer(), 'GLNXA64')
    % data_dir contains most of the preprocessed data
    data_dir = '/home/surly-raid3/dguest-data/Guestetal2021_data/';

    % nsd_dir contains the preprocessed NSD data
    nsd_dir = '/home/surly-raid4/kendrick-data/nsd/nsddata/';
    nsdbeta_dir = '/home/surly-raid4/kendrick-data/nsd/nsddata_betas/';

    % addtl_dir contains the path to some additional files needed for full replication
    % NOTE: currently, this is not available in the public NSD release?
    addtl_dir = '~kendrick/ext/figurefiles/nsd';

    % plot_dir and plot_dir_bulk are where most intermediate figure raster files are saved;
    % plot_dir is also where .svg files used to build the figures are saved.
    plot_dir = 'figures';
    plot_dir_bulk = 'figs_output_bulk';
else
    % data_dir contains most of the preprocessed data
    data_dir = 'data';

    % nsd_dir contains the preprocessed NSD data
    nsd_dir = 'nsd/nsddata';
    nsdbeta_dir = 'nsd/nsddata_betas';

    % addtl_dir contains the path to some additional files needed for full replication
    % NOTE: currently, this is not available in the public NSD release?
    addtl_dir = '~kendrick/ext/figurefiles/nsd';

    % plot_dir and plot_dir_bulk are where most intermediate figure raster files are saved;
    % plot_dir is also where .svg files used to build the figures are saved.
    plot_dir = 'figures';
    plot_dir_bulk = 'figures_bulk';
end