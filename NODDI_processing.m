%% include NODDI toolbox in directory
addpath(genpath('/usr/local/NODDI_tool'))

%% set FSL environment
setenv('PATH', [getenv('PATH') ':/Users/alexw/fsl/bin']);

%% open dataset directory
cd('/Users/alexw/Downloads/NODDI Project/wetransfer_noddi_2023-06-02_0338/SAH_NODDI/DICOM/')

%% Pre-processing
