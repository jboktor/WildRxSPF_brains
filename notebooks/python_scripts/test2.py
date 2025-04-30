import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from matplotlib import cm
from matplotlib.lines import Line2D
import os
from os.path import exists, split, join, splitext
from os import makedirs
import glob
import requests
from collections import defaultdict
import nrrd
import torch
from torch.nn.functional import grid_sample
import tornado
from STalign import STalign
import copy


base_path = "/central/groups/mthomson/jboktor/spatial_genomics"
ref_path = os.path.join(base_path, 'allen_registration_ref')
analysis_path = os.path.join(base_path, 'jess_2024-01-23')
analysis_reg_path = os.path.join(analysis_path, 'data/interim/registration')

cmd_file = os.path.join(analysis_reg_path, 'cell-metadata_2024-05-10.csv')
cmd = pd.read_csv(cmd_file, index_col=0)
# Load file

df = cmd[cmd.sample_id == 'roi2_run2']
#Load xy positions and convert to um
x = np.array(df['center_x'])[1:].astype(float) * 0.106
y = np.array(df['center_y'])[1:].astype(float) * 0.106
# load rna counts
cell_reads = np.array(df['nCount_RNA'])[1:].astype(float)
# parameters
dx = 50
blur = 0.5

url = 'http://api.brain-map.org/api/v2/data/query.csv?criteria=model::Structure,rma::criteria,[ontology_id$eq1],rma::options[order$eq%27structures.graph_order%27][num_rows$eqall]'
ontology_name, namesdict = STalign.download_aba_ontology(url, 'allen_ontology.csv')  # url for adult mouse
imageurl = 'http://download.alleninstitute.org/informatics-archive/current-release/mouse_ccf/ara_nissl/ara_nissl_50.nrrd'
labelurl = 'http://download.alleninstitute.org/informatics-archive/current-release/mouse_ccf/annotation/ccf_2017/annotation_50.nrrd'
imagefile, labelfile = STalign.download_aba_image_labels(imageurl, labelurl, 'aba_nissl.nrrd', 'aba_annotation.nrrd')

# Rasterize Image
X_,Y_,W = STalign.rasterize(x,y,g=cell_reads,dx=dx, blur = blur,draw=False)

# Plot unrasterized/rasterized images
fig, ax = plt.subplots(1, 2)
ax[0].scatter(x, y, s=0.5, alpha=0.25)
ax[0].invert_yaxis()
ax[0].set_title('List of cells')
ax[0].set_aspect('equal')

W = W[0]
extent = (X_[0], X_[-1], Y_[0], Y_[-1])
ax[1].imshow(W, origin='lower')
ax[1].invert_yaxis()
ax[1].set_title('Rasterized')

# save figure
fig.canvas.draw()
fig.savefig('slice_Rasterized.png')

# Matching slice to ABA
# find slice
# peruse through images in atlas
# Loading the atlas
slice = 135

vol, hdr = nrrd.read(imagefile)
A = vol
vol, hdr = nrrd.read(labelfile)
L = vol

dxA = np.diag(hdr['space directions'])
nxA = A.shape
xA = [np.arange(n) * d - (n - 1) * d / 2.0 for n, d in zip(nxA, dxA)]
XA = np.meshgrid(*xA, indexing='ij')

fig, ax = plt.subplots(1, 2)
extentA = STalign.extent_from_x(xA[1:])
ax[0].imshow(A[slice], extent=extentA)
ax[0].set_title('Atlas Slice')

ax[1].imshow(W, extent=extentA)
ax[1].set_title('Target Image')
# plt.show()

fig.canvas.draw()
fig.savefig('slice_reference_matched.png')

from scipy.ndimage import rotate

theta_deg = 0

fig, ax = plt.subplots(1, 2)
extentA = STalign.extent_from_x(xA[1:])
ax[0].imshow(rotate(A[slice], angle=theta_deg), extent=extentA)
ax[0].set_title('Atlas Slice')

ax[1].imshow(W, extent=extentA)
ax[1].set_title('Target Image')
fig.canvas.draw()

fig.savefig('slice_reference_matched_rotated.png')

points_target = np.array([[0,0], [0, 5000]]) #(-y,-x)
points_atlas = np.array([[0,-5000], [0,5000]]) # (-y,-x)
Li,Ti = STalign.L_T_from_points(points_atlas,points_target)


xJ = [Y_, X_]
J = W[None] / np.mean(np.abs(W))
xI = xA
I = A[None] / np.mean(np.abs(A), keepdims=True)
I = np.concatenate((I, (I - np.mean(I)) ** 2))
Inorm = STalign.normalize(I, t_min=0, t_max=1)
Jnorm = STalign.normalize(J, t_min=0, t_max=1)

sigmaA = 0.1  # standard deviation of artifact intensities
sigmaB = 0.1  # standard deviation of background intensities
sigmaM = 0.1  # standard deviation of matching tissue intensities
muA = torch.tensor([0.7, 0.7, 0.7], device='cuda')  # average of artifact intensities
muB = torch.tensor([0, 0, 0], device='cuda')  # average of background intensities

fig, ax = plt.subplots()
ax.hist(Jnorm.ravel())
plt.xlabel('Intensity')
plt.ylabel('Number of Pixels')
plt.title('Intensity Histogram of Target Image')
fig.savefig('slice_intensity_histogram.png')

# param initialization
# initialize variables
scale_x = 0.9  # default = 0.9
scale_y = 0.9  # default = 0.9
scale_z = 0.9  # default = 0.9
theta0 = (np.pi / 180) * theta_deg

# get an initial guess
if 'Ti' in locals():
    T = np.array([-xI[0][slice], np.mean(xJ[0]) - (Ti[0] * scale_y), np.mean(xJ[1]) - (Ti[1] * scale_x)])
else:
    T = np.array([-xI[0][slice], np.mean(xJ[0]), np.mean(xJ[1])])
# T = np.array([-xI[0][slice], 0, 0])

scale_atlas = np.array([[scale_z, 0, 0],
                        [0, scale_x, 0],
                        [0, 0, scale_y]])
L = np.array([[1.0, 0.0, 0.0],
              [0.0, np.cos(theta0), -np.sin(theta0)],
              [0.0, np.sin(theta0), np.cos(theta0)]])
L = np.matmul(L, scale_atlas)  # np.identity(3)

# %%time
# returns mat = affine transform, v = velocity, xv = pixel locations of velocity points
transform = STalign.LDDMM_3D_to_slice(
    xI, Inorm, xJ, Jnorm,
    T=T, L=L,
    nt=4, niter=1000,
    a=250,
    device='cuda',
    sigmaA=sigmaA,  # standard deviation of artifact intensities
    sigmaB=sigmaB,  # standard deviation of background intensities
    sigmaM=sigmaM,  # standard deviation of matching tissue intensities
    muA=muA,  # average of artifact intensities
    muB=muB  # average of background intensities
)

A = transform['A']
v = transform['v']
xv = transform['xv']
Xs = transform['Xs']

df = STalign.analyze3Dalign(labelfile, xv, v, A, xJ, dx, scale_x=scale_x, scale_y=scale_y, x=x, y=y, X_=X_, Y_=Y_, namesdict=namesdict, device='cuda')

It = torch.tensor(I, device='cuda', dtype=torch.float64)
AI = STalign.interp3D(xI, It, Xs.permute(3, 0, 1, 2), padding_mode="border")
Ishow_source = ((AI - torch.amin(AI, (1, 2, 3))[..., None, None]) / (torch.amax(AI, (1, 2, 3)) - torch.amin(AI, (1, 2, 3)))[..., None, None, None]).permute(1, 2, 3, 0).clone().detach().cuda()
Jt = torch.tensor(J, device='cuda', dtype=torch.float64)
Ishow_target = Jt.permute(1, 2, 0).cuda() / torch.max(Jt).item()

# Move tensors from GPU to CPU and convert to numpy arrays for plotting
Ishow_source_cpu = Ishow_source.cpu().numpy()
Ishow_target_cpu = Ishow_target.cpu().numpy()

import matplotlib as mpl
fig, ax = plt.subplots(1, 3, figsize=(15, 5))
ax[0].imshow(Ishow_target_cpu, cmap=mpl.cm.Blues, alpha=0.9)
ax[0].set_title('STARmap Slice')
ax[1].imshow(Ishow_source_cpu[0, :, :, 0], cmap=mpl.cm.Reds, alpha=0.2)
ax[1].set_title('z=0 slice of Aligned 3D Allen Brain Atlas')
ax[2].imshow(Ishow_target_cpu, cmap=mpl.cm.Blues, alpha=0.9)
ax[2].imshow(Ishow_source_cpu[0, :, :, 0], cmap=mpl.cm.Reds, alpha=0.3)
ax[2].set_title('Overlayed')

# plt.show()
fig.savefig('slice_aligned.png')
