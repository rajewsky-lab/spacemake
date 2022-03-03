#!/usr/bin/env python
# coding: utf-8

# In[318]:


import pandas as pd
import numpy as np
import cv2
import scanpy as sc
from sklearn.mixture import GaussianMixture


def create_visium_bead_img(adata):
    width = height = 6500
    bead_diameter = 55
    width = width + bead_diameter
    height = height + bead_diameter
    
    bead_img = np.zeros((height, width), np.uint8)
    
    max_x_pos = adata.obs.x_pos.max()
    max_y_pos = adata.obs.y_pos.max()
    
    adata.obs.x_pos = adata.obs.x_pos * height / max_x_pos + int(bead_diameter/2)
    adata.obs.y_pos = adata.obs.y_pos * width / max_y_pos + int(bead_diameter/2)
    
    for i, row in adata.obs.iterrows():
        #if row['tissue']:
        cv2.circle(bead_img, (int(row['y_pos']), int(row['x_pos'])), int(bead_diameter/2), (255,255,255), -1)
        
    cv2.imwrite('tmp.png', bead_img)
     
    bead_img = 255 - bead_img
    
    bead_img = cv2.resize(bead_img, (1000, 1000), cv2.INTER_AREA)
    
    bead_img[bead_img<255] = 0
    
    return bead_img


def create_seq_scope_bead_img(adata, filter_percentage=70):
    seq_scope_df = adata.obs.copy()
    # load the clusters in their places as per coordinates
    max_x = seq_scope_df.x_pos.max()
    max_y = seq_scope_df.y_pos.max()

    img = np.zeros((int(max_x)+1, int(max_y)+1), np.uint8)

    img[seq_scope_df.x_pos, seq_scope_df.y_pos] = seq_scope_df.total_counts
    
    
    scale_f = int(img.shape[0]/500)
    
    img_scaled=np.add.reduceat(img, range(0, img.shape[0],scale_f))
    img_scaled=np.add.reduceat(img_scaled, range(0, img.shape[1],scale_f), axis=1)

    # get only top 20%
    filter_val = np.percentile(img_scaled.flatten(), filter_percentage)
    
    img_scaled_bw = img_scaled.copy()
    
    img_scaled_bw[img_scaled_bw > filter_val] = 255
    img_scaled_bw[img_scaled_bw <=filter_val] = 0

    img_scaled_bw=img_scaled_bw.astype(np.uint8)
    
    # fill black
    img_scaled_bw = img_scaled_bw | fill_image(img_scaled_bw, 15, 5)
    # fill white
    img_scaled_bw = ~img_scaled_bw
    img_scaled_bw = img_scaled_bw | fill_image(img_scaled_bw, 21, 5)
    img_scaled_bw = ~img_scaled_bw
    # fill black

    img_scaled = img_scaled * 255 / img_scaled.max()
    img_scaled = img_scaled.astype(np.uint8)
    
    return ~img_scaled, ~img_scaled_bw

def fill_image(X, n_neighbors=3, infer_square_side=3):
    # create a zero array with + 1 row/column at the start and at the end
    # completely white image
    width, height = X.shape
    Y = np.zeros((width + infer_square_side, height + infer_square_side), np.uint16)

    for i in range(infer_square_side):
        for j in range(infer_square_side):
            Y[i:(width+i), j:(height+j)] = Y[i:(width+i), j:(height+j)] + X
            
    # get the middle
    Y = Y[1:(width+1), 1:(height+1)] - X
    
    # get only the pixels which have at least n_neighbors black neighbors and set them black
    Y[np.where(Y < 255*n_neighbors)] = 0
    Y[np.where(Y >= 255*n_neighbors)] = 255
    Y = Y.astype(np.uint8)

    return Y


def fill_holes_by_neighbors(dat, dist=2):
    top = np.zeros(dat.shape, dtype=np.uint8)
    bottom = np.zeros(dat.shape, dtype=np.uint8)
    left = np.zeros(dat.shape, dtype=np.uint8)
    right = np.zeros(dat.shape, dtype=np.uint8)
    h,w = dat.shape
    
    for i in range(dist):
        top[0:(h-i), :] = top[0:(h-i), :] | dat[i:h, :]
        bottom[i:h, :] = bottom[i:h, :] | dat[0:(h-i), :]
        right[:, 0:(w-i)] = right[:, 0:(w-i)] | dat[:, i:w]
        left[:, i:w] = right[:, i:w] | dat[:, 0:(w-i)]
    
    return top & bottom & left & right

def load_he_img(he_path, bw_threshold=None):
    he = cv2.imread(he_path, cv2.IMREAD_COLOR)
    he_gray = cv2.cvtColor(he, cv2.COLOR_BGR2GRAY)

    # create binary image
    if bw_threshold is None:
        thresh, he_bw = cv2.threshold(he_gray, 127, 255, cv2.THRESH_BINARY | cv2.THRESH_OTSU)
    else:
        he_bw = cv2.threshold(he_gray, bw_threshold, 255, cv2.THRESH_BINARY)[1]
        
    return he, he_gray, he_bw

def match_he_img(he_path, bead_img, bw_threshold=None, use_bw=True):
    # get contour of bead_img
    bead_img_cnt = np.where(bead_img < 255)
    bead_img_cropped = bead_img[bead_img_cnt[0].min():bead_img_cnt[0].max(),
                                 bead_img_cnt[1].min():bead_img_cnt[1].max()]
    
    he, he_gray, he_bw = load_he_img(he_path, bw_threshold=bw_threshold)
    cv2.imwrite('he_gray.png', he_gray)
    cv2.imwrite('he_bw.png', he_bw)
    cv2.imwrite('he.png', he)
    
    he_orig = he.copy()
    
    # find the scale by which we resize the images
    bead_img_cropped_height, bead_img_cropped_width = bead_img_cropped.shape
    he_height, he_width = he_gray.shape
    
    max_zoom = 3
    
    height_ratio = he_height / (max_zoom*bead_img_cropped_height)
    width_ratio = he_width / (max_zoom*bead_img_cropped_width)
    
    resize_type = cv2.INTER_NEAREST if use_bw else cv2.INTER_AREA
    
    if height_ratio > 1.0 and width_ratio > 1.0:
        # we scale the he image
        if height_ratio > width_ratio:
            scale_f = 1.0/width_ratio
        else:
            scale_f = 1.0/height_ratio
        
        dim = (int(he_width*scale_f)+1,
               int(he_height*scale_f)+1)
        
        he = cv2.resize(he, dim, 0,0,cv2.INTER_AREA)
        he_gray = cv2.resize(he_gray, dim, 0,0,cv2.INTER_AREA)
        he_bw = cv2.resize(he_bw, dim, 0,0,cv2.INTER_NEAREST)
        
    else:
        if height_ratio > width_ratio:
            scale_f = width_ratio
        else:
            scale_f = height_ratio
        
        bead_img_cropped = cv2.resize(bead_img_cropped,
                                      (int(bead_img_cropped_width*scale_f),
                                       int(bead_img_cropped_height*scale_f)),
                                      resize_type)
    
    he_bw = ~he_bw
    
    for i in range(10):
        he_bw = fill_holes_by_neighbors(he_bw, int(max(he_bw.shape)/200) )
        
    #for i in range(2):
        #he_bw = he_bw | fill_image(he_bw, 25, 10)
    
    he_bw = ~he_bw
    
    highest_cor = 0.0
    scale_dim = None
    match_res = None
    rotate = 0
    flip = False
    
    he_res =None
        
    for make_flip in [False, True]:
        he_copy = he_bw.copy() if use_bw else he_gray.copy()
        
        if make_flip:
            he_copy = cv2.flip(he_copy, 0)
        
        for rotate_n in [0,1,2,3]:
            he_rotated = he_copy
            for r in range(rotate_n):
                he_rotated = cv2.rotate(he_rotated, cv2.ROTATE_90_CLOCKWISE).copy()    
            
            for scale_x in np.linspace(1.0, 1.0/max_zoom + 0.01, 50):
                scales_y = (min(1.0, scale_x*1.1), max(1.0/max_zoom + 0.01, scale_x*0.9))
                for scale_y in np.linspace(scales_y[0], scales_y[1], 10):
                    dim = (int(he_bw.shape[1] * scale_x),
                           int(he_bw.shape[0] * scale_y))
                    
                    he_scaled = cv2.resize(he_rotated, dim, 0,0, resize_type)
                    m_res = cv2.matchTemplate(he_scaled, bead_img_cropped, cv2.TM_CCOEFF_NORMED)
                    
                    if m_res.max() > highest_cor:
                        highest_cor = m_res.max()
                        scale_dim = dim
                        match_res = m_res
                        rotate = rotate_n
                        flip = make_flip
    
    # find the boundaries of the match
    min_val, max_val, min_loc, max_loc = cv2.minMaxLoc(match_res)
    top_left = max_loc
    bottom_right = (top_left[0] + bead_img_cropped.shape[1],
                    top_left[1] + bead_img_cropped.shape[0])
    
    #bead_img_cropped[bead_img_cropped > 0] = 255
    bead_img_cropped_clr = cv2.cvtColor(bead_img_cropped, cv2.COLOR_GRAY2BGR)
    
    if flip:
        he = cv2.flip(he, 0)
        he_gray = cv2.flip(he_gray, 0)
        he_bw = cv2.flip(he_bw, 0)
        he_orig = cv2.flip(he_orig, 0)
    
    for r in range(rotate):
        he = cv2.rotate(he, cv2.ROTATE_90_CLOCKWISE)
        he_gray = cv2.rotate(he_gray, cv2.ROTATE_90_CLOCKWISE)
        he_bw = cv2.rotate(he_bw, cv2.ROTATE_90_CLOCKWISE)
        he_orig = cv2.rotate(he_orig, cv2.ROTATE_90_CLOCKWISE)
    
    he = cv2.resize(he, scale_dim, 0,0,cv2.INTER_AREA)
    he_bw_res = cv2.resize(he_bw, scale_dim, 0,0, cv2.INTER_NEAREST)
    he_gray_res = cv2.resize(he_gray, scale_dim, 0,0, cv2.INTER_AREA)

    return highest_cor, he, he_orig, top_left, bottom_right


# In[ ]:





# In[319]:


def match_he_visium(adata, he_path):
    bead_img = create_visium_bead_img(adata)

    bead_img_cnt = np.where(bead_img < 255)
    bead_img = bead_img[bead_img_cnt[0].min():bead_img_cnt[0].max(),
                        bead_img_cnt[1].min():bead_img_cnt[1].max()]

    highest_cor, he_res, he_orig, tl, br = match_he_img(
        he_path,
        bead_img = bead_img,
        use_bw = True)

    he_res_ratio = he_res.shape[1] / he_res.shape[0]

    he_orig = cv2.resize(he_orig,
                     (int(he_orig.shape[0] * he_res_ratio),he_orig.shape[0]),
                     cv2.INTER_AREA)

    horizontal_ratio = he_orig.shape[0] / he_res.shape[0]
    vertical_ratio = he_orig.shape[1] / he_res.shape[1]
    
    box_tl = int(tl[0] * horizontal_ratio), int(tl[1] * vertical_ratio)
    box_br = int(br[0] * horizontal_ratio), int(br[1] * vertical_ratio)

    bead_img_resized = cv2.resize(bead_img,
                                  (box_br[0] - box_tl[0], box_br[1] - box_tl[1]),
                                   cv2.INTER_NEAREST)

    match_bw = cv2.cvtColor(bead_img_resized, cv2.COLOR_GRAY2BGR)

    he_match = he_orig[box_tl[1]:box_br[1], box_tl[0]:box_br[0]].copy()
    he_orig[box_tl[1]:box_br[1], box_tl[0]:box_br[0]] &= match_bw

    cv2.rectangle(he_orig, box_tl, box_br,  (180, 233, 86), 3)
    
    return he_orig, he_match


# In[ ]:





# In[320]:


from spacemake.config import ConfigFile
from spacemake.project_df import ProjectDF

class SpaceMake:
    def __init__(self, root):
        self.root = root
        self.config = ConfigFile.from_yaml(f'{root}/config.yaml')
        self.project_df = ProjectDF(f'{root}/project_df.csv', config=self.config)
    
    def load_processed_adata(self,
        project,
        sample,
        run_mode_name,
        umi_cutoff 
    ):
        run_mode = self.config.get_run_mode(run_mode_name)
        
        adata = sc.read(f'{self.root}/projects/{project}/processed_data/{sample}/'                   f'illumina/complete_data/automated_analysis/{run_mode_name}/'                   f'umi_cutoff_{umi_cutoff}/results.h5ad')
        
        adata.uns['run_mode_variables'] = run_mode.variables
        adata.uns['puck_variables'] = self.project_df.get_puck_variables(
            project_id = project,
            sample_id = sample)
        
        return adata
    
    def load_raw_spatial_adata(
        self,
        project,
        sample,
        run_mode_name
    ):
        run_mode = self.config.get_run_mode(run_mode_name)

        dge_type = ""
        dge_cleaned = ""
        polyA_adapter_trimmed = ""
        mm_included = ""

        if run_mode.variables["polyA_adapter_trimming"]:
            polyA_adapter_trimmed = ".polyA_adapter_trimmed"

        if run_mode.variables["count_intronic_reads"]:
            dge_type = ".all"
        else:
            dge_type = ".exon"

        if run_mode.variables["count_mm_reads"]:
            mm_included = ".mm_included"

        if run_mode.variables["clean_dge"]:
            dge_cleaned = ".cleaned"

        adata = sc.read(f'{self.root}/projects/{project}/processed_data/{sample}/'            f'illumina/complete_data/dge/dge{dge_type}{dge_cleaned}'            f'{polyA_adapter_trimmed}{mm_included}.spatial_beads.h5ad')
        
        adata.uns['run_mode_variables'] = run_mode.variables
        adata.uns['puck_variables'] = self.project_df.get_puck_variables(
            project_id = project,
            sample_id = sample)
        
        return adata


# In[ ]:





# In[321]:


def attach_he_adata(
    adata,
    matched_he,
    adata_raw=None,
    push_by_bead_diameter=True,
):
    import math
    # get the width of the puck
    puck_width_um = adata.uns['puck_variables']['width_um']
    bead_diameter_um = adata.uns['puck_variables']['spot_diameter_um']
    
    if adata.uns['run_mode_variables']['mesh_data']:
        bead_diameter_um = math.sqrt(3) * adata.uns['run_mode_variables']['mesh_spot_diameter_um']
    
    if adata_raw is None:
        adata_raw = adata
    
    bead_width = adata_raw.obsm['spatial'].max(axis=0)[0] - adata_raw.obsm['spatial'].min(axis=0)[0]

    bead_distance_um = puck_width_um / bead_width


    # get the width of the processed adata
    adata_width_um = bead_distance_um *        (adata.obsm['spatial'].max(axis=0)[0] - adata.obsm['spatial'].min(axis=0)[0])

    # rotate he to align with coordinate system of scanpy
    rotated_he = cv2.flip(cv2.rotate(matched_he, cv2.ROTATE_90_CLOCKWISE), 1)
    
    h_px, w_px = rotated_he.shape[:2]

    px_per_um = w_px / adata_width_um 

    # convert so that origin at 0, 0
    locations = adata.obsm['spatial']
    locations = locations - locations.min(axis=0)

    bead_diameter_px = bead_diameter_um * px_per_um

    locations = locations * [(w_px ) / locations.max(axis=0)[0],
                             (h_px ) / locations.max(axis=0)[1]]
    
    if push_by_bead_diameter:
        locations = [bead_diameter_px, bead_diameter_px] + locations
        
    locations = locations.astype('int')

    adata.uns['spatial']={
        'img': {'images': {
                    'hires': rotated_he
                    },
                'scalefactors': {'tissue_hires_scalef': 1,
                                 'spot_diameter_fullres': bead_diameter_px}} 
    }
    adata.obsm['spatial'] = locations
    
    return adata


# In[322]:


import scanpy as sc
sc.set_figure_params(dpi=300)

spmk = SpaceMake('/data/rajewsky/home/tsztank/projects/spatial/sts-paper')

cluster_clrs = ["#0000FF","#FF0000","#00FF00","#000033","#FF00B6","#005300","#FFD300",
"#009FFF","#9A4D42","#00FFBE","#783FC1","#1F9698","#FFACFD","#B1CC71",
"#F1085C","#FE8F42","#DD00FF","#201A01","#720055","#766C95","#02AD24",
"#C8FF00","#886C00","#FFB79F","#858567","#A10300","#14F9FF","#00479E",
"#DC5E93","#93D4FF","#004CFF","#F2F318"]

cluster_clrs = {str(i): cluster_clrs[i] for i in range(len(cluster_clrs))}

adata = spmk.load_processed_adata(
    project = 'visium',
    sample = 'public_1',
    run_mode_name = 'exon',
    umi_cutoff = 10000
)
adata_raw = spmk.load_raw_spatial_adata(
    'visium', 'public_1', 'exon'
)

he, matched_he = match_he_visium(adata, './V1_Adult_Mouse_Brain_image_small.tif')
cv2.imwrite('mouse_coronal.png', he)

adata = attach_he_adata(adata, matched_he, adata_raw)


# In[ ]:


plt = sc.pl.spatial(adata, color='leiden_1.2', palette=cluster_clrs, return_fig=True, show=False, title='')
                 
plt[0].invert_yaxis()
plt[0].figure
plt[0].set_xlabel('')
plt[0].set_ylabel('')
ticks = [x * 19.5 for x in range(1, 7, 2)]
ticks_labels = ['2mm', '4mm', '6mm']

ticks_x = [x * 7 for x in ticks]
ticks_y = [x * 6.5 for x in ticks]
plt[0].set_xticks(ticks_x)
plt[0].set_xticklabels(ticks_labels)
plt[0].set_yticks(ticks_y)
plt[0].set_yticklabels(ticks_labels)
plt[0].grid(False)


# In[ ]:




plt = sc.pl.spatial(adata, color='leiden_1.2', groups=['21', '22', '23'],  palette=cluster_clrs, return_fig=True,
                   show = False, title='')
plt[0].invert_yaxis()
plt[0].figure

plt[0].set_xlabel('')
plt[0].set_ylabel('')
ticks = [x * 19.5 for x in range(1, 7, 2)]
ticks_labels = ['2mm', '4mm', '6mm']

ticks_x = [x * 7 for x in ticks]
ticks_y = [x * 6.5 for x in ticks]
plt[0].set_xticks(ticks_x)
plt[0].set_xticklabels(ticks_labels)
plt[0].set_yticks(ticks_y)
plt[0].set_yticklabels(ticks_labels)
plt[0].grid(False)


# In[ ]:


plt = sc.pl.spatial(adata, color='leiden_1.2', palette=cluster_clrs,
              return_fig=True, show=False, title='', img_key=None)

plt[0].set_facecolor((240/256,240/256,240/256))
plt[0].invert_yaxis()
plt[0].figure
plt[0].set_xlabel('')
plt[0].set_ylabel('')
ticks = [x * 19.5 for x in range(1, 7, 2)]
ticks_labels = ['2mm', '4mm', '6mm']

ticks_x = [x * 7 for x in ticks]
ticks_y = [x * 6.5 for x in ticks]
plt[0].set_xticks(ticks_x)
plt[0].set_xticklabels(ticks_labels)
plt[0].set_yticks(ticks_y)
plt[0].set_yticklabels(ticks_labels)
plt[0].grid(False)


# In[ ]:





# In[ ]:


xs


# In[ ]:


px_to_pos = adata.obsm['spatial'].max(axis=0)[1]/(tmp.obs.y_pos.max() - tmp.obs.y_pos.min())
px_to_pos


# In[ ]:





# In[ ]:


adata = spmk.load_processed_adata(
    project = 'visium',
    sample = 'visium_mouse_brain_sagittal_anterior',
    run_mode_name = 'exon_mm',
    umi_cutoff = 12000
)

adata_raw = spmk.load_raw_spatial_adata(
    'visium', 'visium_mouse_brain_sagittal_anterior', 'exon_mm'
)

he, matched_he = match_he_visium(adata, './V1_Mouse_Brain_Sagittal_Anterior_image_sm.tif')
cv2.imwrite('mouse_saggital.png', matched_he)

adata = attach_he_adata(adata, matched_he, adata_raw)

sc.pl.spatial(adata, color='leiden_1.2', palette=cluster_clrs)


# In[ ]:





# In[ ]:


adata = load_processed_adata(
    root = root,
    project = 'visium',
    sample = 'visium_mouse_kidney',
    run_mode = 'exon_mm',
    umi_cutoff = 12000
)

matched_he = match_he_visium(adata, './V1_Mouse_Kidney_image_sm.tif')
cv2.imwrite('mouse_kidney.png', matched_he)


# In[340]:


def match_he_seq_scope(adata, he_path, suffix='',
                       bw_threshold=200, filter_percentage=70,
                       box_size=0.5):
    bead_img, bead_img_bw = create_seq_scope_bead_img(adata, filter_percentage)

    h, w = bead_img.shape

    bottom_right = (int(w * ((1+box_size)/2)),
                    int(h * ((1+box_size)/2)))
    top_left = (int(w * ((1-box_size) / 2)),
                int(h * ((1-box_size) / 2)))

    # save intermediate files
    tmp_img = cv2.cvtColor(bead_img.copy(), cv2.COLOR_GRAY2BGR)
    cv2.imwrite(f'bead_img_{suffix}.png', tmp_img)
    cv2.rectangle(tmp_img, bottom_right, top_left, (180, 233, 86), 2)
    cv2.imwrite(f'bead_img_box_{suffix}.png', tmp_img)
    
    tmp_img = cv2.cvtColor(bead_img_bw.copy(), cv2.COLOR_GRAY2BGR)
    cv2.imwrite(f'bead_img_bw_{suffix}.png', tmp_img)
    cv2.rectangle(tmp_img, bottom_right, top_left, (180, 233, 86), 2)
    cv2.imwrite(f'bead_img_bw_box_{suffix}.png', tmp_img)

    img= bead_img_bw[
        top_left[1]:bottom_right[1],
        top_left[0]:bottom_right[0]    
    ]


    highest_cor, he_res, he_orig, tl, br = match_he_img(
        he_path,
        bead_img=img,
        bw_threshold=bw_threshold,
        use_bw = True
    )

    match_w = br[0]-tl[0]
    match_h = br[1]-tl[1]
    multiplier = (1-box_size) / (box_size * 2)
    box_tl = (int(tl[0]-match_w * multiplier), int(tl[1]-match_h * multiplier))
    #box_tl = h,w
    box_br = (int(br[0]+match_w * multiplier), int(br[1]+match_h * multiplier))

    he_res_ratio = he_res.shape[1] / he_res.shape[0]

    he_orig = cv2.resize(he_orig,
                     (int(he_orig.shape[0] * he_res_ratio),he_orig.shape[0]),
                     cv2.INTER_AREA)

    horizontal_ratio = he_orig.shape[0] / he_res.shape[0]
    vertical_ratio = he_orig.shape[1] / he_res.shape[1]

    # scale boxes
    tl = int(tl[0] * horizontal_ratio), int(tl[1] * vertical_ratio)
    br = int(br[0] * horizontal_ratio), int(br[1] * vertical_ratio)
    box_tl = (int(box_tl[0] * horizontal_ratio),
              int(box_tl[1] * vertical_ratio))
    box_br = (int(box_br[0] * horizontal_ratio),
              int(box_br[1] * vertical_ratio))

    bead_img_resized = cv2.resize(bead_img_bw,
                                  (box_br[0] - box_tl[0], box_br[1] - box_tl[1]),
                                   cv2.INTER_NEAREST)

    match_bw = cv2.cvtColor(bead_img_resized, cv2.COLOR_GRAY2BGR)

    he_match = he_orig[box_tl[1]:box_br[1], box_tl[0]:box_br[0]].copy()

    cv2.rectangle(he_orig, tl, br,  (180, 233, 86, 255), 5)
    cv2.rectangle(he_orig, box_tl, box_br,  (180, 233, 86, 255), 5)
    
    overlay = he_orig.copy()
    overlay[box_tl[1]:box_br[1], box_tl[0]:box_br[0]] &= match_bw
    
    cv2.addWeighted(overlay, 0.4, he_orig, 0.6, 0, he_orig)
    
    return he_orig, he_match


# In[ ]:


adata = load_raw_spatial_adata(root,
    'seq_scope', 'seq_scope_liver_2105', 'seq_scope')


# In[ ]:


spmk.project_df.get_metadata(
            "puck", project_id='seq_scope', sample_id='seq_scope_liver_2105'
        )


# In[343]:


adata = spmk.load_processed_adata(
    project = 'seq_scope',
    sample = 'seq_scope_liver_2106',
    run_mode_name = 'seq_scope',
    umi_cutoff = 300
)

adata_raw = spmk.load_raw_spatial_adata(
    'seq_scope', 'seq_scope_liver_2106', 'seq_scope'
)

he, matched_he = match_he_seq_scope(adata_raw,
                                './wt_4X_1.jpg',
                                '2106', bw_threshold=180, box_size=0.7)

cv2.imwrite('matched_he_2106.png', he)


# In[ ]:


matched_he.shape


# In[ ]:


#adata = attach_he_adata(adata, matched_he, adata_raw)
sc.pl.spatial(adata,color='total_counts',palette=cluster_clrs)
sc.pl.spatial(adata,palette=cluster_clrs)


# In[ ]:


adata = load_raw_adata(root, 'seq_scope', 'seq_scope_liver_2106')

matched_he = match_he_seq_scope(adata,
                                './wt_4X_1.jpg',
                                '2106',
                                bw_threshold=190,
                                box_size=0.7)
cv2.imwrite('matched_he_2106.png', matched_he)


# In[ ]:


adata = load_raw_adata(root, 'seq_scope', 'seq_scope_liver_2107')

matched_he = match_he_seq_scope(adata,
                                './wt_4X_1.jpg',
                                '2107',
                                bw_threshold=200,
                                box_size=0.6)
cv2.imwrite('matched_he_2107.png', matched_he)

