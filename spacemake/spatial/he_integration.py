import pandas as pd
import numpy as np
import cv2
import scanpy as sc
import logging

from spacemake.errors import SpacemakeError

logger_name = "spacemake.spatial"
logger = logging.getLogger(logger_name)

def create_spot_expression_img(
        adata,
        binary=True,
    ):
    df = adata.obs.copy()
    x_pos_max = df.x_pos.max()
    y_pos_max = df.y_pos.max()
    coord_by_um = adata.uns['puck_variables']['coord_by_um']

    width = int(x_pos_max / coord_by_um)
    height = int(y_pos_max / coord_by_um)

    spot_diameter = int(adata.uns['puck_variables']['spot_diameter_um'])

    width = width + spot_diameter
    height = height + spot_diameter
    
    spot_img = np.zeros((height, width), np.uint8)

    df.total_counts = df.total_counts * 255 / df.total_counts.max()
    
    df.x_pos = df.x_pos * height / x_pos_max + int(spot_diameter/2)
    df.y_pos = df.y_pos * width / y_pos_max + int(spot_diameter/2)
    
    for i, row in df.iterrows():
        #if row['tissue']:
        if binary:
            color = 255
        else:
            color = int(row['total_counts'])

        # draw the circle
        cv2.circle(spot_img,
                   (int(row['y_pos']), int(row['x_pos'])),
                   int(spot_diameter/2),
                   (color, color, color), -1)
        
    spot_img = 255 - spot_img
    
    spot_img = cv2.resize(spot_img, (1000, 1000), cv2.INTER_AREA)
    
    if binary:
        spot_img[spot_img<255] = 0
    
    return spot_img

def create_aggregated_expression_img(adata, filter_percentage=70):
    df = adata.obs.copy()
    # load the clusters in their places as per coordinates
    df.x_pos = df.x_pos.astype('int')
    df.y_pos = df.y_pos.astype('int')

    max_x = df.x_pos.max()
    max_y = df.y_pos.max()

    img = np.zeros((max_x+1, max_y+1), np.uint8)

    img[df.x_pos, df.y_pos] = df.total_counts
    
    if max_x <= 500 or max_y <= 500:
        raise SpacemakeError('Underlying spatial coordinates cannot be aggregated\n'+
            'Only coordinates with maximum values for x_pos and y_pos ' +
            'both higher than 500 can be aggregated\n' +
            'Use `spacemake spatial create_spot_expression_img` instead')

    scale_f = int(img.shape[0]/500)

    # creating a 500x500 pixel image
    img_scaled=np.add.reduceat(img, range(0, img.shape[0],scale_f))
    img_scaled=np.add.reduceat(img_scaled, range(0, img.shape[1],scale_f), axis=1)

    # get only top 30%, by default
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

def match_he_img(he_path, spot_img, bw_threshold=None, use_bw=True):
    # get contour of spot_img
    spot_img_cnt = np.where(spot_img < 255)
    spot_img_cropped = spot_img[spot_img_cnt[0].min():spot_img_cnt[0].max(),
                                 spot_img_cnt[1].min():spot_img_cnt[1].max()]
    
    he, he_gray, he_bw = load_he_img(he_path, bw_threshold=bw_threshold)
    
    he_orig = he.copy()
    
    # find the scale by which we resize the images
    spot_img_cropped_height, spot_img_cropped_width = spot_img_cropped.shape
    he_height, he_width = he_gray.shape
    
    max_zoom = 3
    
    height_ratio = he_height / (max_zoom*spot_img_cropped_height)
    width_ratio = he_width / (max_zoom*spot_img_cropped_width)
    
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
        
        spot_img_cropped = cv2.resize(spot_img_cropped,
                                      (int(spot_img_cropped_width*scale_f),
                                       int(spot_img_cropped_height*scale_f)),
                                      resize_type)
    
    he_bw = ~he_bw
    
    for i in range(10):
        he_bw = fill_holes_by_neighbors(he_bw, int(max(he_bw.shape)/200))
    
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
                    m_res = cv2.matchTemplate(he_scaled, spot_img_cropped, cv2.TM_CCOEFF_NORMED)
                    
                    if m_res.max() > highest_cor:
                        highest_cor = m_res.max()
                        scale_dim = dim
                        match_res = m_res
                        rotate = rotate_n
                        flip = make_flip
    
    # find the boundaries of the match
    min_val, max_val, min_loc, max_loc = cv2.minMaxLoc(match_res)
    top_left = max_loc
    bottom_right = (top_left[0] + spot_img_cropped.shape[1],
                    top_left[1] + spot_img_cropped.shape[0])
    
    #spot_img_cropped[spot_img_cropped > 0] = 255
    spot_img_cropped_clr = cv2.cvtColor(spot_img_cropped, cv2.COLOR_GRAY2BGR)
    
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

def match_he_spot_img(adata, he_path):
    spot_img = create_expression_spot_img(adata)

    spot_img_cnt = np.where(spot_img < 255)
    spot_img = spot_img[spot_img_cnt[0].min():spot_img_cnt[0].max(),
                        spot_img_cnt[1].min():spot_img_cnt[1].max()]

    highest_cor, he_res, he_orig, tl, br = match_he_img(
        he_path,
        spot_img = spot_img,
        use_bw = True)

    he_res_ratio = he_res.shape[1] / he_res.shape[0]

    he_orig = cv2.resize(he_orig,
                     (int(he_orig.shape[0] * he_res_ratio),he_orig.shape[0]),
                     cv2.INTER_AREA)

    horizontal_ratio = he_orig.shape[0] / he_res.shape[0]
    vertical_ratio = he_orig.shape[1] / he_res.shape[1]
    
    box_tl = int(tl[0] * horizontal_ratio), int(tl[1] * vertical_ratio)
    box_br = int(br[0] * horizontal_ratio), int(br[1] * vertical_ratio)

    spot_img_resized = cv2.resize(spot_img,
                                  (box_br[0] - box_tl[0], box_br[1] - box_tl[1]),
                                   cv2.INTER_NEAREST)

    match_bw = cv2.cvtColor(spot_img_resized, cv2.COLOR_GRAY2BGR)

    he_match = he_orig[box_tl[1]:box_br[1], box_tl[0]:box_br[0]].copy()
    he_orig[box_tl[1]:box_br[1], box_tl[0]:box_br[0]] &= match_bw

    cv2.rectangle(he_orig, box_tl, box_br,  (180, 233, 86), 3)
    
    return he_orig, he_match

def match_he_aggregated_img(
        adata, he_path,
        bw_threshold=200, filter_percentage=70,
        box_size=0.5
    ):
    spot_img, spot_img_bw = create_aggregated_expression_img(adata, filter_percentage)

    h, w = spot_img.shape

    bottom_right = (int(w * ((1+box_size)/2)),
                    int(h * ((1+box_size)/2)))
    top_left = (int(w * ((1-box_size) / 2)),
                int(h * ((1-box_size) / 2)))

    # save intermediate files
    tmp_img = cv2.cvtColor(spot_img.copy(), cv2.COLOR_GRAY2BGR)
    cv2.rectangle(tmp_img, bottom_right, top_left, (180, 233, 86), 2)
    
    tmp_img = cv2.cvtColor(spot_img_bw.copy(), cv2.COLOR_GRAY2BGR)
    cv2.rectangle(tmp_img, bottom_right, top_left, (180, 233, 86), 2)

    img= spot_img_bw[
        top_left[1]:bottom_right[1],
        top_left[0]:bottom_right[0]    
    ]

    highest_cor, he_res, he_orig, tl, br = match_he_img(
        he_path,
        spot_img=img,
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

    spot_img_resized = cv2.resize(spot_img_bw,
                                  (box_br[0] - box_tl[0], box_br[1] - box_tl[1]),
                                   cv2.INTER_NEAREST)

    match_bw = cv2.cvtColor(spot_img_resized, cv2.COLOR_GRAY2BGR)

    he_match = he_orig[box_tl[1]:box_br[1], box_tl[0]:box_br[0]].copy()

    cv2.rectangle(he_orig, tl, br,  (180, 233, 86, 255), 5)
    cv2.rectangle(he_orig, box_tl, box_br,  (180, 233, 86, 255), 5)
    
    overlay = he_orig.copy()
    overlay[box_tl[1]:box_br[1], box_tl[0]:box_br[0]] &= match_bw
    
    cv2.addWeighted(overlay, 0.4, he_orig, 0.6, 0, he_orig)
    
    return he_orig, he_match

def attach_he_adata(
    adata,
    matched_he,
    adata_raw=None,
    push_by_spot_diameter=True,
    raw_matched=False
):
    import math
    # get the width of the puck
    puck_width_um = adata.uns['puck_variables']['width_um']
    spot_diameter_um = adata.uns['puck_variables']['spot_diameter_um']
    
    if adata.uns['run_mode_variables']['mesh_data']:
        spot_diameter_um = math.sqrt(3) * adata.uns['run_mode_variables']['mesh_spot_diameter_um']
    
    if adata_raw is None:
        adata_raw = adata
    
    spot_width = adata_raw.obsm['spatial'].max(axis=0)[0]

    spot_distance_um = puck_width_um / spot_width

    # get the width of the processed adata
    adata_width_um = spot_distance_um *\
        (adata.obsm['spatial'].max(axis=0)[0])

    # rotate he to align with coordinate system of scanpy
    rotated_he = cv2.flip(cv2.rotate(matched_he, cv2.ROTATE_90_CLOCKWISE), 1)
    
    h_px, w_px = rotated_he.shape[:2]

    px_per_um = w_px / adata_width_um 

    # convert so that origin at 0, 0
    locations_raw = adata_raw.obsm['spatial']
    locations = adata.obsm['spatial']
    
    if not raw_matched:
        locations = locations - locations.min(axis=0)

    spot_diameter_px = spot_diameter_um * px_per_um

    if raw_matched:
        locations = locations * [(w_px ) / locations_raw.max(axis=0)[0],
                             (h_px ) / locations_raw.max(axis=0)[1]]
    else:
        locations = locations * [(w_px ) / locations.max(axis=0)[0],
                             (h_px ) / locations.max(axis=0)[1]]
    
    if push_by_spot_diameter:
        locations = [spot_diameter_px, spot_diameter_px] + locations
        
    locations = locations.astype('int')

    adata.uns['spatial']={
        'img': {'images': {
                    'hires': rotated_he
                    },
                'scalefactors': {'tissue_hires_scalef': 1,
                                 'spot_diameter_fullres': spot_diameter_px}} 
    }
    adata.obsm['spatial'] = locations
    
    return adata
