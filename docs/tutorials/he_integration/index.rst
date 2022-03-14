H&E integration tutorial
========================

Downloading tutorial data
-------------------------

First make sure you have installed spacemake as specified :ref:`here <installation>`.

Second, we need to download the pre-processed spacemake data.

To do this, use example datasets from here:
https://github.com/rajewsky-lab/spacemake-test-data

To download the data we do:

.. code:: ipython3

    !git clone https://github.com/rajewsky-lab/spacemake-test-data


.. parsed-literal::

    Cloning into 'spacemake-test-data'...
    remote: Enumerating objects: 21, done.[K
    remote: Counting objects: 100% (21/21), done.[K
    remote: Compressing objects: 100% (19/19), done.[K
    remote: Total 21 (delta 2), reused 16 (delta 1), pack-reused 0[K
    Receiving objects: 100% (21/21), 71.47 MiB | 28.94 MiB/s, done.
    Resolving deltas: 100% (2/2), done.
    Filtering content: 100% (4/4), 219.64 MiB | 27.28 MiB/s, done.


Then we uncompress the datafiles, as scanpy can’t load compressed data:

.. code:: ipython3

    !find ./spacemake-test-data -type f -wholename '*.h5ad.gz' -exec unpigz {} \;

After we have downloaded the data, and uncompressed the objects, we load
scanpy and set some default parameters and colors.

.. code:: ipython3

    import scanpy as sc
    import cv2
    
    root_dir = 'spacemake-test-data'
    
    # we set the figure to have higher dpi
    sc.set_figure_params(dpi=300)
    
    
    cluster_clrs = ["#0000FF","#FF0000","#00FF00","#000033",
    "#FF00B6","#005300","#FFD300","#009FFF","#9A4D42",
    "#00FFBE","#783FC1","#1F9698","#FFACFD","#B1CC71",
    "#F1085C","#FE8F42","#DD00FF","#201A01","#720055",
    "#766C95","#02AD24","#C8FF00","#886C00","#FFB79F",
    "#858567","#A10300","#14F9FF","#00479E","#DC5E93",
    "#93D4FF","#004CFF","#F2F318"]
    
    cluster_clrs = {str(i): cluster_clrs[i] for i in range(len(cluster_clrs))}

Matching Visium data
--------------------

This Visium dataset showed here is an adult mouse brain coronal section
dataset downloaded from here
https://www.10xgenomics.com/resources/datasets, processed with spacemake
using the visium run mode as specified
:ref:`here <provided run\\_mode(s)>`.

First we load the processed and filtered visium data:

.. code:: ipython3

    adata_visium = sc.read(f'{root_dir}/visium/adata.h5ad')

Then, we can match the H&E image with our data. To do this we use the :func:`spacemake.spatial.he_integration.match_he_spot_img` function from spacemake. 

.. code:: ipython3

    from spacemake.spatial.he_integration import match_he_spot_img
    
    he, matched_he = match_he_spot_img(adata_visium,
                                       f'{root_dir}/visium/V1_Adult_Mouse_Brain_image_small.tif')
    
    cv2.imwrite('mouse_coronal.png', he)


.. parsed-literal::

    True



After we matched the data, we can attach the result image to the
original data which then we can use for plotting. For this, we use the
:func:`spacemake.spatial.he_integration.attach_he_adata` function.

.. code:: ipython3

    from spacemake.spatial.he_integration import attach_he_adata
    
    adata_visium_attached = attach_he_adata(adata_visium.copy(), matched_he)
    
    plt = sc.pl.spatial(adata_visium_attached,
                        color='leiden_1.2',
                        palette=cluster_clrs,
                        return_fig=True,
                        show=False,
                        title='')
    
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



.. image:: output_11_0.png
   :width: 1423px


Match Seq-scope data
--------------------

The Seq-scope dataset
(https://www.sciencedirect.com/science/article/abs/pii/S0092867421006279)
shown here is tile 2107. It was processed with spacemake using the
seq_scope run mode as specified :ref:`here <provided run\\_mode(s)>`.
In this run mode the data will be meshed into a hexagonal meshgrid of 10
micron side equal hexagons.

Similar to Visium, we first need to load the data:

.. code:: ipython3

    adata_seq_scope = sc.read(f'{root_dir}/seq_scope/adata_raw.h5ad')

Notice, that instead of loading the already filtered and analysed
dataset like for Visium, here we load the raw dataset (pre-processed by
spacemake) which is unfiltered.

Then we use the
:func:``spacemake.spatial.he_integration.match_he_aggregated_img``
function to match the H&E data with Seq-scope:

.. code:: ipython3

    from spacemake.spatial.he_integration import match_he_aggregated_img
    
    he, matched_he = match_he_aggregated_img(adata_seq_scope,
                                    './wt_4X_1.jpg', bw_threshold=200, box_size=0.6)
    
    cv2.imwrite('matched_he_seq_scope.png', he)




.. parsed-literal::

    True



After we matched the H&E image, we attach it with the same function
:func:`spacemake.spatial.he_integration.attach_he_adata` as for
visium. However here we have to specify that the raw dataset was
matched, and also the matched H&E should not be pushed by one spot
diameter.

.. code:: ipython3

    adata_seq_scope_attached = attach_he_adata(adata_seq_scope.copy(), matched_he,
                             push_by_spot_diameter=False, raw_matched=True)
    
    sc.pl.spatial(adata2, color='total_counts')



.. image:: output_17_0.png
   :width: 1180px

.. image:: output_17_1.png
   :width: 1180px


To see better how the match looks like, we plot only spots with at least
700 UMIs.

.. code:: ipython3

    sc.pl.spatial(adata2[adata2.obs.total_counts > 700,], color='total_counts')
