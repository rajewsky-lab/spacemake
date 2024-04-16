import argparse
import pandas as pd
import networkx as nx
from sklearn.neighbors import kneighbors_graph

"""
Generate a graph for the spatial adjacency of tiles/pucks of a 
NovaSeq 6000 S4 flow cell, from global coordinates in a coordinate system

Usage:
    python create_novaseq_S4_puck_adjacency.py \
        --coordinate-system-in <openst_coordinate_system.csv> \
        --puck-adjacency-out <openst_puck_adjacency.edgelist>

Author:
    Daniel León-Periñán
"""

def setup_parser(parser):
    parser.add_argument(
        "--coordinate-system-in",
        type=str,
        help="where to read the file with puck names and global (x,y,z) coordinates from",
        required=True,
    )

    parser.add_argument(
        "--puck-adjacency-out",
        type=str,
        help="where to write the puck adjacency edgelist into",
        required=True,
    )

    return parser


def create_tile_adjacency(coordinate_system: pd.DataFrame):
    """
    Generates the tile adjacency graph as the k-nearest neighbor graph 
    (k=6, excluding self node), where node names are the puck_names. 

    :param coordinate_system: DataFrame containing puck_id, x_offset, and y_offset.
    :type coordinate_system: pd.DataFrame
    :return: Tile adjacency graph.
    :rtype: nx.Graph
    """

    grouped = coordinate_system.groupby(['lane', 'side'])

    final_graph = nx.Graph()

    for (lane, side), group_df in grouped:
        knn_graph = kneighbors_graph(group_df[['x_offset', 'y_offset']],
                                    n_neighbors=6, mode='connectivity', include_self=False)

        tile_adjacency = nx.Graph(knn_graph)
        puck_names = coordinate_system['puck_id'].tolist()
        node_mapping = {i: puck_names[i] for i in range(len(puck_names))}
        tile_adjacency = nx.relabel_nodes(tile_adjacency, node_mapping)
        final_graph = nx.compose(final_graph, tile_adjacency)

    return final_graph


def cmdline():
    """cmdline."""
    parser = argparse.ArgumentParser(
        allow_abbrev=False,
        description="Spatial adjacency of NovaSeq 6000 S4 flow cell",
    )
    parser = setup_parser(parser)
    args = parser.parse_args()

    coordinate_system = pd.read_csv(args.coordinate_system_in)

    puck_names_coords = create_tile_adjacency(
        coordinate_system
    )

    nx.write_edgelist(puck_names_coords, args.puck_adjacency_out)


if __name__ == "__main__":
    cmdline()
