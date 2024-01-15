import argparse
import pandas as pd

"""
Global Coordinate System Generator for NovaSeq S4 Flow Cell

This Python script is designed to create a global coordinate system 
for a NovaSeq S4 flow cell. 

It generates a DataFrame with puck names and their corresponding global 
(x, y, z) coordinates and saves it to a CSV file.

Usage:
    python create_openst_coordinate_system.py --output <output_file> [options]

Example:
    python create_openst_coordinate_system.py \
        --output output.csv \
        --format-string fc_1_L{lane}{side_letter}_tile_{side_number}{column}{row:02d} \
        --x-offset 33739 \
        --y-offset 36282 \
        --swath-offset-odd 0 \
        --swath-offset-even 6201 \
        --rows 78 \
        --columns 6 \
        --n_lanes 4 \
        --zero-coded

Author:
    Daniel León-Periñán
"""

def setup_parser(parser):
    parser.add_argument(
        "--output",
        type=str,
        help="where to store the output file with puck names and global (x,y,z) coordinates",
        required=True,
    )

    parser.add_argument(
        "--format-string",
        type=str,
        help="this the format for puck names. There are 4 attributes that can be chosen:"
        + "\{lane\} (int), \{column\} (int), \{row\} (int), \{side_letter\} (str), \{side_number\} (int).\n"
        + "For instance, a valid string format would be: \n"
        + "fc_1_L{lane}{side_letter}_tile_{side_number}{column}{row:02d}\n"
        + "This name must be used, as is, when creating a new sample in spacemake.",
        default="L{lane}{side_letter}_tile_{side_number}{column}{row:02d}",
    )

    parser.add_argument(
        "--x-offset",
        type=int,
        help="the offset in the x axis. Units are important during puck collection generation.",
        default=33739,
    )

    parser.add_argument(
        "--y-offset",
        type=int,
        help="the offset of the y axis. Units are important during puck collection generation.",
        default=36282,
    )

    parser.add_argument(
        "--swath-offset-odd",
        type=int,
        help="the swath offset for odd columns",
        default=0,
    )

    parser.add_argument(
        "--swath-offset-even",
        type=int,
        help="the swath offset for even columns",
        default=6201,
    )

    parser.add_argument(
        "--rows",
        type=int,
        help="number of rows",
        default=78,
    )

    parser.add_argument(
        "--columns",
        type=int,
        help="number of columns",
        default=6,
    )

    parser.add_argument(
        "--n_lanes",
        type=int,
        help="number of lanes",
        default=4,
    )

    parser.add_argument(
        "--zero-coded",
        default=False,
        action="store_true",
        help="whether row and column indices should start at 0, instead of 1",
    )

    return parser


def create_coordinate_system(
    n_lanes: int,
    n_cols: int,
    n_rows: int,
    x_offset: int,
    y_offset: int,
    swath_offsets_odd: int,
    swath_offsets_even: int,
    zero_coded: bool,
    format_string: str,
) -> pd.DataFrame:
    """
    Create a global coordinate system for a NovaSeq S4 flow cell.

    :param n_lanes: Number of lanes in the flow cell.
    :type n_lanes: int
    :param n_cols: Number of columns in the flow cell.
    :type n_cols: int
    :param n_rows: Number of rows in the flow cell.
    :type n_rows: int
    :param x_offset: Offset in the x-axis for coordinate calculations.
    :type x_offset: int
    :param y_offset: Offset in the y-axis for coordinate calculations.
    :type y_offset: int
    :param swath_offsets_odd: Swath offset for odd columns.
    :type swath_offsets_odd: int
    :param swath_offsets_even: Swath offset for even columns.
    :type swath_offsets_even: int
    :param zero_coded: Whether row and column indices should start at 0, instead of 1.
    :type zero_coded: bool
    :param format_string:The format for puck names.
    :type format_string: str
    :returns: DataFrame with puck names and their corresponding global coordinates.
    :rtype: pd.DataFrame
    """

    one_coded_offset = 0 if zero_coded else 1
    swath_offsets = [swath_offsets_even, swath_offsets_odd]
    sides_letter = {1: "a", 2: "b"}
    l = []
    for lane in range(one_coded_offset, n_lanes + one_coded_offset):
        for side in [1, 2]:
            for col in range(n_cols + one_coded_offset):
                for row in range(one_coded_offset, n_rows + one_coded_offset):
                    puck_id = format_string.format(
                        lane=lane,
                        side_letter=sides_letter[side],
                        side_number=side,
                        column=col,
                        row=row,
                    )

                    x_ofs = int(col) * x_offset

                    swath_offset = swath_offsets[int(col) % 2]
                    swath_offset = -swath_offset if side == 1 else swath_offset

                    y_ofs = int(row) * y_offset + swath_offset

                    z_ofs = 0

                    l.append(
                        pd.DataFrame(
                            {
                                "puck_id": [puck_id],
                                "x_offset": [x_ofs],
                                "y_offset": [y_ofs],
                                "z_offset": [z_ofs],
                            }
                        )
                    )

    puck_names_coords = pd.concat(l)

    return puck_names_coords


def cmdline():
    """cmdline."""
    parser = argparse.ArgumentParser(
        allow_abbrev=False,
        description="Global Coordinate System Generator for NovaSeq S4 Flow Cell",
    )
    parser = setup_parser(parser)
    args = parser.parse_args()

    puck_names_coords = create_coordinate_system(
        n_lanes=args.n_lanes,
        n_cols=args.columns,
        n_rows=args.rows,
        x_offset=args.x_offset,
        y_offset=args.y_offset,
        swath_offsets_odd=args.swath_offset_odd,
        swath_offsets_even=args.swath_offset_even,
        zero_coded=args.zero_coded,
        format_string=args.format_string,
    )

    puck_names_coords.to_csv(args.output, index=False)


if __name__ == "__main__":
    cmdline()
