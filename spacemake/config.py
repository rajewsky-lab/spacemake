import os

# import yaml
import argparse
import re
import logging

from spacemake.errors import *
from spacemake.util import str2bool, assert_file, bool_in_str, message_aggregation
from spacemake.util import check_star_index_compatibility

logger_name = "spacemake.config"


__global_config = None


def get_global_config(root="."):
    global __global_config
    if __global_config is None:
        __global_config = ConfigFile.from_yaml(f"{root}/config.yaml")

    return __global_config

def get_puck_parser(required=True):
    parser = argparse.ArgumentParser(allow_abbrev=False, add_help=False)
    parser.add_argument("--name", help="name of the puck", type=str, required=True)
    parser.add_argument(
        "--width_um", type=float, required=required, help="width of the puck in microns"
    )
    parser.add_argument(
        "--spot_diameter_um",
        type=float,
        required=required,
        help="diameter of the spots in this puck, in microns",
    )
    parser.add_argument(
        "--barcodes",
        type=str,
        required=False,
        help="path to barcode file. if not provided, the --puck_barcode_file variable"
        + " of `spacemake projects add_sample` has to be set",
    )
    parser.add_argument(
        "--coordinate_system",
        type=str,
        required=False,
        help="path to coordinate system file. When specified, spacemake will 'stitch'"
        + " pucks into a single file, with corresponding global coordinates",
    )

    return parser

def get_run_mode_parser(required=True):
    parser = argparse.ArgumentParser(
        allow_abbrev=False,
        formatter_class=argparse.RawTextHelpFormatter,
        description="add/update run_mode parent parser",
        add_help=False,
    )
    parser.add_argument(
        "--name", type=str, help="name of the run_mode to be added", required=True
    )
    parser.add_argument(
        "--parent_run_mode",
        type=str,
        help="Name of the parent run_mode. All run_modes will fall back to 'default'",
    )
    parser.add_argument(
        "--umi_cutoff",
        type=int,
        nargs="+",
        help="umi_cutoff for this run_mode."
        + "the automated analysis will be run with these cutoffs",
    )
    parser.add_argument(
        "--n_beads", type=int, help="number of expected beads for this run_mode"
    )
    parser.add_argument(
        "--clean_dge",
        required=False,
        choices=bool_in_str,
        type=str,
        help="if True, the DGE will be cleaned of barcodes which overlap with primers",
    )
    parser.add_argument(
        "--detect_tissue",
        required=False,
        choices=bool_in_str,
        type=str,
        help="By default only beads having at least umi_cutoff UMI counts are analysed "
        + "during the automated analysis, all other beads are filtered out. If this "
        + "parameter is set, contiguous islands within umi_cutoff passing beads will "
        + "also be included in the analysis",
    )
    # parser.add_argument(
    #     "--polyA_adapter_trimming",
    #     required=False,
    #     choices=bool_in_str,
    #     type=str,
    #     help="if set, reads will have polyA stretches and adapter sequence overlaps trimmed "
    #     + "BEFORE mapping.",
    # )
    parser.add_argument(
        "--count_intronic_reads",
        required=False,
        choices=bool_in_str,
        type=str,
        help="if set, INTRONIC reads will also be countsed (apart from UTR and CDS)",
    )
    parser.add_argument(
        "--count_mm_reads",
        required=False,
        choices=bool_in_str,
        type=str,
        help="if True, multi-mappers will also be counted. For every "
        + "multimapper only reads which have one unique read mapped"
        + "to a CDS or UTR region will be counted",
    )

    parser.add_argument(
        "--mesh_data",
        required=False,
        choices=bool_in_str,
        type=str,
        help="if True, this data will be 'mehsed': a hexagonal structured "
        + "meshgrid will be created, where each new spot will have diameter"
        + " of --mesh_spot_diameter_um micron diameter and the spots will "
        + "be spaced --mesh_spot_distance_um microns apart",
    )
    parser.add_argument(
        "--mesh_type",
        required=False,
        choices=["circle", "hexagon"],
        type=str,
        help="circle: circles with diameter of --mesh_spot_diameter_um will be placed"
        + " in a hex grid, where he distance between any two circles will be "
        + "--mesh_spot_distance_um\n"
        + "hexagon: a mesh of touching hexagons will be created, with their "
        + "centers being --mesh_spot_distance_um apart. This will cover all "
        + "the data without any holes",
    )
    parser.add_argument(
        "--mesh_spot_diameter_um",
        type=float,
        required=False,
        help="diameter of mesh spot, in microns. to create a visium-style "
        + "mesh, use 55um",
    )
    parser.add_argument(
        "--mesh_spot_distance_um",
        type=float,
        required=False,
        help="distance between mesh spots in um. to create a visium-style "
        + "mesh use 100um",
    )
    parser.add_argument(
        "--spatial_barcode_min_matches",
        type=float,
        required=False,
        help="minimum ratio (0, 1] of spatial barcode matches to further consider a puck"
        + "across the rest of the pipeline",
    )

    return parser


def get_species_parser(required=True):
    "a parser that allows to add a reference sequence and annotation, belonging to some species"
    parser = argparse.ArgumentParser(allow_abbrev=False, add_help=False)
    parser.add_argument(
        "--reference",
        help="name of the reference (default=genome)",
        type=str,
        default="genome",
    )
    parser.add_argument("--name", help="name of the species", type=str, required=True)
    parser.add_argument(
        "--sequence",
        help="path to the sequence (.fa) file for the species/reference to be added (e.g. the genome)",
        type=str,
        required=required,
    )
    parser.add_argument(
        "--genome",
        help="[DEPRECATED] path to the genome (.fa) file for the species to be added. --genome=<arg> is a synonym for --reference=genome --sequence=<arg>",
        type=str,
        required=False,
    )

    parser.add_argument(
        "--annotation",
        help="path to the genome annotation (.gtf) file for the species to be added",
        type=str,
        default="",
        required=False,
    )
    parser.add_argument(
        "--STAR_index_dir",
        help="path to STAR index directory",
        type=str,
        required=False,
    )
    parser.add_argument(
        "--BT2_index",
        help="path to BOWTIE2 index",
        type=str,
        required=False,
    )
    parser.add_argument(
        "--BT2_flags",
        help="bt2 mapping arguments for this reference (default=mapping.smk:default_BT2_MAP_FLAGS) ",
        type=str,
        default="",
        required=False,
    )
    parser.add_argument(
        "--STAR_flags",
        help="STAR mapping arguments for this reference (default=mapping.smk:default_STAR_MAP_FLAGS)",
        type=str,
        default="",
        required=False,
    )

    return parser


def get_barcode_flavor_parser(required=True):
    parser = argparse.ArgumentParser(
        allow_abbrev=False, description="add/update barcode-flavor", add_help=False
    )
    parser.add_argument(
        "--name", help="name of the barcode flavor", type=str, required=True
    )
    parser.add_argument(
        "--umi",
        help="structure of UMI, using python's list syntax. Example: to set UMI to "
        + "13-20 NT of Read1, use --umi r1[12:20]. It is also possible to use the first 8nt of "
        + "Read2 as UMI: --umi r2[0:8]",
        type=str,
        required=required,
    )
    parser.add_argument(
        "--cell-barcode",
        help="structure of CELL BARCODE, using python's list syntax. Example: to set"
        + " the cell_barcode to 1-12 nt of Read1, use --cell-barcode r1[0:12]. It is also possible "
        + " to reverse the CELL BARCODE, for instance with r1[0:12][::-1] (reversing the first 12nt of"
        + " Read1, and assigning them as CELL BARCODE).",
        type=str,
        required=required,
    )
    parser.add_argument(
        "--cell_barcode",
        help=argparse.SUPPRESS,
        type=str,
    )

    return parser


def get_variable_action_subparsers(config, parent_parser, variable):
    if variable == "barcode_flavors":
        variable_singular = variable[:-1]
        variable_add_update_parser = get_barcode_flavor_parser
    elif variable == "pucks":
        variable_singular = variable[:-1]
        variable_add_update_parser = get_puck_parser
    elif variable == "run_modes":
        variable_singular = variable[:-1]
        variable_add_update_parser = get_run_mode_parser
    elif variable == "species":
        variable_singular = variable
        variable_add_update_parser = get_species_parser

    command_help = {
        "list": f"list {variable} and their settings",
        "delete": f"delete {variable_singular}",
        "add": f"add a new {variable_singular}",
        "update": f"update an existing {variable_singular}",
    }

    # list command
    list_parser = parent_parser.add_parser(
        f"list_{variable}", description=command_help["list"], help=command_help["list"]
    )
    list_parser.set_defaults(
        func=lambda args: list_variables_cmdline(config, args), variable=variable
    )

    func = lambda args: add_update_delete_variable_cmdline(config, args)

    # delete command
    delete_parser = parent_parser.add_parser(
        f"delete_{variable_singular}",
        description=command_help["delete"],
        help=command_help["delete"],
    )
    delete_parser.add_argument(
        "--name",
        help=f"name of the {variable_singular} to be deleted",
        type=str,
        required=True,
    )
    if variable == "species":
        delete_parser.add_argument(
            "--reference",
            help=f"name of the reference to be deleted (genome, rRNA, ...)",
            type=str,
            required=True,
        )

    delete_parser.set_defaults(func=func, action="delete", variable=variable)

    # add command
    add_parser = parent_parser.add_parser(
        f"add_{variable_singular}",
        parents=[variable_add_update_parser()],
        description=command_help["add"],
        help=command_help["add"],
    )
    add_parser.set_defaults(func=func, action="add", variable=variable)

    # update command
    update_parser = parent_parser.add_parser(
        f"update_{variable_singular}",
        parents=[variable_add_update_parser(False)],
        description=command_help["update"],
        help=command_help["update"],
    )
    update_parser.set_defaults(func=func, action="update", variable=variable)


def setup_config_parser(config, parent_parser_subparsers):
    parser_config = parent_parser_subparsers.add_parser(
        "config", help="configure spacemake"
    )
    parser_config_subparsers = parser_config.add_subparsers(
        help="config sub-command help"
    )

    for variable in ConfigFile.main_variables_pl2sg.keys():
        get_variable_action_subparsers(config, parser_config_subparsers, variable)

    return parser_config


@message_aggregation(logger_name)
def add_update_delete_variable_cmdline(config, args):
    # set the name and delete from dictionary
    name = args["name"]
    variable = args["variable"]
    action = args["action"]

    config.assert_main_variable(variable)

    # remove the args from the dict
    del args["action"]
    del args["variable"]
    del args["name"]

    if action == "add":
        func = config.add_variable
    elif action == "update":
        func = config.update_variable
    elif action == "delete":
        func = config.delete_variable

    var_variables = func(variable, name, **args)
    # print and dump config file
    config.logger.info(yaml.dump(var_variables, sort_keys=False))
    config.dump()


@message_aggregation(logger_name)
def list_variables_cmdline(config, args):
    variable = args["variable"]
    del args["variable"]
    import yaml
    config.logger.info(f"Listing {variable}")
    config.logger.info(yaml.dump(config.variables[variable]))


class ConfigMainVariable:
    def __init__(self, name, **kwargs):
        self.name = name
        self.variables = {}

        # loop over kwargs
        for variable_key, variable in kwargs.items():
            if variable_key not in self.variable_types:
                raise UnrecognisedConfigVariable(
                    f"{variable_key}", list(self.variable_types.keys())
                )
            elif self.variable_types[variable_key] == "int_list":
                self.variables[variable_key] = [int(x) for x in kwargs[variable_key]]
            else:
                new_type = self.variable_types[variable_key]
                self.variables[variable_key] = new_type(kwargs[variable_key])

    def __str__(self):
        class_name = self.__class__.__name__
        return f"{class_name}: {self.name}. variables: {self.variables}"

    def update(self, other):
        self.variables.update(other.variables)


class RunMode(ConfigMainVariable):
    variable_types = {
        "parent_run_mode": str,
        "n_beads": int,
        "umi_cutoff": "int_list",
        "clean_dge": bool,
        "detect_tissue": bool,
        "polyA_adapter_trimming": bool,
        "count_mm_reads": bool,
        "count_intronic_reads": bool,
        "mesh_data": bool,
        "mesh_type": str,
        "mesh_spot_diameter_um": int,
        "mesh_spot_distance_um": int,
        "spatial_barcode_min_matches": float,
    }

    def has_parent(self):
        if "parent_run_mode" in self.variables.keys():
            return True
        else:
            return False

    @property
    def parent_name(self):
        if self.has_parent():
            return self.variables["parent_run_mode"]
        else:
            return None


class Puck(ConfigMainVariable):
    variable_types = {"barcodes": str, "spot_diameter_um": float, "width_um": int, "coordinate_system": str}

    @property
    def has_barcodes(self):
        return (
            "barcodes" in self.variables
            and self.variables["barcodes"]
            and self.variables["barcodes"] != "None"
        )
    
    @property
    def has_coordinate_system(self):
        return (
            "coordinate_system" in self.variables
            and self.variables["coordinate_system"]
            and self.variables["coordinate_system"] != "None"
        )


class ConfigFile:
    initial_config_path = os.path.join(
        os.path.dirname(__file__), "data/config/config.yaml"
    )

    main_variables_pl2sg = {
        "pucks": "puck",
        "barcode_flavors": "barcode_flavor",
        "adapters": "adapter",
        "adapter_flavors": "adapter_flavor",
        "run_modes": "run_mode",
        "species": "species",
        "quant": "quant",
    }

    main_variables_sg2pl = {value: key for key, value in main_variables_pl2sg.items()}

    main_variable_sg2type = {
        "puck": str,
        "barcode_flavor": str,
        "adapter": str,
        "adapter_flavor": str,
        "run_mode": "str_list",
        "species": str,
        "map_strategy": str,
    }

    def __init__(self):
        self.variables = {
            "root_dir": ".",
            "temp_dir": "/tmp",
            "species": {},
            "barcode_flavors": {},
            "adapter_flavors": {},
            "adapters": {},
            "run_modes": {},
            "pucks": {},
            "quant": {},
        }
        self.file_path = "config.yaml"
        self.logger = logging.getLogger(logger_name)

    @classmethod
    def from_yaml(cls, file_path="config.yaml"):
        cf = cls()
        import yaml

        config_yaml_variables = None
        with open(file_path, "r") as f:
            config_yaml_variables = yaml.load(f, Loader=yaml.FullLoader)

        if config_yaml_variables is not None:
            cf.variables.update(config_yaml_variables)

        cf.file_path = file_path

        if file_path != cf.initial_config_path:
            initial_config = ConfigFile.from_yaml(cf.initial_config_path)

            # correct variables to ensure backward compatibility
            cf.correct()
            # check which variables do not exist, if they dont,
            # copy them from initial config
            for main_variable in cf.main_variables_pl2sg:
                # update new main_variables
                if main_variable not in cf.variables:
                    cf.variables[main_variable] = initial_config.variables[
                        main_variable
                    ]

            # deduce variables which have 'default' value. this is to ensure spacemake
            # always runs w/o errors downstream: ie when barcode flavor, run_mode or puck
            # is set to default
            cf.vars_with_default = [
                key
                for key, value in initial_config.variables.items()
                if "default" in value
            ]

            for var_with_default in cf.vars_with_default:
                default_val = initial_config.variables[var_with_default]["default"]
                if "default" not in cf.variables[var_with_default]:
                    cf.variables[var_with_default]["default"] = default_val
                else:
                    # update default run mode with missing values
                    default_val.update(cf.variables[var_with_default]["default"])
                    cf.variables[var_with_default]["default"] = default_val

        cf.expand_strings()
        return cf

    def expand_strings(self, **kw):
        kw['spacemake_dir'] = os.path.abspath(os.path.join(os.path.dirname(__file__), "../"))
        # print(f"expanding strings in config with context={kw}")
        import collections.abc as abc
        def recurse(d):
            if type(d) == str:
                try:
                    return d.format(**kw)
                except KeyError:
                    return d

            if isinstance(d, abc.Mapping):
                for key, value in d.items():
                    d[key] = recurse(value)

            elif isinstance(d, abc.Collection):
                d = [recurse(x) for x in d]

            return d

        self.variables = recurse(self.variables)

    def assert_main_variable(self, variable):
        if variable not in self.main_variables_pl2sg.keys():
            raise UnrecognisedConfigVariable(
                variable, list(self.main_variables_pl2sg.keys())
            )

    def correct(self):
        # ensures backward compatibility
        if "pucks" not in self.variables:
            self.variables["pucks"] = {}
            if "puck_data" in self.variables:
                if "pucks" in self.variables["puck_data"]:
                    self.variables["pucks"] = self.variables["puck_data"]["pucks"]
                    del self.variables["puck_data"]["pucks"]

        if "barcode_flavors" not in self.variables:
            self.variables["barcode_flavors"] = {}
            if "knowledge" in self.variables:
                if "barcode_flavor" in self.variables["knowledge"]:
                    self.variables["barcode_flavors"] = self.variables["knowledge"][
                        "barcode_flavor"
                    ]

        if "species" not in self.variables:
            # get all the species and save them in the right place
            # if species is empty, create a species dictionary
            self.variables["species"] = {}

            # if "knowledge" in self.variables:
            #     # TODO: check if this still applies and make compatible with new
            #     # two-layer species->reference model
            #     # extract all annotation info, if exists
            #     for species in self.variables["knowledge"].get("annotations", {}):
            #         if species not in self.variables["species"]:
            #             self.variables["species"][species] = {}

            #         self.variables["species"][species]["annotation"] = self.variables[
            #             "knowledge"
            #         ]["annotations"][species]

            #     for species in self.variables["knowledge"].get("genomes", {}):
            #         if species not in self.variables["species"]:
            #             self.variables["species"][species] = {}

            #         self.variables["species"][species]["genome"] = self.variables[
            #             "knowledge"
            #         ]["genomes"][species]

            #     for species in self.variables["knowledge"].get("rRNA_genomes", {}):
            #         if species not in self.variables["species"]:
            #             self.variables["species"][species] = {}

            #         self.variables["species"][species]["rRNA_genome"] = self.variables[
            #             "knowledge"
            #         ]["rRNA_genomes"][species]

        for name, species_data in self.variables["species"].items():
            if "annotation" in species_data:
                # we have a pre-bowtie2-support config.yaml with following structur
                # species:
                #   human:
                #     genome: ....
                #     annotation: ....
                self.logger.warning(
                    f"converting old-style config.yaml species section for '{name}' {species_data}..."
                )
                new_data = {}
                new_data["genome"] = {
                    "sequence": species_data["genome"],
                    "annotation": species_data["annotation"],
                }
                if "rRNA_genome" in species_data:
                    new_data["rRNA"] = {
                        "sequence": species_data["rRNA_genome"],
                        "annotation": "",
                    }

                self.variables["species"][name] = new_data

        if "knowledge" in self.variables:
            del self.variables["knowledge"]

        # correct run modes
        for run_mode_name, run_mode_variables in self.variables["run_modes"].items():
            variables = run_mode_variables.copy()
            for var in run_mode_variables:
                if not var in RunMode.variable_types:
                    del variables[var]

            if ("polyA_adapter_trimming" in variables) and (
                variables["polyA_adapter_trimming"] == False
            ):
                import logging

                logger = logging.getLogger("spacemake.config")
                logger.warning(
                    f"WARNING: run_mode {run_mode_name} lists polyA_adapter_trimming=false. This is no longer supported and will be overriden with true"
                )
                variables["polyA_adapter_trimming"] = True

            # print(f"assigning run mode {run_mode_name}: {variables}")
            self.variables["run_modes"][run_mode_name] = variables

    def dump(self):
        import yaml

        self.logger.debug(f"writing config.yaml to '{self.file_path}")
        with open(self.file_path, "w") as fo:
            fo.write(yaml.dump(self.variables))

    def set_file_path(self, file_path):
        self.file_path = file_path

    @property
    def puck_data(self):
        return self.variables["puck_data"]

    def variable_exists(self, variable_name, variable_key):
        return variable_key in self.variables[variable_name]

    def assert_variable(self, variable_name, variable_key):
        self.assert_main_variable(variable_name)
        if not isinstance(variable_key, list):
            variable_key = [variable_key]

        for key in variable_key:
            if not self.variable_exists(variable_name, key):
                variable_name_sg = self.main_variables_pl2sg[variable_name]
                raise ConfigVariableNotFoundError(variable_name_sg, key)

    def delete_variable(self, variable_name, variable_key, **kw):
        self.assert_variable(variable_name, variable_key)

        if variable_name in self.vars_with_default and variable_key == "default":
            raise EmptyConfigVariableError(variable_name)

        variable_data = self.variables[variable_name][variable_key]
        if variable_name == "species":
            ref_name = kw['reference']
            if ref_name in self.variables[variable_name][variable_key]:
                del self.variables[variable_name][variable_key][kw['reference']]
            else:
                logging.warning(f"reference {ref_name} is not registered under species {variable_key}")

            if len(self.variables[variable_name][variable_key]) == 0:
                # this was the last reference entry. 
                # Let's remove the species altogether
                del self.variables[variable_name][variable_key]
        else:
            del self.variables[variable_name][variable_key]

        return variable_data

    def process_run_mode_args(self, name=None, **kwargs):
        # typeset boolean values of run_mode
        default_run_mode = self.get_variable("run_modes", "default")

        for key, value in default_run_mode.items():
            if isinstance(value, bool) and key in kwargs.keys():
                kwargs[key] = str2bool(kwargs[key])

        return kwargs

    def process_barcode_flavor_args(self, **kw):
        bam_tags = "CR:{cell},CB:{cell},MI:{UMI},RG:{assigned}"
        
        # r(1|2) and then string slice
        # this is too restrictive. What about umi = "r1[12:20] + r2[:8]" ? 
        # would be perfectly reasonable 
        # to_match = r"r(1|2)(\[((?=-)-\d+|\d)*\:((?=-)-\d+|\d*)(\:((?=-)-\d+|\d*))*\])+$"

        UMI = kw.get("UMI", None)
        # if umi is not None and re.match(to_match, umi) is None:
        #     raise InvalidBarcodeStructureError("umi", to_match)

        cell_barcode = kw.get("cell", None)
        # if cell_barcode is not None and re.match(to_match, cell_barcode) is None:
        #     raise InvalidBarcodeStructureError("cell", to_match)

        barcode_flavor = {"bam_tags": bam_tags}

        if UMI is not None:
            barcode_flavor["UMI"] = UMI

        if cell_barcode is not None:
            barcode_flavor["cell"] = cell_barcode

        for key, value in kw.items():
            if value is not None:
                barcode_flavor[key] = value

        return barcode_flavor

    def process_adapter_args(self, name=None, seq=None):
        # print(f"processing {name} {seq}")
        if (name is None) or (seq is None):
            raise ValueError(f"need a name and a sequence: name={name} seq={seq}")

        return {name: seq}

    def process_adapter_flavor_args(self, cut_left=[], cut_right=[], name=None, **kw):
        def parse_adapters(cut):
            type_d = {
                "max_errors": float,
                "min_overlap": int,
                "min_base_qual": float,
            }
            from collections import OrderedDict

            out = []
            for c in cut:
                parts = c.split(":")
                adap_name = parts[0]
                if adap_name != "Q" and not adap_name in self.variables["adapters"]:
                    raise KeyError(
                        (
                            f"unknown adapter {adap_name}! We have "
                            f"{','.join(self.variables['adapters'])} "
                            " you can add new adapter with --add_adapter "
                        )
                    )
                d = {}
                if len(parts) > 1:
                    for param_str in parts[1:]:
                        key, value = param_str.split("=")
                        d[key] = type_d.get(key, str)(value)

                out.append({adap_name: d})

            return out

        adapter_flavor = {
            "cut_left": parse_adapters(cut_left),
            "cut_right": parse_adapters(cut_right),
        }

        for key, value in kw.items():
            if value is not None:
                adapter_flavor[key] = value

        # print(f"processed adapter_flavors {adapter_flavor}")
        return {name: adapter_flavor}

    def process_species_args(
        self,
        name=None,
        sequence=None,
        annotation=None,
        reference=None,
        STAR_index_dir=None,
        BT2_index=None,
        BT2_flags=None,
        STAR_flags=None,
    ):
        assert_file(sequence, default_value=None, extension=[".fa", ".fa.gz"])
        if annotation:
            assert_file(annotation, default_value=None, extension=[".gtf", ".gtf.gz"])

        d = dict(
            sequence=sequence,
            annotation=annotation,
        )
        # assert_file(STAR_genome, default_value=None, extension=[".fa"])
        if STAR_index_dir:
            check_star_index_compatibility(STAR_index_dir)
            d["STAR_index_dir"] = STAR_index_dir

        if BT2_index:
            d["BT2_index"] = BT2_index

        if BT2_flags:
            d["BT2_flags"] = BT2_flags

        if STAR_flags:
            d["STAR_flags"] = STAR_flags

        species_refs = self.variables["species"].get(name, {})
        species_refs[reference] = d
        self.variables["species"][name] = species_refs
        return species_refs

    def process_puck_args(self, width_um=None, spot_diameter_um=None, barcodes=None, coordinate_system=None, name=None):
        assert_file(barcodes, default_value=None, extension="all")
        assert_file(coordinate_system, default_value=None, extension="all")

        puck = {}
        if width_um is not None:
            puck["width_um"] = float(width_um)

        if spot_diameter_um is not None:
            puck["spot_diameter_um"] = float(spot_diameter_um)

        if barcodes is not None:
            puck["barcodes"] = barcodes

        if coordinate_system is not None:
            puck["coordinate_system"] = coordinate_system


        return puck

    def process_variable_args(self, variable, **kwargs):
        if variable == "barcode_flavors":
            return self.process_barcode_flavor_args(**kwargs)
        elif variable == "adapters":
            return self.process_adapter_args(**kwargs)
        elif variable == "adapter_flavors":
            return self.process_adapter_flavor_args(**kwargs)
        elif variable == "run_modes":
            return self.process_run_mode_args(**kwargs)
        elif variable == "pucks":
            return self.process_puck_args(**kwargs)
        elif variable == "species":
            return self.process_species_args(**kwargs)
        else:
            ValueError(f"Invalid variable: {variable}")

    def add_variable(self, variable, name, **kwargs):
        kwargs["name"] = name
        # It's not very clear that a cmdline arg
        # --name is absolutely REQUIRED and its value has to map somehow onto
        # an internal function name
        # @TAMAS: can you help?
        # print(f"add_variable() called with variable={variable} name={name} kw={kwargs}")

        if variable == "species":
            # for the species command, collision check is on the reference name, not the species name
            if "genome" in kwargs:
                # deprecated cmdline option --genome ... was used. Translate to
                # --sequence ... --reference=genome
                kwargs["sequence"] = kwargs["genome"]
                kwargs["reference"] = "genome"

            ref = kwargs["reference"]
            if ref in self.variables[variable].get(name, {}):
                raise DuplicateConfigVariableError(variable, f"{name}.{ref}")
            else:
                values = self.process_variable_args(variable, **kwargs)
                self.variables[variable][name] = values
        elif variable == "adapters":
            if not self.variable_exists(variable, name):
                values = self.process_variable_args(variable, **kwargs)
                self.variables[variable].update(values)
            else:
                raise DuplicateConfigVariableError(variable, name)
        else:
            if not self.variable_exists(variable, name):
                values = self.process_variable_args(variable, **kwargs)
                self.variables[variable][name] = values
            else:
                if variable in ["run_modes", "pucks", "barcode_flavors"]:
                    # drop the last s
                    variable = variable[:-1]

                raise DuplicateConfigVariableError(variable, name)

        return kwargs

    def update_variable(self, variable, name, **kwargs):
        if self.variable_exists(variable, name):
            # print(f"updating {variable} {name} {kwargs}")
            values = self.process_variable_args(variable, name=name, **kwargs)
            if variable == "adapters":
                self.variables[variable].update(values)
            elif variable == "adapter_flavors":
                raise NotImplementedError(
                    (
                        "Updating ofadapter flavors is not implemented! "
                        "please delete the flavor and then add aggain."
                    )
                )
            else:
                self.variables[variable][name].update(values)
            variable_data = self.variables[variable][name]
        else:
            if variable in ["run_modes", "pucks", "barcode_flavors"]:
                # drop the last s
                variable = variable[:-1]

            raise ConfigVariableNotFoundError(variable, name)

        return variable_data

    def get_variable(self, variable, name):
        if not self.variable_exists(variable, name):
            raise ConfigVariableNotFoundError(variable, name)
        else:
            return self.variables[variable][name]

    def get_run_mode(self, name):
        # first load the default values
        rm = RunMode(name, **self.get_variable("run_modes", name))

        default_rm = RunMode("default", **self.get_variable("run_modes", "default"))

        if rm.has_parent():
            parent_rm = self.get_run_mode(rm.parent_name)
            parent_rm.update(rm)
            rm = parent_rm

        # update with self
        default_rm.update(rm)
        rm = default_rm

        return rm

    def get_puck(self, name, return_empty=False):
        try:
            return Puck(name, **self.get_variable("pucks", name))
        except ConfigVariableNotFoundError as e:
            if not return_empty:
                raise
            else:
                return Puck(name)

    def get_barcode_flavor(self, flavor):
        return self.get_variable("barcode_flavors", flavor)

    def get_adapter_flavor(
        self, flavor, add_adapter_sequences=["cut_left", "cut_right"]
    ):
        adapter_sequences = self.variables["adapters"]
        af = self.get_variable("adapter_flavors", flavor)
        # populate each adapter clip definition in this flavor with the full sequence of the adapter
        for section, entries in af.items():
            if section in add_adapter_sequences:
                for (name, d_adap) in entries.items():
                    if name == "Q":
                        continue
                    d_adap["seq"] = adapter_sequences.get(
                        name,
                        f"'{name}' NOT FOUND in 'adapters:' section of config.yaml",
                    )

        return af


def get_run_mode_parser(required=True):
    parser = argparse.ArgumentParser(
        allow_abbrev=False,
        formatter_class=argparse.RawTextHelpFormatter,
        description="add/update run_mode parent parser",
        add_help=False,
    )
    parser.add_argument(
        "--name", type=str, help="name of the run_mode to be added", required=True
    )
    parser.add_argument(
        "--parent_run_mode",
        type=str,
        help="Name of the parent run_mode. All run_modes will fall back to 'default'",
    )
    parser.add_argument(
        "--umi_cutoff",
        type=int,
        nargs="+",
        help="umi_cutoff for this run_mode."
        + "the automated analysis will be run with these cutoffs",
    )
    parser.add_argument(
        "--n_beads", type=int, help="number of expected beads for this run_mode"
    )
    parser.add_argument(
        "--clean_dge",
        required=False,
        choices=bool_in_str,
        type=str,
        help="if True, the DGE will be cleaned of barcodes which overlap with primers",
    )
    parser.add_argument(
        "--detect_tissue",
        required=False,
        choices=bool_in_str,
        type=str,
        help="By default only beads having at least umi_cutoff UMI counts are analysed "
        + "during the automated analysis, all other beads are filtered out. If this "
        + "parameter is set, contiguous islands within umi_cutoff passing beads will "
        + "also be included in the analysis",
    )
    parser.add_argument(
        "--polyA_adapter_trimming",
        required=False,
        choices=bool_in_str,
        type=str,
        help="if set, reads will have polyA stretches and adapter sequence overlaps trimmed "
        + "BEFORE mapping.",
    )
    parser.add_argument(
        "--count_intronic_reads",
        required=False,
        choices=bool_in_str,
        type=str,
        help="if set, INTRONIC reads will also be countsed (apart from UTR and CDS)",
    )
    parser.add_argument(
        "--count_mm_reads",
        required=False,
        choices=bool_in_str,
        type=str,
        help="if True, multi-mappers will also be counted. For every "
        + "multimapper only reads which have one unique read mapped"
        + "to a CDS or UTR region will be counted",
    )

    parser.add_argument(
        "--mesh_data",
        required=False,
        choices=bool_in_str,
        type=str,
        help="if True, this data will be 'mehsed': a hexagonal structured "
        + "meshgrid will be created, where each new spot will have diameter"
        + " of --mesh_spot_diameter_um micron diameter and the spots will "
        + "be spaced --mesh_spot_distance_um microns apart",
    )
    parser.add_argument(
        "--mesh_type",
        required=False,
        choices=["circle", "hexagon"],
        type=str,
        help="circle: circles with diameter of --mesh_spot_diameter_um will be placed"
        + " in a hex grid, where he distance between any two circles will be "
        + "--mesh_spot_distance_um\n"
        + "hexagon: a mesh of touching hexagons will be created, with their "
        + "centers being --mesh_spot_distance_um apart. This will cover all "
        + "the data without any holes",
    )
    parser.add_argument(
        "--mesh_spot_diameter_um",
        type=float,
        required=False,
        help="diameter of mesh spot, in microns. to create a visium-style "
        + "mesh, use 55um",
    )
    parser.add_argument(
        "--mesh_spot_distance_um",
        type=float,
        required=False,
        help="distance between mesh spots in um. to create a visium-style "
        + "mesh use 100um",
    )

    return parser


def get_species_parser(required=True):
    "a parser that allows to add a reference sequence and annotation, belonging to some species"
    parser = argparse.ArgumentParser(allow_abbrev=False, add_help=False)
    parser.add_argument(
        "--reference",
        help="name of the reference (default=genome)",
        type=str,
        default="genome",
    )
    parser.add_argument("--name", help="name of the species", type=str, required=True)
    parser.add_argument(
        "--sequence",
        help="path to the sequence (.fa) file for the species/reference to be added (e.g. the genome)",
        type=str,
        required=required,
    )
    parser.add_argument(
        "--genome",
        help="[DEPRECATED] path to the genome (.fa) file for the species to be added. --genome=<arg> is a synonym for --reference=genome --sequence=<arg>",
        type=str,
        required=False,
    )

    parser.add_argument(
        "--annotation",
        help="path to the genome annotation (.gtf) file for the species to be added",
        type=str,
        default="",
        required=False,
    )
    parser.add_argument(
        "--STAR_index_dir",
        help="path to STAR index directory",
        type=str,
        required=False,
    )
    parser.add_argument(
        "--BT2_index",
        help="path to BOWTIE2 index",
        type=str,
        required=False,
    )
    parser.add_argument(
        "--BT2_flags",
        help="bt2 mapping arguments for this reference (default=mapping.smk:default_BT2_MAP_FLAGS) ",
        type=str,
        default="",
        required=False,
    )
    parser.add_argument(
        "--STAR_flags",
        help="STAR mapping arguments for this reference (default=mapping.smk:default_STAR_MAP_FLAGS)",
        type=str,
        default="",
        required=False,
    )

    return parser


def get_barcode_flavor_parser(required=True):
    parser = argparse.ArgumentParser(
        allow_abbrev=False, description="add/update barcode_flavor", add_help=False
    )
    parser.add_argument(
        "--name", help="name of the barcode flavor", type=str, required=True
    )
    parser.add_argument(
        "--umi",
        dest="UMI",
        help="structure of UMI, using python's list syntax. Example: to set UMI to "
        + "13-20 NT of Read1, use --UMI r1[12:20]. It is also possible to use the first 8nt of "
        + "Read2 as UMI: --UMI r2[0:8]",
        type=str,
        required=required,
    )
    # parser.add_argument(
    #     "--read1",
    #     help="can be used to modify or replace 'read1' sequence",
    #     type=str,
    # )

    parser.add_argument(
        "--cell-barcode",
        dest="cell",
        help="structure of CELL BARCODE, using python's list syntax. Example: to set"
        + " the cell_barcode to 1-12 nt of Read1, use --cell-barcode r1[0:12]. It is also possible "
        + " to reverse the CELL BARCODE, for instance with r1[0:12][::-1] (reversing the first 12nt of"
        + " Read1, and assigning them as CELL BARCODE).",
        type=str,
        required=required,
    )
    parser.add_argument(
        "--cell_barcode",
        dest="cell",
        help=argparse.SUPPRESS,
        type=str,
        required=False,
    )

    parser.add_argument(
        "--seq",
        help=(
            "if set, allows to modify the sequence, e.g. "
            "--seq='r2[4:]' will trim the first 4 bases. "
            "NOTE: has to be matched with --qual='r2_qual[4:]' !"
        ),
        type=str,
        default=None,
    )
    parser.add_argument(
        "--qual",
        help=(
            "if set, allows to modify the base qualities, e.g. "
            "--qual='r2_qual[4:]' will trim the first 4 quality values"
            "NOTE: has to be matched with --seq='r2[4:]' !"
        ),
        type=str,
        default=None,
    )
    parser.add_argument(
        "--min-qual-trim",
        help=(
            "Especially if you need access to bases 'at the end' of the read, "
            "that may require trimming bases that bcl2fastq already "
            "marked as non-existent (i.e. sequencing adapter). "
            "a good option for that is to clip bases from the 3' end if they "
            "have very low quality scores (default=disabled)"
        ),
        type=int,
        default=None,
    )
    parser.add_argument(
        "--min_qual_trim",
        help=argparse.SUPPRESS,
        type=int,
        default=None,
    )

    return parser


def get_adapter_parser(required=True):
    parser = argparse.ArgumentParser(
        allow_abbrev=False,
        add_help=False,
        description=("add/update adapter sequences"),
    )
    parser.add_argument("--name", help="name of the adapter", type=str, required=True)
    parser.add_argument(
        "--seq",
        help="sequence of the adapter",
        type=str,
        required=required,
    )
    return parser


def get_adapter_flavor_parser(required=True):
    parser = argparse.ArgumentParser(
        allow_abbrev=False,
        add_help=False,
        description="add/update adapter_flavor",
    )
    parser.add_argument(
        "--name",
        help=(
            "name of the adapter flavor."
            "an 'adapter_flavor' describes a recipe for trimming raw reads. "
            "this can entail clipping sequences from either end (such as TSO "
            "and/or polyA tails) as well as quality-score based trimming "
            "of bases on the reads 3' end. \n"
            "NOTE: adapters need to be first added via \n"
            "   add_adapter --name <name> --seq <sequence>\n"
            "before they can be used in an adapter_flavor. To customize the "
            "parameters for determining a match, key-value pairs can be passed "
            "to --cut_left and --cut_right, for example:\n\n"
            "   --cut_left SMART:min_overlap=10:max_errors=0.1 "
            "   --cut_right Q:min_base_qual=30 polyA:min_overlap=3:max_errors=0.25 "
        ),
        type=str,
        required=True,
    )
    parser.add_argument(
        "--cut_left",
        help=(
            "list the adapter sequences that should "
            "be trimmed from the left side of the read,"
            "e.g. --cut_left TSO_SMART SMART"
        ),
        type=str,
        nargs="+",
    )
    parser.add_argument(
        "--cut_right",
        help=(
            "list the adapter sequences that should "
            "be trimmed from the right side of the read,"
            "e.g. --cut_right Q polyG polyA "
            "NOTE: 'Q' means quality-score trimming"
        ),
        type=str,
        nargs="+",
    )
    parser.add_argument(
        "--paired_end",
        help=(
            "if consecutive reads are mates of a paired-end read,"
            "and trimming leads to a read that is shorter than "
            "min_read_length, we need to discard either both reads,"
            "or, instead of discarding, the short read can be replaced"
            " with a single 'N' [EXPERIMENTAL]"
        ),
        type=str,
        choices=["single_end", "replace_N", "discard_both"],
    )
    parser.add_argument(
        "--min_read_length",
        help=(
            "if triming reduces the length of the read below this"
            " threshold, the read should be discarded (default=18)"
        ),
        type=int,
        default=None,
    )

    return parser


def get_quant_parser(required=True):
    parser = argparse.ArgumentParser(
        allow_abbrev=False,
        add_help=False,
        description="add/update counting flavors (quant)",
    )
    parser.add_argument(
        "--name",
        help=(
            "name of the counting flavor."
        ),
        type=str,
        required=True,
    )
    return parser

def get_variable_action_subparsers(parent_parser, variable):

    if variable == "barcode_flavors":
        variable_singular = variable[:-1]
        variable_add_update_parser = get_barcode_flavor_parser
    elif variable == "adapters":
        variable_singular = variable[:-1]
        variable_add_update_parser = get_adapter_parser
    elif variable == "adapter_flavors":
        variable_singular = variable[:-1]
        variable_add_update_parser = get_adapter_flavor_parser
    elif variable == "pucks":
        variable_singular = variable[:-1]
        variable_add_update_parser = get_puck_parser
    elif variable == "run_modes":
        variable_singular = variable[:-1]
        variable_add_update_parser = get_run_mode_parser
    elif variable == "species":
        variable_singular = variable
        variable_add_update_parser = get_species_parser
    elif variable == "quant":
        variable_singular = variable
        variable_add_update_parser = get_quant_parser

    else:
        print(f"suspicious variable {variable}")

    command_help = {
        "list": f"list {variable} and their settings",
        "delete": f"delete {variable_singular}",
        "add": f"add a new {variable_singular}",
        "update": f"update an existing {variable_singular}",
    }

    # list command
    list_parser = parent_parser.add_parser(
        f"list-{variable.replace('_', '-')}",
        description=command_help["list"],
        help=command_help["list"],
    )
    list_parser.set_defaults(func=list_variables_cmdline, variable=variable)
    # snake_case for backward compatibility. remove in a future update
    list_parser_legacy = parent_parser.add_parser(f"list_{variable}",)
    list_parser_legacy.set_defaults(func=list_variables_cmdline, variable=variable)

    func = add_update_delete_variable_cmdline

    # delete command
    delete_parser = parent_parser.add_parser(
        f"delete-{variable_singular.replace('_', '-')}",
        description=command_help["delete"],
        help=command_help["delete"],
    )
    delete_parser_legacy = parent_parser.add_parser(f"delete_{variable_singular}",) 
    delete_parser.add_argument(
        "--name",
        help=f"name of the {variable_singular} to be deleted",
        type=str,
        required=True,
    )
    delete_parser_legacy.add_argument(
        "--name",
        help=f"name of the {variable_singular} to be deleted",
        type=str,
        required=True,
    )
    if variable == "species":
        delete_parser.add_argument(
            "--reference",
            help=f"name of the reference to be deleted (genome, rRNA, ...)",
            type=str,
            required=True,
        )
        delete_parser_legacy.add_argument(
            "--reference",
            help=f"name of the reference to be deleted (genome, rRNA, ...)",
            type=str,
            required=True,
        )
    delete_parser.set_defaults(func=func, action="delete", variable=variable)
    delete_parser_legacy.set_defaults(func=func, action="delete", variable=variable)

    # add command
    add_parser = parent_parser.add_parser(
        f"add-{variable_singular.replace('_', '-')}",
        parents=[variable_add_update_parser()],
        description=command_help["add"],
        help=command_help["add"],
    )
    add_parser_legacy = parent_parser.add_parser(
        f"add_{variable_singular}",
        parents=[variable_add_update_parser()],
    )
    add_parser.set_defaults(func=func, action="add", variable=variable)
    add_parser_legacy.set_defaults(func=func, action="add", variable=variable)

    # update command
    update_parser = parent_parser.add_parser(
        f"update-{variable_singular.replace('_', '-')}",
        parents=[variable_add_update_parser(False)],
        description=command_help["update"],
        help=command_help["update"],
    )
    update_parser_legacy = parent_parser.add_parser(
        f"update_{variable_singular}",
        parents=[variable_add_update_parser(False)],)
    update_parser.set_defaults(func=func, action="update", variable=variable)
    update_parser_legacy.set_defaults(func=func, action="update", variable=variable)


def setup_config_parser(parent_parser_subparsers):
    parser_config = parent_parser_subparsers.add_parser(
        "config", help="configure spacemake"
    )
    parser_config_subparsers = parser_config.add_subparsers(
        help="config sub-command help"
    )

    for variable in ConfigFile.main_variables_pl2sg.keys():
        get_variable_action_subparsers(parser_config_subparsers, variable)

    return parser_config


@message_aggregation(logger_name)
def add_update_delete_variable_cmdline(args):
    # set the name and delete from dictionary
    name = args["name"]
    variable = args["variable"]
    action = args["action"]

    import yaml

    config = get_global_config()
    config.assert_main_variable(variable)

    # remove the args from the dict
    del args["action"]
    del args["variable"]
    del args["name"]

    if action == "add":
        func = config.add_variable
    elif action == "update":
        func = config.update_variable
    elif action == "delete":
        func = config.delete_variable

    var_variables = func(variable, name, **args)
    # print and dump config file
    config.logger.info(yaml.dump(var_variables, sort_keys=False))
    config.dump()


@message_aggregation(logger_name)
def list_variables_cmdline(args):
    variable = args["variable"]
    del args["variable"]
    import yaml

    config = get_global_config()

    config.logger.info(f"Listing {variable}")
    config.logger.info(yaml.dump(config.variables[variable]))
