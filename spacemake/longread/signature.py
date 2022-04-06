import os
import logging
from collections import OrderedDict, defaultdict
from spacemake.util import rev_comp

"""
Small helper class to load longread signature definitions (see docs/tutorials/longreads)
and make them accessible to the various cmdline tools.
"""

logger = logging.getLogger("spacemake.longread.signature")


class SignatureDB:
    def __init__(self, blocks=OrderedDict(), **kw):
        self.blocks = blocks
        self.lkup = {}
        self.fields = sorted(kw.keys())
        for f in self.fields:
            self.lkup[f] = kw[f]

    @classmethod
    def from_YAML(cls, fname="samples.yaml"):
        import yaml

        logger = logging.getLogger("spacemake.longread.SignatureDB.from_YAML")
        logger.info(f"reading longread signature definitions from '{fname}'")

        groups = yaml.load(open(fname), Loader=yaml.SafeLoader)
        signatures = groups["signatures"]
        default = signatures[groups["default"]]

        # load all building block oligo sequences and their reverse complements
        blocks = OrderedDict()
        for fa_id, seq in groups["blocks"].items():
            blocks[fa_id] = seq
            blocks[fa_id + "_RC"] = rev_comp(seq)

        logger.info(f"load_oligos(): loaded {len(blocks)} sequences from '{fname}'")

        # load the signature definitions and split into separate dictionaries
        field_lkups = {}
        for name, d in signatures.items():
            # print(f"name={name} d={d}")
            for f in d.keys():
                if f not in field_lkups:
                    field_lkups[f] = defaultdict(lambda: default[f])

                field_lkups[f][name] = d[f]

        logger.info(
            f"found {len(signatures)} signature definitions:"
            f"{sorted(signatures.keys())}."
        )
        return cls(blocks, **field_lkups)

    def __getattr__(self, attr):
        return self.lkup[attr]

    def sort_samples(self, samples, signatures):
        """
        Sort samples by the priority assigned in the signature definitions first,
        then lexicographically. Used for overview plots combining multiple longread
        sample results to group samples sharing a signature.
        """
        return sorted(
            zip(samples, signatures), key=lambda x: (self.prio.get(x[1], np.inf), x[0])
        )


def get_signature_db(try_path):
    """
    try to load a YAML file with longread signature definitions from <try_path>.
    If that fails, default to spacemake/data/config/longread.yaml
    """
    if os.access(try_path, os.R_OK):
        cfg = try_path
    else:
        cfg = os.path.join(os.path.dirname(__file__), "../data/config/longread.yaml")

    return SignatureDB.from_YAML(cfg)


def process_intact_signature(complete_signature, prefixes=["P5"], suffixes=["N70X"]):
    complete = complete_signature.split(",")
    while complete and complete[0] in prefixes:
        complete.pop(0)

    while complete and complete[-1] in suffixes:
        complete.pop()

    complete_order = dict(x[::-1] for x in enumerate(complete))
    # print(f"complete={complete}")

    return tuple(complete), complete_order


def digest_signatures(
    sig_counts,
    bead_related="bead_start",
    complete_signature="P5,bead_start,OP1,polyT,N70X",
    prefixes=[
        "P5",
    ],
    suffixes=[
        "N70X",
    ],
):
    bead_counts = defaultdict(int)
    ov_counts = defaultdict(int)
    n_bead_related = 0

    complete, complete_order = process_intact_signature(
        complete_signature, prefixes, suffixes
    )
    complete_set = set(complete)
    found_part_counts = defaultdict(int)

    def describe(found_set):
        missing = complete_set - found_set
        if not missing:
            descr = "complete"
        elif len(missing) < len(found_set):
            descr = f"missing_{','.join(sorted(missing))}"
        else:
            descr = f"only_{','.join(sorted(found_set))}"

        return descr

    def bead_relation(parts):
        search = list(complete)
        at = 0

        try:
            i = parts.index(search[0])  # look for first part, e.g. bead_start
        except ValueError:
            i = 0

        found = []
        for part in parts[i:]:
            # find co-linear matches,
            # ignore extra inserted segments
            # (for now)
            if part in search[at:]:
                found.append(part)
                at = search.index(part)

        found_set = set(found)
        found_tup = tuple(sorted(found_set, key=lambda x: complete_order[x]))

        return describe(found_set), found_tup

    for sig, count in sig_counts.items():
        parts = sig.split(",")
        if bead_related in parts:
            br, found_tup = bead_relation(parts)
            bead_counts[br] += count
            n_bead_related += count

            for i in range(1, len(found_tup) + 1):
                found_part_counts[found_tup[:i]] += count
        else:
            ov_counts[sig] = count

    ov_counts["bead-related"] = n_bead_related
    return ov_counts, bead_counts, found_part_counts, complete
