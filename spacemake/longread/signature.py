import os
import logging
from collections import defaultdict

"""
Small helper class to load longread signature definitions (see docs/tutorials/longreads)
and make them accessible to the various cmdline tools.
"""


class SignatureDB:
    def __init__(self, **kw):
        self.lkup = {}
        self.fields = sorted(kw.keys())
        for f in self.fields:
            self.lkup[f] = kw[f]

    @classmethod
    def from_YAML(cls, fname="samples.yaml"):
        import yaml

        logger = logging.getLogger("spacemake.longread.SampleDB.from_YAML")
        logger.info(f"reading longread signature definitions from '{fname}'")

        groups = yaml.load(open(fname), Loader=yaml.SafeLoader)
        grp = groups["signatures"]
        default = grp[groups["default"]]

        field_lkups = {}
        for name, d in groups["signatures"].items():
            # print(f"name={name} d={d}")
            for f in d.keys():
                if f not in field_lkups:
                    field_lkups[f] = defaultdict(lambda: default[f])

                field_lkups[f][name] = d[f]

        logger.info(
            f"found {len(grp)} signature definitions:"
            f"{sorted(groups['signatures'].keys())}."
        )
        return cls(**field_lkups)

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
