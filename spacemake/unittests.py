import unittest
import os
import os.path
import sys
import subprocess
import yaml
import pandas as pd

dropseq_tools = "/data/rajewsky/shared_bins/Drop-seq_tools-2.4.0"
base_dir = os.path.join(os.path.dirname(__file__), "..")

test_species_data = [
    (
        "test_hsa",
        "genome",
        f"{base_dir}/test_data/test_genome.fa.gz",
        f"{base_dir}/test_data/test_genome.gtf.gz",
    ),
    (
        "test_hsa",
        "rRNA",
        f"{base_dir}/test_data/rRNA_hsa.fa.gz",
        f"",
    ),
    (
        "test_hsa",
        "miRNA",
        f"{base_dir}/test_data/mirgenedb.hsa.mature.fa.gz",
        f"{base_dir}/test_data/mirgenedb.hsa.mature.gtf.gz",
    ),
]
test_project_data = [
    (
        "test_hsa",
        "test",
        "test_01",
        f"{base_dir}/test_data/test_reads.R1.fastq.gz",
        f"{base_dir}/test_data/test_reads.R2.fastq.gz",
        "--map_strategy='STAR:genome:final'",
    ),
    # (
    #     "test_hsa",
    #     "test",
    #     "test_02",
    #     f"{base_dir}/test_data/test_reads.R1.fastq.gz",
    #     f"{base_dir}/test_data/test_reads.R2.fastq.gz",
    #     "--map_strategy='bowtie2:rRNA->bowtie2:miRNA->STAR:genome:final'",
    # ),
]


class SpaceMakeCmdlineTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        os.system("mkdir -p _tests")
        os.chdir("_tests")

    @classmethod
    def tearDownClass(cls):
        # os.chdir("..")
        # os.system("rm -rf _tests")
        pass

    def run_spacemake(self, *args, expect_fail=False, **kwargs):
        print(f"running spacemake: cwd='{os.getcwd()}' args={args} kw={kwargs}")
        p = subprocess.run(*args, shell=True, capture_output=True, **kwargs)
        print(f"\nrun_spacemake() results:\nstdout={p.stdout}\nstderr={p.stderr}")
        if expect_fail:
            self.assertFalse(p.stdout.endswith(b"SUCCESS!\n"))
        else:
            self.assertEqual(p.returncode, 0)
            self.assertTrue(p.stdout.endswith(b"SUCCESS!\n"))
            self.assertEqual(len(p.stderr), 0)

        return p

    def load_config(self):
        with open("config.yaml") as yamlfile:
            y = yaml.safe_load(yamlfile.read())

        return y

    def load_project_df(self):
        return pd.read_csv("project_df.csv")

    def add_species(self, name, ref, seq, ann, **kw):
        p = self.run_spacemake(
            f"spacemake config add_species"
            f" --name={name}"
            f" --reference={ref}"
            f" --sequence={seq}"
            f" --annotation={ann}",
            **kw,
        )
        return self.load_config()

    def add_sample(self, species, pid, sid, r1, r2, options, **kw):
        p = self.run_spacemake(
            f"spacemake projects add_sample"
            f" --species={species}"
            f" --project_id={pid}"
            f" --sample_id={sid}"
            f" --R1={r1}"
            f" --R2={r2} {options}",
            **kw,
        )
        return self.load_project_df()

    def test_0_init(self):
        p = self.run_spacemake(f"spacemake init --dropseq_tools={dropseq_tools}")
        self.assertTrue(os.access("config.yaml", os.R_OK))

    def test_1_add_species(self):
        for name, ref, seq, ann in test_species_data:
            y = self.add_species(name, ref, seq, ann)

            self.assertTrue("species" in y)
            self.assertTrue(name in y["species"])
            self.assertTrue(ref in y["species"][name])
            self.assertEqual(y["species"][name][ref]["sequence"], seq)
            self.assertEqual(y["species"][name][ref]["annotation"], ann)

            # expect failure if trying to add same species, reference combination again
            y2 = self.add_species(name, ref, seq, ann, expect_fail=True)

            # expect unchanged config.yaml
            y1_str = yaml.dump(y)
            y2_str = yaml.dump(y)
            self.assertEqual(y1_str, y2_str)

    def test_2_add_project(self):
        for species, pid, sid, r1, r2, options in test_project_data:
            df = self.add_sample(species, pid, sid, r1, r2, options)
            self.assertTrue(os.access("project_df.csv", os.R_OK))
            x = df.set_index("sample_id").loc[sid]
            self.assertEqual(x.species, species)
            self.assertEqual(x.project_id, pid)
            self.assertEqual(x.R1, str([r1]))
            self.assertEqual(x.R2, str([r2]))

            # expect failure if trying to add the same sample again
            df2 = self.add_sample(species, pid, sid, r1, r2, options, expect_fail=True)
            # expect unchanged project_df
            self.assertTrue(df.equals(df2))

    def test_3_run(self):
        self.run_spacemake("spacemake run --cores=8")

        # test correct BAM content
        # test correct DGE content


if __name__ == "__main__":
    unittest.main(verbosity=2)
