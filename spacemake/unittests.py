#!/usr/bin/env python
import unittest
import os
import os.path
import sys
import subprocess
import yaml
import pandas as pd

dropseq_tools = "/data/rajewsky/shared_bins/Drop-seq_tools-2.4.0"
base_dir = os.path.join(os.path.dirname(__file__), "..")

# set to true to track code coverage of the tests here using 'coverage'
code_coverage = False

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
    (
        "test_hsa",
        "test",
        "test_02",
        f"{base_dir}/test_data/test_reads.R1.fastq.gz",
        f"{base_dir}/test_data/test_reads.R2.fastq.gz",
        "--map_strategy='bowtie2:rRNA->bowtie2:miRNA->STAR:genome:final'",
    ),
    (
        "test_hsa",
        "test",
        "test_03_nofinal",
        f"{base_dir}/test_data/test_reads.R1.fastq.gz",
        f"{base_dir}/test_data/test_reads.R2.fastq.gz",
        "--map_strategy='bowtie2:rRNA->bowtie2:miRNA->STAR:genome'",
    ),
]


def which_spacemake():
    p = subprocess.run("which spacemake", shell=True, capture_output=True)
    return p.stdout.decode("ascii").rstrip()


spacemake_cmd = which_spacemake()

if code_coverage:
    # gathering coverage information
    spacemake_cmd = f"coverage run -p {spacemake_cmd}"


def shell_bam_to_md5(path):
    "prints BAM content as SAM, sorts lexicographically and computes md5 hash"
    p = subprocess.run(
        f"samtools view {path} | sort | md5sum", shell=True, capture_output=True
    )
    md5 = p.stdout.split()[0]
    return (path, md5.decode("ascii"))


def gather_bam_hashes(path):
    results = {}
    for root, dirs, files in os.walk(path):
        for fname in files:
            if fname.endswith("bam"):
                bpath, md5 = shell_bam_to_md5(os.path.join(root, fname))
                results[bpath] = md5

    return results


def print_bam_hashes(results):
    for bpath, md5 in sorted(results.items()):
        print(f"{bpath}\t{md5}")


def load_bam_hashes(path):
    results = {}
    with open(path) as hfile:
        for line in hfile:
            bpath, md5 = line.rstrip().split("\t")
            results[bpath] = md5

    return results


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

    def run_spacemake(
        self,
        *args,
        expect_fail=False,
        check_returncode=True,
        check_success=True,
        check_stderr=True,
        **kwargs,
    ):
        print(
            f"running spacemake: cwd='{os.getcwd()}' args={args} kw={kwargs}", end=" "
        )
        p = subprocess.run(*args, shell=True, capture_output=True, **kwargs)
        print(f"returncode={p.returncode}")
        # print(f"\nrun_spacemake() results:\nstdout={p.stdout}\nstderr={p.stderr}")
        with open("run_spacemake.out.log", "ab") as out:
            if expect_fail:
                self.assertFalse(p.stdout.endswith(b"SUCCESS!\n"))
            else:
                out.write(bytes(f"RUN {args} \n", "ascii"))
                out.write(b"\nSTDOUT\n")
                out.write(p.stdout)
                out.write(b"\nSTDERR\n")
                out.write(p.stderr)
                out.write(bytes(f"RETURNCODE {p.returncode}\n", "ascii"))
                if check_returncode:
                    self.assertEqual(p.returncode, 0)

                if check_success:
                    self.assertTrue(p.stdout.endswith(b"SUCCESS!\n"))

                if check_stderr:
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
            f"{spacemake_cmd} config add_species"
            f" --name={name}"
            f" --reference={ref}"
            f" --sequence={seq}"
            f" --annotation={ann}",
            **kw,
        )
        return self.load_config()

    def add_genome_old(self, name, seq, ann, **kw):
        p = self.run_spacemake(
            f"{spacemake_cmd} config add_species"
            f" --name={name}"
            f" --genome={seq}"
            f" --annotation={ann}",
            **kw,
        )
        return self.load_config()

    def add_sample(self, species, pid, sid, r1, r2, options, **kw):
        p = self.run_spacemake(
            f"{spacemake_cmd} projects add_sample"
            f" --species={species}"
            f" --project_id={pid}"
            f" --sample_id={sid}"
            f" --R1={r1}"
            f" --R2={r2} {options}",
            check_stderr=False,
            **kw,
        )
        return self.load_project_df()

    def test_0_init(self):
        self.run_spacemake(f"{spacemake_cmd} init --dropseq_tools={dropseq_tools}", check_stderr=False)
        self.assertTrue(os.access("config.yaml", os.R_OK))

    def test_1_add_species(self):
        for name, ref, seq, ann in test_species_data:
            if name == "genome":
                # test backward-compatible --genome option
                y = self.add_genome_old(name, seq, ann, check_stderr=False)
            else:
                y = self.add_species(name, ref, seq, ann, check_stderr=False)

            self.assertTrue("species" in y)
            self.assertTrue(name in y["species"])
            self.assertTrue(ref in y["species"][name])
            self.assertEqual(y["species"][name][ref]["sequence"], seq)
            self.assertEqual(y["species"][name][ref]["annotation"], ann)

            # expect failure if trying to add same species, reference combination again
            y2 = self.add_species(name, ref, seq, ann, expect_fail=True, check_stderr=False)

            # expect unchanged config.yaml
            y1_str = yaml.dump(y)
            y2_str = yaml.dump(y2)
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

    def test_6_update_sample(self):
        self.run_spacemake(
            f"{spacemake_cmd} projects update_sample "
            "--project_id=test --sample_id=test_01 "
            "--map_strategy='bowtie2:rRNA->STAR:genome'",
            check_stderr=False
        )
        self.run_spacemake(
            f"{spacemake_cmd} projects update_sample "
            "--project_id=test --sample_id=test_01 "
            "--map_strategy='STAR:genome'",
            check_stderr=False
        )

    def test_99_run(self):
        self.run_spacemake(
            f"{spacemake_cmd} run --cores=8", check_returncode=False, check_stderr=False
        )

#    def test_4_bamcheck(self):
#        # test correct BAM content
#        expect = load_bam_hashes("../test_data/test_bam_md5.txt")
#        for bpath, md5 in sorted(gather_bam_hashes(".").items()):
#            if bpath in expect:
#                print(f"checking '{bpath}'")
#                self.assertEqual(md5, expect[bpath])
#            else:
#                print(f"missing reference checksum for '{bpath}'- skipping test")
#
#        # TODO: test correct DGE content


if __name__ == "__main__":
    ## run this line once, together with output redirect to create
    ## reference md5 hashes from a run you deem correct
    # print_bam_hashes(gather_bam_hashes("."))
    unittest.main(verbosity=2)
