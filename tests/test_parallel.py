import pytest
import os
from spacemake.parallel import *

try:
    from pytest_cov.embed import cleanup_on_sigterm
except ImportError:
    pass
else:
    cleanup_on_sigterm()

spacemake_dir = os.path.dirname(__file__) + '/..'

def simple_worker(Qin, Qout, Qerr, abort_flag, stat_list, **kw):
    with ExceptionLogging("spacemake.test", Qerr=Qerr, exc_flag=abort_flag, set_proc_title=False) as el:
        for ref, sam in AlignedSegmentsFromQueue(Qin, abort_flag):
            print('.')
            put_or_abort(Qout, (ref, sam.qname), abort_flag)
        
def simple_collector(Qin, Qerr, abort_flag, stat_list, **kw):
    for ref, qname in queue_iter(Qin, abort_flag):
        print(ref, qname)
        Qin.task_done()

def test_small_bam():
    parallel_BAM_workflow([f"{spacemake_dir}/test_data/rRNA.bowtie2.bam",], simple_worker, simple_collector, buffer_size=2)
