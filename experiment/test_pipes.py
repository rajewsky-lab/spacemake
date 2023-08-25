import sys
import os
import multiprocessing as mp
import pyximport; pyximport.install()
import fast_loop
import fcntl
F_SETPIPE_SZ = 1031  # Linux 2.6.35+
F_GETPIPE_SZ = 1032  # Linux 2.6.35+

def open_fifo(fifo, mode='rb+'):
    try:
        fifo_fd = open(fifo, mode)
        print("Pipe size            : "+str(fcntl.fcntl(fifo_fd, F_GETPIPE_SZ)))
        fcntl.fcntl(fifo_fd, F_SETPIPE_SZ, 1000000)
        print("Pipe (modified) size : "+str(fcntl.fcntl(fifo_fd, F_GETPIPE_SZ)))
        return fifo_fd
    except Exception as e:
        print("Unable to create fifo, error: "+str(e))

fcntl.fcntl(sys.stdin, F_SETPIPE_SZ, 1000000)

def sum_up_from_fifo(fname):
    print(f"fname={fname}")
    s = 0
    fifo = open_fifo(fname, 'rt')
    for line in fifo:
        s += 1
    
    print(s)
    print(f"worker {fname}: done")


def make_named_pipe(name):
    if os.path.exists(name):
        os.remove(name)
    
    os.mkfifo(name)
    
fnames = ['fifo1', 'fifo2', 'fifo3', 'fifo4']

processes = []
for name in fnames:
    make_named_pipe(name)
    p = mp.Process(target=sum_up_from_fifo, args=(name,))
    p.start()
    processes.append(p)

#make_named_pipe('igzip_input')
fast_loop.distribute('/dev/stdin', fnames, chunk_size=4)
# fast_loop.distribute(sys.stdin, ['fifo1', 'fifo2'], rec_size=4, chunk_size=1)

# for i, line in enumerate(sys.stdin):
#     x = (i // 10000) % 2
#     fifos[x].write(line)

#fifo1_w.close()
#fifo2_w.close()

for p in processes:
    p.join()

print("main: done")

