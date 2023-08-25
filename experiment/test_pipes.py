import sys
import os
import multiprocessing as mp
import pyximport; pyximport.install()
import fast_loop

def sum_up_from_fifo(fname):
    print(f"fname={fname}")
    s = 0
    fifo = open(fname, 'rt')
    for line in fifo:
        s += 1
    
    print(s)
    print(f"worker {fname}: done")


os.mkfifo('fifo1')
os.mkfifo('fifo2')

p1 = mp.Process(target=sum_up_from_fifo, args=('fifo1',))
p1.start()

p2 = mp.Process(target=sum_up_from_fifo, args=('fifo2',))
p2.start()


fifo1_w = open('fifo1', 'wt')
fifo2_w = open('fifo2', 'wt')

fifos = [fifo1_w, fifo2_w]
buffers = [[], []]

fast_loop.distribute(sys.stdin, fifos, rec_size=4, chunk_size=1)

# for i, line in enumerate(sys.stdin):
#     x = (i // 10000) % 2
#     fifos[x].write(line)

fifo1_w.close()
fifo2_w.close()

p1.join()
p2.join()
print("main: done")


