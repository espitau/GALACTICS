import sys, os
os.environ['TF_CPP_MIN_LOG_LEVEL'] = '3'
import threading
import numpy as np
import tensorflow as tf
import matplotlib.pyplot as plt
tf.logging.set_verbosity(tf.logging.FATAL)

sigma   = 271
verbose = 1

LIM      = 50
# What is the right value... seems to depend on the size..
LIM_NORM = 20
################################################################
#                   Hyperparameters of descent
################################################################

# Testing here
b      = 20
eps    = 0.05
#eps    = 0.015
maxmou = 0.0001
#maxmou = 0.01
maxmou *=sigma**2




# Reload from std input
b = int(sys.argv[2])
eps = float(sys.argv[3])
maxmou = float(sys.argv[4])
################################################################
#                       Load and prepare IV
################################################################

print("------- Loading IV ----------")
fic = np.load(sys.argv[1])
A       = fic['A']
y       = fic['y']
m,n     = A.shape[0],A.shape[1]
print (n,m)
eps     /= np.linalg.norm(fic['s'],ord=1)
signe   =- np.int(np.sign(np.dot(fic['z'],fic['s'])))


################################################################
#                       Setup tensorflow
################################################################

# Working variables
s = tf.Variable(np.float32(fic['z']),name="s")
a = tf.placeholder(tf.float32,[None,n],name="a")
k = tf.placeholder(tf.float32,[None],name="k");

prod    = tf.abs(tf.einsum('ij,j->i',a,s))
proba   = (1+tf.exp(-tf.sqrt(prod**2+maxmou)/sigma/sigma))/2-1e-5
erreur  = -k*tf.log(1-proba) -tf.log(proba)
total   = tf.reduce_mean(erreur)
# Objective function
regul   = total+eps*tf.norm(s,ord=1)
train   = tf.train.AdagradOptimizer(0.07).minimize(regul)
init    = tf.global_variables_initializer()

np.set_printoptions(precision=2,suppress=True)

def calc(j=1):
    sess.run(train,feed_dict={a:np.float32(A2[j]),k:np.float32(y2[j])});

################################################################
#                       Multithreded descent
################################################################
with tf.Session() as sess:

    sess.run(init)
    A2 = np.vsplit(A,b)
    y2 = np.split(y,b)
    dist  = []

    for i in range(LIM):
        # Fix the current estimation to integers if the descent is too
        # slow (avoid to have a divergence if there is not enough samples
        if i%LIM_NORM==0 and i>1:
            curs = sess.run(s)
            for j in range(n):
                fact=1 if j<n/2 else 2
                if abs(round(curs[j]/fact)-curs[j]/fact)<(0.25 if n<1000 else 0.35)/2*fact:
                    curs[j]=round(curs[j]/fact)*fact
            print(curs+signe*fic['s'])
            s.load(curs,sess)
            print(sess.run(total,feed_dict={a:np.float32(A2[0]),k:np.float32(y2[0])}))

        # Set the step to each thread
        thr = []
        for j in range(len(A2)):
            thr.append(threading.Thread(target=calc,kwargs={'j':j}))
        for t in thr: t.start()
        for t in thr: t.join()
        p = sess.run(total,feed_dict={a:np.float32(A2[0]),k:np.float32(y2[0])})
        curs    = sess.run(s)
        cdist   = np.linalg.norm(curs+signe*fic['s'])**2
        dist.append(cdist)
        print (i,"\t\t", cdist)

    np.set_printoptions(precision=2,suppress=True)

    # Retrieve candidate
    rz = sess.run(s)
    f = open("res_des","wb")
    np.savez(f,s=fic['s'], rz=rz);
    f.close()
