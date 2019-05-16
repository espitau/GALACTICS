import sys
import numpy as np
import subprocess


def dist_to_2Z(x):
    return min([min(abs(x-i), abs(x+i)) for i in [0,2,4]])

########################################################################
#                            Main routine                              #
########################################################################
B = [20,50,100]
EPS = [0.05, 0.1, 0.005]
MAXMOU = [0.0001, 0.001, 0.000001]

if __name__=='__main__':

    if len(sys.argv) < 2:
        sys.stdout.write("Error: please provide a threshold\n")
        sys.exit(-1)
    name = sys.argv[1]
    threshold = float(sys.argv[2])

    estimated_errors = []

    for b in B:
        for eps in EPS:
            for maxmou in MAXMOU:
                ## Launch the descent
                opt = name+" "+str(b)+" "+str(eps)+" "+str(maxmou)
                cmd = "python3 descente.py "+opt
                print ("-- Launching: "+opt)
                process = subprocess.Popen(cmd.split(), stdout=subprocess.PIPE)
                output, error = process.communicate()

                ## Load and make the localization
                fic = np.load("res_des")
                s = fic['s']
                candidate = fic['rz']
                sign = -np.int(np.sign(np.dot(fic['rz'],fic['s'])))

                # Only look at s_2 part of the secret
                s           = s[513:]
                candidate   = candidate[513:]

                # Average of the point-wise distance between the candidate and
                # the lattice (2Z)^n, where the secret lies.
                dist = map(dist_to_2Z, candidate)

                #g = list_plot(dist)
                #save(g,'dom.pdf')

                # Separate the errors according to their distance to 2Z:
                # heuristically the further away from 2Z, the worst the
                # candidate coefficient is.
                for pos, d in enumerate(dist):
                    if d>threshold: estimated_errors.append(pos)
                estimated_errors = list(set(estimated_errors))
                print ("-- Number of possible errors: ", len(estimated_errors))

                # Construct the natural estimator from the candidate
                # by rounding.
                rz = np.array(list(map(round,list(candidate/2))))
                rz = rz.astype(int)*2
                errfinal=s+sign*rz

                # Yield the errors made by the estimator minus the
                # set of errors estiamted by the heuristics.
                errors, undetected = [], []
                for i in range(len(candidate)):
                    if errfinal[i] != 0:
                        errors.append(i)
                        if i not in estimated_errors: undetected.append(i)

                print("-- Errors at positions: ", errors)
                print("-- Undetected errors: ",   undetected, "(",
                        len(undetected), ")")
