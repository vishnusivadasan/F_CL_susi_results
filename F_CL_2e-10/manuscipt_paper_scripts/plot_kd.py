from __future__ import division
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

phivals = np.array([20,25,30,35,40,45,48,50,52,53,54,55,56])
paths = ["rhor1_phi{}_yd1.0/result/nu_t.dat".format(i) for i in phivals]
def get_and_plot_data(path):
	db = pd.read_csv(path)
	nu_rel = np.mean(db[db.columns[4]][100:])
	std = np.std(db[db.columns[4]][100:])
	return [nu_rel,std]
def kd(phi,phimax):
	return (1-phi/phimax)**(-2.5*phimax)
nu_rel_vals = [get_and_plot_data(i)[0] for i in paths]
std_vals = [get_and_plot_data(i)[1] for i in paths]
#for phimaxval in [0.63,0.635,0.64,0.645,0.65,0.655,0.66,0.665,0.67]:
phimaxval = 0.647
plt.figure(dpi=300)
plt.errorbar(phivals/100.0, nu_rel_vals,std_vals,label="Simulation",fmt="o",elinewidth=0.5,capthick=0.5,capsize=5)
plt.plot(phivals/100.0, [kd(i/100.0,phimaxval) for i in phivals],"g",label="Krieger Dougerty",linewidth=3,)
plt.legend(loc="upper left")
plt.xlabel("$\phi$",fontsize=16)
plt.ylabel("$\eta_{r}$",fontsize=16)
#plt.title(str(phimaxval))
plt.savefig("kd_compare",dpi=300)
