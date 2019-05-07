import pylab
import numpy as np
import csv
from scipy.signal import savgol_filter

pylab.rcParams['mathtext.fontset'] = 'stix'
pylab.rcParams['font.family'] = 'STIXGeneral'
pylab.rcParams['font.size'] = 16
params = {'legend.fontsize': 16}
pylab.rcParams.update(params)

# Some linestyle params
l1=(0, (1, 1))
l2=(0, (5, 1))
l3=(0, (3, 1, 1, 1))

# Physical constants
mp = 1.672e-27
AMU = 39.948
mi_to_me = AMU*1836
mi = mp*AMU
B0 = 0.1
cs = 4896.759617436305
B0 = 0.1 # Background magnetic field
eV = 1.602e-19

Lambda = np.log(np.sqrt(mi_to_me/(2*np.pi)))

txtExp = np.loadtxt('Exp_Data/ground/tElc.txt')
x = txtExp[:,0]
tGndB = txtExp[:,1]
vGnd = -np.gradient(tGndB*Lambda,x)/B0
gGnd = np.gradient(vGnd,x)

txtExp = np.loadtxt('Exp_Data/m40/tElc.txt')
tMinusB = txtExp[:,1]
vMinus = -np.gradient(tMinusB*Lambda,x)/B0
gMinus = np.gradient(vMinus,x)

txtExp = np.loadtxt('Exp_Data/p10/tElc.txt')
tPlusB = txtExp[:,1]
vPlus = -np.gradient(tPlusB*Lambda,x)/B0
gPlus = np.gradient(vPlus,x)

txtExp = np.loadtxt('Exp_Data/ground/phi.txt')
xPhi = txtExp[:,0]
phi_gnd = txtExp[:,1]
vGnd = -np.gradient(phi_gnd,xPhi)/B0
gGnd = np.abs(np.gradient(vGnd,xPhi))
gGnd = savgol_filter(gGnd, 51, 3)

txtExp = np.genfromtxt('Exp_Data/m40/phi.txt')
x = txtExp[:,0]
phi_m40 = txtExp[:,1]
vMinus = -np.gradient(phi_m40,xPhi)/B0
gMinus = np.abs(np.gradient(vMinus,xPhi))
gMinus = savgol_filter(gMinus, 51, 3)

txtExp = np.genfromtxt('Exp_Data/p10/phi.txt')
x = txtExp[:,0]
phi_p10 = txtExp[:,1]
vPlus = -np.gradient(phi_p10,xPhi)/B0
gPlus = np.abs(np.gradient(vPlus,xPhi))
gPlus = savgol_filter(gPlus, 51, 3)

# xstart = int((0.8-0.6)*len(x))
# xend = int((1.4-0.6)*len(x))
# x = x[xstart:xend]
# phi_gnd = phi_gnd[xstart:xend]
# phi_m25 = phi_m25[xstart:xend]
# phi_p10 = phi_p10[xstart:xend]

#Compare phi profiles
fig1 = pylab.figure(figsize=(12,4))
fig1.suptitle('Experiment',y=1,x=.52)
ax1 = fig1.add_subplot(1,2,1)
ax1.axvspan(0.86, 1.06, facecolor='lightgray', alpha=0.5)

ax1.plot(x,phi_gnd,label='$V_b = 0$ V',color='C9',linewidth=2,linestyle=l1)
ax1.plot(x,phi_p10,label='$V_b = +10$ V',color='C3',linewidth=2,linestyle=l3)
ax1.plot(x,phi_m40,label='$V_b = -40$ V',color='C4',linewidth=2,linestyle=l2)
#pylab.title('Gkeyll plasma potential \n $L_c$ = 40 m')
#ax1.get_xaxis().set_ticks([])
#ax1.xaxis.set_major_formatter(pylab.NullFormatter())
ax1.set_xlabel('$R$ (m)')
ax1.set_ylabel('$\phi$ (V)')
ax1.legend()
ax1.legend(loc=1)
ax1.set_title("(c)", position=(0.1, 0.85))

# ExB flow
ax2 = fig1.add_subplot(1,2,2)
ax2.axvspan(0.86, 1.06, facecolor='lightgray', alpha=0.5)

ax2.plot(x,vGnd/cs,label='$V_b = 0$ V',color='C9',linewidth=2,linestyle=l1)
ax2.plot(x,vPlus/cs,label='$V_b = +10$ V',color='C3',linewidth=2,linestyle=l3)
ax2.plot(x,vMinus/cs,label='$V_b = -40$ V',color='C4',linewidth=2,linestyle=l2)
#pylab.title('Gkeyll plasma potential \n $L_c$ = 40 m')
ax2.set_xlabel('$R$ (m)')
ax2.set_ylabel(r'$V_E/c_s$')
#ax2.legend(loc=4)
ax2.set_title("(d)", position=(0.1, 0.85))
pylab.tight_layout()
pylab.savefig('thesis-images/phi-bias-exp.pdf')
#pylab.show()
pylab.close()

# compare eq profiles
txtExp = np.loadtxt('Exp_Data/ground/nElc.txt')
x = txtExp[:,0]
nGndB = txtExp[:,1]
#nGndT = txtExp[:,2]

txtExp = np.loadtxt('Exp_Data/m40/nElc.txt')
nMinusB = txtExp[:,1]
#nMinusT = txtExp[:,2]

txtExp = np.loadtxt('Exp_Data/p10/nElc.txt')
nPlusB = txtExp[:,1]
#nPlusT = txtExp[:,2]

# Calculate interchange growth rate
csv = np.genfromtxt('Gk_Data/Lc40p-P150/gammaE_1000to1601.csv', delimiter=",")
gPlus = csv[:,1]

csv = np.genfromtxt('Gk_Data/Lc40g-P150/gammaI_1000to1601.csv', delimiter=",")
xSim = csv[:,0]
gIGnd = csv[:,1]

csv = np.genfromtxt('Gk_Data/Lc40m-P150/gammaI_1000to1601.csv', delimiter=",")
gIMinus = csv[:,1]

csv = np.genfromtxt('Gk_Data/Lc40p-P150/gammaI_1000to1601.csv', delimiter=",")
gIPlus = csv[:,1]

fig1 = pylab.figure(figsize=(12,4))
#fig1.suptitle('Interchange growth rate',x=.4)
ax1 = fig1.add_subplot(1,2,1)
ax1.axvspan(0.86, 1.06, facecolor='gray', alpha=0.5)
#ax1.axvline(x=0.86,color='black')
#ax1.axvline(x=1.06,color='black')
ax1.plot(xSim,gIGnd,label='$V_b = 0$ V ',linewidth=2)
ax1.plot(xSim,gIPlus,label='$V_b = +10$ V ',linewidth=2,linestyle='-.')
ax1.plot(xSim,gIMinus,label='$V_b = -40$ V ',linestyle='--')
ax1.set_xlabel('$R$ (m)')
ax1.set_ylabel('$\gamma_I$ (1/s)')
ax1.ticklabel_format(style='sci',scilimits=(0,0),axis='y')
ax1.set_title("(a) Simulation", position=(0.2, 0.85))
ax1.legend()

cs_loc = np.sqrt(eV*tGndB/mi)
Ln = 1/(-np.gradient(nGndB,x)/nGndB)
sign = Ln/np.abs(Ln)
gamIgnd = sign*cs_loc/np.sqrt(x*abs(Ln))
gamIgnd = np.where(gamIgnd < 0, 0, gamIgnd)
gamIgnd = np.nan_to_num(gamIgnd)
gamIgnd = savgol_filter(gamIgnd, 51, 3)
gamIgnd = np.where(gamIgnd < 0, 0, gamIgnd)

cs_loc = np.sqrt(eV*tPlusB/mi)
Ln = 1/(-np.gradient(nPlusB,x)/nPlusB)
sign = Ln/np.abs(Ln)
gamIplus = sign*cs_loc/np.sqrt(x*abs(Ln))
gamIplus = np.where(gamIplus < 0, 0, gamIplus)
gamIplus = np.nan_to_num(gamIplus)
gamIplus = savgol_filter(gamIplus, 51, 3)
gamIplus = np.where(gamIplus < 0, 0, gamIplus)

cs_loc = np.sqrt(eV*tMinusB/mi)
Ln = 1/(-np.gradient(nMinusB,x)/nMinusB)
sign = Ln/np.abs(Ln)
gamIminus = sign*cs_loc/np.sqrt(x*abs(Ln))
gamIminus = np.where(gamIminus < 0, 0, gamIminus)
gamIminus = np.nan_to_num(gamIminus)
gamIminus = savgol_filter(gamIminus, 51, 3)
gamIminus = np.where(gamIminus < 0, 0, gamIminus)

ax1 = fig1.add_subplot(1,2,2)
ax1.axvspan(0.86, 1.06, facecolor='lightgray', alpha=0.5)
ax1.plot(x,gamIgnd,label='$V_b = 0$ V ',color='C9',linewidth=2,linestyle=l1)
ax1.plot(x,gamIplus,label='$V_b = +10$ V ',color='C3',linewidth=2,linestyle=l3)
ax1.plot(x,gamIminus,label='$V_b = -40$ V ',color='C4',linewidth=2,linestyle=l2)
ax1.set_xlabel('$R$ (m)')
ax1.set_ylabel('$\gamma_I$ (1/s)')
ax1.set_title('Experimental interchange growth rates')
ax1.ticklabel_format(style='sci',scilimits=(0,0),axis='y')
ax1.set_title("(b) Experiment", position=(0.2, 0.85))
ax1.legend()

pylab.tight_layout()
pylab.savefig('thesis-images/gammaI_all.pdf')
pylab.show()

# Plot eq profiles

fig1 = pylab.figure(figsize=(12,4))
fig1.suptitle('Experiment',fontsize=22)
ax1 = fig1.add_subplot(1,2,1)
ax1.axvspan(0.86, 1.06, facecolor='gray', alpha=0.5)
ax1.axvline(x=0.86,color='black')
ax1.axvline(x=1.06,color='black')
ax1.plot(x,nGndB,label='$V_b = 0$ V ',linewidth=2,linestyle=l1)
ax1.plot(x,nPlusB,label='$V_b = +10$ V ',linewidth=2,linestyle=l3)
ax1.plot(x,nMinusB,label='$V_b = -40$ V ',linewidth=2,linestyle=l2)
ax1.set_xlabel('$R$ (m)')
ax1.set_title("Density", position=(0.75, 0.05))
ax1.set_ylabel(r'$n$ (m$^{-3}$)')
ax1.legend()
#ax1.legend()

# plot temp profiles
ax2 = fig1.add_subplot(1,2,2)
ax2.axvspan(0.86, 1.06, facecolor='gray', alpha=0.5)
ax2.axvline(x=0.86,color='black')
ax2.axvline(x=1.06,color='black')
ax2.plot(x,tGndB,label='$V_b = 0$ V ',linewidth=2,linestyle=l1)
ax2.plot(x,tPlusB,label='$V_b = +10$ V ',linewidth=2,linestyle=l3)
ax2.plot(x,tMinusB,label='$V_b = -40$ V ',linewidth=2,linestyle=l2)
ax2.set_xlabel('$R$ (m)')
# increase number of yticks
ax2.set_ylabel(r'$T_e$ (eV)')
ax2.set_title("Electron temperature", position=(.75, 0.05))
#ax2.legend()
pylab.tight_layout()
pylab.savefig('nt-exp-ttf.pdf')
pylab.show()
pylab.close()
exit()

# # Plot shear
# fig2 = pylab.figure(figsize=(12,4))
# fig2.suptitle('Experiment')
# ax1 = fig2.add_subplot(1,2,1)
# ax1.axvspan(0.86, 1.06, facecolor='lightgray', alpha=0.5)
# #ax1.plot(xPhi,gGnd,label='$\gamma_E$, $V_b = 0$ V',color='C9',linewidth=2,linestyle=l1)
# #ax1.plot(xPhi,gMinus,label='$\gamma_E$, $V_b = -40$ V',color='C4',linewidth=2,linestyle=l2)
# #ax1.plot(xPhi,gPlus,label='$\gamma_E$, $V_b = +10$ V',color='C3',linewidth=2,linestyle=l3)
# ax1.plot(xPhi,gGnd,label='$V_b = 0$ V',color='C9',linewidth=2,linestyle=l1)
# ax1.plot(xPhi,gPlus,label='$V_b = +10$ V',color='C3',linewidth=2,linestyle=l3)
# ax1.plot(xPhi,gMinus,label='$V_b = -40$ V',color='C4',linewidth=2,linestyle=l2)
# ax1.plot(x,gamIgnd, label='$\gamma_I (V_b=0)$',color='gray',zorder=1)
# ax1.ticklabel_format(style='sci',scilimits=(0,0),axis='y')
# ax1.set_xlabel('$R$ (m)')
# ax1.set_ylabel(r'shear (1/s)')
# ax1.set_xlim(0.75,1.49)
# ax1.legend(loc=1)
# ax1.set_title("(c)", position=(0.1, 0.85))

# #Compare dn profiles
# txtExp = np.loadtxt('Exp_Data/ground/dn.txt')
# x = txtExp[:,0]
# dnGnd = txtExp[:,1]
# gFac = nGndB/np.mean(nGndB)

# txtExp = np.loadtxt('Exp_Data/m40/dn.txt')
# dnMinus = txtExp[:,1]
# mFac = nMinusB/np.mean(nMinusB)

# txtExp = np.loadtxt('Exp_Data/p10/dn.txt')
# dnPlus = txtExp[:,1]
# pFac = nPlusB/np.mean(nPlusB)

# ax2 = fig2.add_subplot(1,2,2)
# ax2.axvspan(0.86, 1.06, facecolor='lightgray', alpha=0.5)
# ax2.plot(x,dnGnd*gFac,label='$V_b = 0$ V',color='C9',linewidth=2,linestyle=l1)
# ax2.plot(x,dnPlus*pFac,label='$V_b=+10$ V',color='C3',linewidth=2,linestyle=l3)
# ax2.plot(x,dnMinus*mFac,label='$V_b = -40$ V',color='C4',linewidth=2,linestyle=l2)
# #pylab.title('Gkeyll plasma potential \n $L_c$ = 40 m')
# ax2.set_xlabel('$R$ (m)')
# ax2.set_ylabel(r'$\frac{\~{n}_{rms}}{\langle n \rangle_\mathrm{vol}}$',rotation=0,fontsize=24,labelpad=18)
# #ax2.legend(loc=1)
# ax2.set_title("(d)", position=(0.1, 0.85))
# ax2.set_xlim(0.75,1.49)
# pylab.tight_layout()
# pylab.savefig('thesis-images/shear-dn-exp.pdf')
# pylab.show()
# pylab.close()

# Plot dT profiles
txtExp = np.loadtxt('Exp_Data/ground/dT.txt')
x = txtExp[:,0]
print(x)
dTGnd = txtExp[:,1]
txtExp = np.loadtxt('Exp_Data/ground/tBaff.txt')
gTmean = txtExp[:,1]
gFac = gTmean/np.mean(gTmean)

txtExp = np.loadtxt('Exp_Data/m40/dT.txt')
dTMinus = txtExp[:,1]
txtExp = np.loadtxt('Exp_Data/m40/tBaff.txt')
mTmean = txtExp[:,1]
mFac = mTmean/np.mean(mTmean)

txtExp = np.loadtxt('Exp_Data/p10/dT.txt')
dTPlus = txtExp[:,1]
txtExp = np.loadtxt('Exp_Data/p10/tBaff.txt')
pTmean = txtExp[:,1]
pFac = pTmean/np.mean(pTmean)

fig2 = pylab.figure(figsize=(12,4))
fig2.suptitle('Experiment',y=1)
ax1 = fig2.add_subplot(1,2,1)
ax1.axvspan(0.86, 1.06, facecolor='lightgray', alpha=0.5)
ax1.plot(x,dTGnd*gFac,label='$V_b = 0$ V',color='C9',linewidth=2,linestyle=l1)
ax1.plot(x,dTPlus*pFac,label='$V_b = +10$ V',color='C3',linewidth=2,linestyle=l3)
ax1.plot(x,dTMinus*mFac,label='$V_b = -40$ V',color='C4',linewidth=2,linestyle=l2)
ax1.set_ylabel(r'$\frac{\~{T}_{e,rms}}{\langle T_e \rangle}$',rotation=0,fontsize=24,labelpad=18)
ax1.set_xlabel('$R$ (m)')
ax1.legend()
ax1.set_title("(c)", position=(0.1, 0.85))

#Compare dPhi profiles
txtExp = np.loadtxt('Exp_Data/ground/dPhi.txt')
#x = txtExp[:,0]
dPhiGnd = txtExp[:,1]

txtExp = np.loadtxt('Exp_Data/m40/dPhi.txt')
dPhiMinus = txtExp[:,1]

txtExp = np.loadtxt('Exp_Data/p10/dPhi.txt')
dPhiPlus = txtExp[:,1]

ax2 = fig2.add_subplot(1,2,2)
ax2.axvspan(0.86, 1.06, facecolor='lightgray', alpha=0.5)
ax2.plot(x,dPhiGnd*gFac,label='$V_b = 0$ V',color='C9',linewidth=2,linestyle=l1)
ax2.plot(x,dPhiPlus*pFac,label='$V_b=+10$ V',color='C3',linewidth=2,linestyle=l3)
ax2.plot(x,dPhiMinus*mFac,label='$V_b = -40$ V',color='C4',linewidth=2,linestyle=l2)
#pylab.title('Gkeyll plasma potential \n $L_c$ = 40 m')
ax2.set_xlabel('$R$ (m)')
ax2.set_ylabel(r'$\frac{\~{\phi}_{rms}}{\langle T_e \rangle}$',rotation=0,fontsize=24,labelpad=18)
#ax2.legend(loc=1)
ax2.set_title("(d)", position=(0.1, 0.85))
pylab.tight_layout()
pylab.savefig('thesis-images/dT-dphi-exp.pdf')
pylab.show()
pylab.close()


# Plot flux profiles
txtExp = np.loadtxt('Exp_Data/ground/nFlux.txt')
x = txtExp[:,0]
nFluxGnd = txtExp[:,1]

txtExp = np.loadtxt('Exp_Data/m40/nFlux.txt')
nFluxMinus = txtExp[:,1]

txtExp = np.loadtxt('Exp_Data/p10/nFlux.txt')
nFluxPlus = txtExp[:,1]

fig2 = pylab.figure(figsize=(12,4))
fig2.suptitle('Experiment',y=1)
ax1 = fig2.add_subplot(1,2,1)
ax1.axvspan(0.86, 1.06, facecolor='lightgray', alpha=0.5)
ax1.plot(x,nFluxGnd,label='$V_b = 0$ V',color='C9',linewidth=2,linestyle=l1)
ax1.plot(x,nFluxPlus,label='$V_b = +10$ V',color='C3',linewidth=2,linestyle=l3)
ax1.plot(x,nFluxMinus,label='$V_b = -40$ V',color='C4',linewidth=2,linestyle=l2)
ax1.set_ylabel(r'$\tilde{\Gamma}_n/\langle n \rangle$  (m/s)')
ax1.set_xlabel('$R$ (m)')
ax1.legend(loc=1)
ax1.set_title("(c)", position=(0.1, 0.85))

#Compare dPhi profiles
txtExp = np.loadtxt('Exp_Data/ground/tFlux.txt')
#x = txtExp[:,0]
tFluxGnd = txtExp[:,1]

txtExp = np.loadtxt('Exp_Data/m40/tFlux.txt')
tFluxMinus = txtExp[:,1]

txtExp = np.loadtxt('Exp_Data/p10/tFlux.txt')
tFluxPlus = txtExp[:,1]

ax2 = fig2.add_subplot(1,2,2)
ax2.axvspan(0.86, 1.06, facecolor='lightgray', alpha=0.5)
ax2.plot(x,tFluxGnd,label='$V_b = 0$ V',color='C9',linewidth=2,linestyle=l1)
ax2.plot(x,tFluxPlus,label='$V_b=+10$ V',color='C3',linewidth=2,linestyle=l3)
ax2.plot(x,tFluxMinus,label='$V_b = -40$ V',color='C4',linewidth=2,linestyle=l2)
#pylab.title('Gkeyll plasma potential \n $L_c$ = 40 m')
ax2.set_xlabel('$R$ (m)')
ax2.set_ylabel(r'$\tilde{\Gamma}_T/(\langle T_e \rangle \langle n \rangle)$  (m/s)')
#ax2.legend(loc=1)
ax2.set_title("(d)", position=(0.1, 0.85))
pylab.tight_layout()
pylab.savefig('thesis-images/fluxes-exp.pdf')
pylab.show()
pylab.close()
exit()

# Plot lrad
txtExp = np.loadtxt('Exp_Data/ground/cc-dx-1000to1600.csv')
x = txtExp[:,0]
lradGnd = txtExp[:,1]

txtExp = np.loadtxt('Exp_Data/m40/cc-dx-1000to1600.csv')
lradMinus = txtExp[:,1]

txtExp = np.loadtxt('Exp_Data/p10/cc-dx-1000to1600.csv')
lradPlus = txtExp[:,1]

pylab.figure(figsize=(8,4))
pylab.axvspan(0.86, 1.06, facecolor='lightgray', alpha=0.5)
pylab.plot(x,lradGnd,label='$V_b = 0$ V',color='C9',linewidth=2,linestyle=l1)
pylab.plot(x,lradPlus,label='$V_b = +10$ V',color='C3',linewidth=2,linestyle=l3)
pylab.plot(x,lradMinus,label='$V_b = -40$ V',color='C4',linewidth=2,linestyle=l2)
ax = pylab.axes()
ax.set_ylim(0,0.105)
pylab.xlabel('$R$ (m)')
pylab.ylabel(r'$L_{rad}$ (m)')
pylab.legend()
pylab.tight_layout()
pylab.savefig('thesis-images/lrad_Lc40_sim.pdf')
pylab.show()
pylab.close()
