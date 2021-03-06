import os

from shutil import *  
from numpy import *
from matplotlib import *
import matplotlib.pyplot as plt


#### 10/30/2012:  particle_tools.py overview.  This is a collection of functions that are used by the 
# main S2E_tools operatings, including all code-specific read-write functions.

# Comments on which direction is the head (positive or negative) for the various particle formats:
# Elegant:  +z (s) = TAIL
# fmt1   
# fmt3
# xpx
# xxs
#
#
#


c0=299792458; mue0=2*pi*1e-7
eps0=1.0/(c0*c0*mue0); Z0=c0*mue0;
m0_ele=9.1081e-31; q0_ele=1.60203e-19
E_ele_eV=(m0_ele*c0*c0)/q0_ele


def AdjustAcceleration(in_file,out_file,V0,phi0,V1,phi1,omega):
    AS=astra(in_file)                           #read astra distribution
    n=size(AS[:,0])
    dfi=AS[:,2].copy()*0; dE=AS[:,2].copy()*0   #memory allocation
    dfi[0]=0;                                   #reference particle
    dfi[1:n]=AS[1:n,2]*omega/c0;                #phase offset
    #plot(AS[1:n,2],dfi[1:n],'.'); show()
    dE=V1*cos(phi1+dfi)-V0*cos(phi0+dfi)        #energy change in eV
    #plot(AS[1:n,2],dE[1:n],'.'); show()
    px=AS[:,3];py=AS[:,4];
    pz=AS[:,5].copy();
    pz[1:n]=pz[1:n]+AS[0,5];                    #add reference particle momentum
    E=dfi                                       # change pointer
    E=sqrt(E_ele_eV*E_ele_eV+px*px+py*py+pz*pz)+dE #new energy
    pz=sqrt(E*E-E_ele_eV*E_ele_eV-px*px-py*py)  #new momentum
    pz[1:n]=pz[1:n]-pz[0]                       # substruct the reference particle
    AS[:,5]=pz
    f=open(out_file,'w')
    writeT(f, AS)
    f.close()

def Momentum(energy):
    return E_ele_eV*sqrt((energy/E_ele_eV)**2-1)

def Energy(momentum):
    return E_ele_eV*sqrt((momentum/E_ele_eV)**2+1)

#            Dh=particle_reader(in_file,format,E_ref,extr_part_1,sort_len)
#            N=Dh[0]; xpx=Dh[1];  xxs=Dh[2]; q=Dh[3]

#
# CRP 10/30/2012 Note: This is the main particle reading functions, that checks what 'fmt' is equal to.  
#
#
#
def particle_reader(name,fmt,E_ref,extr_part_1,sort_len,properties):
#
#   output: 0  Np   4  Eref
#           1  xpx  5  qtot
#           2  xxs  6  zav
#           3  q    7  ct
#
#   uses:   astra_to_xpx           xpx_to_xxs           fmt3_to_xxs
#           xxs_to_xpx             sort                 momenta_2
#           fmt1_to_xpx
#
           OUT=[0]; zav=0; ct=0; N_p=0;
           if fmt=='astra':
                print("Reading astra::")
                D=astra_to_xpx(name)
                xpx=D[0]; q=D[1]; q_tot=D[2]; zav=D[3]
                N_p=size(q)
#                plt.plot (xpx[:,5])
#                plt.show()
                if E_ref<=0:
                    E_ref=P_av=MoM1(xpx[:,5],q,0)
                xxs=xpx_to_xxs(xpx,E_ref)
#                plt.plot (xxs[:,5])
#                plt.show()
#
#PP added 04/13/2017
#            
           if fmt=='astragenerator':
                print("Reading astra generator file::")
                D=astragenerator_to_xpx(name)
                xpx=D[0]; q=D[1]; q_tot=D[2]; zav=D[3]
                N_p=size(q)
#                plt.plot (xpx[:,4],xpx[:,5])
#                plt.show()
                if E_ref<=0:
                    E_ref=P_av=MoM1(xpx[:,5],q,0)
                xxs=xpx_to_xxs(xpx,E_ref)
#
#PP added 11/10/2015            
#
# this assumes the coordinate are dumped in warp as follow:
#                 
           if fmt=='warp':
                print("Reading warp::")
                print("***WARNING: currently assumes cathode dump***")
#               D=warp_to_xpx(name)
                D=warp_to_xpx_cath(name)
                xpx=D[0]; q=D[1]; q_tot=D[2]; zav=D[3]
                N_p=size(q)
                # here E_ref is the average momentum
                if E_ref<=0:
#PP need some correction here 01/03/16        
                    print("getting E_ref...")	 
                    P_av=MoM1(xpx[:,5],q,0)
                    E_ref=E_reference_xpx (xpx,q)
                print(E_ref)
                   
                xxs=xpx_to_xxs(xpx,E_ref)
                
                print("longitudinal momentum information:")
                print("mean [eV]=",mean(D[0][:,5]),", max [eV]=",max(D[0][:,5]),\
                     ", min [eV]=",min(D[0][:,5]))
#
#                print "mean=",mean(D[0][:,5]),"max=",max(D[0][:,5]),\
#                     "min=",min(D[0][:,5])
#
#PP added 11/16/2010            
#        
           if fmt=='impactT':
               print("Reading impactT::")
               D=impactT_to_xpx(name)
               xpx=D[0]; q=D[1]; q_tot=D[2]; zav=D[3]
            
# CRP added 10/17/2011
               bunch_charge=eval(properties['bunch_charge'])
               q=q*bunch_charge/len(xpx)
               q_tot=sum(q)
               N_p=size(q)
# CRP Modified 01/30/12, as there is a meaningful deviation between E_ref and E_average for some bunches that led to issues.
# Modified such that E_ref=-1 gives the new behavior, and E_ref
               if E_ref==-1:
                   E_ref=mean(xpx[:,5])
                   print('E_ref=mean(pz)')
               else:
                       if E_ref<=0:
                              E_ref=P_av=MoM1(xpx[:,5],q,0)
                       print('E_ref=old routine')
               xxs=xpx_to_xxs(xpx,E_ref)
#            continue
           if fmt=='FMT1':
               D=fmt1_to_xpx(name)
               xpx=D[0];q=D[1]; q_tot=D[2]; zav=D[3]; ct=D[4]
               N_p=size(q)
# CRP Modified 01/30/12, as there is a meaningful deviation between E_ref and E_average for some bunches that led to issues.
# Modified such that E_ref=-1 gives the new behavior, and E_ref, a simple average of pz.  However, this MUST be constant weight macroparticles.
               if E_ref==-1:
                   E_ref=mean(xpx[:,5])
                   print('E_ref=mean(pz)')
               else:
                   if E_ref<=0:
                        E_ref=P_av=MoM1(xpx[:,5],q,0)
                        print('E_ref=old routine')
               xxs=xpx_to_xxs(xpx,E_ref)
#               continue
           if fmt=='FMT3':
                   D=fmt3_to_xxs(name)
                   xxs=D[0];q=D[1]; q_tot=D[2]; ct=D[4]; E_ref=D[5]
                   if E_ref>0:
                       N_p=size(q)
                       xpx=xxs_to_xpx(xxs,E_ref)
                   else:
                       N_p=0
#                   continue
           if fmt=='CUSTOM1':
                D=Custom1_to_xxs(name)
                xxs=D[0]; q=D[1];  ct=D[2]; E_ref=D[3]; q_tot=D[4]
                print(E_ref)
                if E_ref>0:
                       N_p=size(xxs[:,0])
                       print(N_p, E_ref)
#                       print xxs
                       xpx=xxs_to_xpx(xxs,E_ref)
                else:
                       N_p=0
                       print(N_p, E_ref)
#                   else:
#                       N_p=0
#                continue
           if fmt=='SDDSOUT':
                 os.system('sddsprintout -noLabel -noTitle -column=x -column=xp -column=y -column=yp -column=t -column=p -width=600 $1 ' + name + ' > TempDistributions/temp.sdds')
                 print('dumped temp.sdds')

                 os.system('pwd')
                 D=SDDSout_to_xxs('TempDistributions/temp.sdds')
                 xxs=D[0]; q=D[1];  ct=D[2]; E_ref=D[3] #; q_tot=D[4]
#### Bunch Charge must be specified in the GlueTrack file, as it is for ImpactT.  Yes, it is possible to extract it from the .out file
                 bunch_charge=eval(properties['bunch_charge'])
                 q=q*bunch_charge/len(xxs)
                 q_tot=sum(q)
                 N_p=size(q)
                 print(E_ref, q_tot, std(xxs[:,0]), std(xxs[:,2]), std(xxs[:,4]))
                 if E_ref>0:
                     N_p=size(xxs[:,0])
                     print(N_p, E_ref)
                     xpx=xxs_to_xpx(xxs,E_ref)
                 else:
                     N_p=0
                     print("sddsout=",N_p, E_ref)
#        else:
#             N_p=0
           if N_p>0:
                if extr_part_1:
            # Kills the first particle.
                    xxs=xxs[1:N_p,0:6]
                    xpx=xpx[1:N_p,0:6]
                    q=q[1:N_p]
                    N_p=N_p-1
                if sort_len:
                    isi=sort1(xxs[:,4])
                    xxs=sort2(xxs,isi)
                    xpx=sort2(xpx,isi)
                    q=sort2(q,isi)
                q=abs(q)
                
           OUT=[N_p,xpx,xxs,q,E_ref,q_tot,zav,ct]
           return OUT
        
# CRP 10/30/2012 Note: ASTRA's 9-column format reading.  
def astra(name):
#       0   1   2   3   4   5   6       7       8       9
#       x   y   z   px  py  pz  clock   charge  index   status
#       m           cV/c        ns      nC
    f=open(name,'r');    AS=readT(f);    f.close()
    # extract invalid particles
    ind_0=[]
    for i in range(AS.shape[0]):
        if AS[i,9]<0: ind_0.append(i)
    AS=remove(AS,ind_0)
    return AS

def astragenerator(name):
#       0   1   2   3   4   5   6       7       8       9
#       x   y   z   px  py  pz  clock   charge  index   status
#       m           cV/c        ns      nC
    f=open(name,'r');    AS=readT(f);    f.close()
    # extract invalid particles
    return AS

# CRP 10/30/2012 Note: Conversion from 9-column ASTRA to xpx.  
def astra_to_xpx(name):
    AS=astra(name)
    n=size(AS[:,0])
    AS[1:n,2]+=AS[0,2]; AS[1:n,5]+=AS[0,5]  #add reference particle parameters
    y=AS[:,1].copy();    AS[:,1]=AS[:,3];    z=AS[:,2].copy()
    zav=mean(z)
    AS[:,2]=y;   AS[:,3]=AS[:,4];    AS[:,4]=z-zav
#PP    plot (AS[:,0],AS[:,5])
#PP    show() 
    return [AS[:,0:6],AS[:,7]*1e-9,abs(sum(AS[:,7])*1e-9),zav]

def astragenerator_to_xpx(name):
    AS=astragenerator(name)
    n=size(AS[:,0])
    AS[1:n,2]+=AS[0,2]; AS[1:n,5]+=AS[0,5]  #add reference particle parameters
    y=AS[:,1].copy();    AS[:,1]=AS[:,3];    z=AS[:,6]*1e-9  # take the time 
    z=z*c0 # convert time to z
    zav=mean(z)
    # the astra generator reference particle is set to the average bunch energy
    pav=mean(AS[1:n,5])
    print(AS[0,5])
    print("average energy:", pav)
    AS[0,5]=pav
    AS[:,2]=y;   AS[:,3]=AS[:,4];    AS[:,4]=z-zav
    return [AS[:,0:6],AS[:,7]*1e-9,abs(sum(AS[:,7])*1e-9),zav]

def warp(name):
#PP 11/10/2015 MB CHANGE NAME
#PP read warp output file (dump from the function beamdiag.py)
#
#    phasespace=numpy.vstack([x,bgx, y, bgy, z, bgz, ke, t, w]).T
#       0   1      2      3      4      5   6    7   8      
#       x   bgx    y     bgy     z    gbz   ke   t  w
#       m   []     m      []     m     []
    f=open(name,'r');    AS=readT(f);    f.close()
    return AS

# PP 11/10/2015 Note: Conversion from WARP format (custom one) to xpx.  
#
def warp_to_xpx_cath(name):
#PP 11/17/2015 this record in AS c*time for z assuming we use a Xcrossing 
#   function to save the file
#  numpy.vstack([x,bgx, y, bgy, z, bgz, ke, t, w]).T in particlediag
#
    AS=warp(name)
    n=size(AS[:,0])

#    x   = AS[:,0]
#    bgx = AS[:,1]
#    y   = AS[:,2]
#    bgy = AS[:,3]
#    z   = AS[:,4]
#    bgz = AS[:,5]
#    ke  = AS[:,6]
#    t   = AS[:,7]
#    q   = AS[:,8]
#pp    print (t)
    
# we assume the beam propagate in x (other case will come later)    
# so we do a circular permutation of warp coordinate to match x with z in astra
    AS[:,1] *=E_ele_eV
    AS[:,3] *=E_ele_eV
    AS[:,4] = AS[:,7]*c0
    AS[:,5] *=E_ele_eV
#    print AS[:,5]
    print(mean(AS[:,5]))
    zav=mean(AS[:,4])
    AS[:,4]=AS[:,4]-zav
    AS[:,8]=AS[:,8]
    return [AS[:,0:6],AS[:,8],abs(sum(AS[:,8])),zav]

# PP 11/10/2015 Note: Conversion from WARP format (custom one) to xpx.  
#
def warp_to_xpx(name):
#PP 11/10/2015 MB CHANGE NAME of function warp() like warp_tsnap()??
    AS=warp(name)
    n=size(AS[:,0])
    z=AS[:,2].copy()
    zav=mean(z)
    AS[:,4]=z-zav
    AS[:,8]=AS[:,8]
    return [AS[:,0:6],AS[:,8],abs(sum(AS[:,8])),zav]

# CRP 10/30/2012 Note: Reads the 6-column Impact-T
#
def impactT(name):
#PP 11/15/2010
#PP read impactT standard output file (dumped with -2)
#       0   1      2      3      4      5         
#       x   bgx    y     bgy     z    gbz
#       m   []     m      []     m     []
    f=open(name,'r');    AS=readT(f);    f.close()
    # extract invalid particles
    ind_0=[]
    for i in range(AS.shape[0]):
        if AS[i,5]==0: ind_0.append(i)
    AS=remove(AS,ind_0)
    return AS

# CRP 10/30/2012 Note: Converts Impact-T to xpx
#
def impactT_to_xpx(name):
#PP 11/15/2010
#
    f=open(name,'r');    AS=readT(f);    f.close()
    n=size(AS[:,0])-1
#    AS=impactT(name)
    n=size(AS[:,0])
    z=AS[:,4].copy()
    zav=mean(z)
    AS[:,1]=AS[:,1]*E_ele_eV;   AS[:,3]=AS[:,3]*E_ele_eV;  
    AS[:,5]=AS[:,5]*E_ele_eV;   AS[:,4]=z-zav
    q=zeros((n))
    q+=1.00; q_tot=sum(q) 
    print(std(AS[:,0]), std(AS[:,2]), std(AS[:,4])) 
# CRP modification 10/26/2011
#    q+=1.00 
#    q+=q_tot/n
#    q_tot=sum(q) 
    return [AS,q,q_tot,zav]




def MoM1(A,C,csum):
    csu=csum
    if csum==0:  csu=sum(C)
    return sum(C*A)/csu

# CRP 10/30/2012 Note: Very important function, switch from xpx to xxs.  Very commonly used, and it is important for new functions that work in xpx to copy to xxs as well.
def xpx_to_xxs(xpx,E_ref):
    h=xsysde(E_ref,xpx[:,1],xpx[:,3],xpx[:,5])
    xxs=xpx.copy()
    xxs[:,1]=h[0]   #xs
    xxs[:,3]=h[1]   #ys
    xxs[:,5]=h[2]   #dE
    return xxs

def xsysde(E_ref,px,py,pz):
    E=sqrt(E_ele_eV*E_ele_eV+px*px+py*py+pz*pz)
#    print (E) 
    return [px/pz, py/pz,(E-E_ref)/E_ref]


# CRP 10/30/2012 Note:  Very important function, switch from xxs to xpx.  
# Very commonly used, and it is important for new functions that work in xxs to copy to xpx as well.
def xxs_to_xpx(xxs,E_ref):
    print('in xxs_to_xxp')
    p_z=Momentum(E_ref,xxs[:,5])/sqrt(1.+xxs[:,1]**2+xxs[:,3]**2)
    print(p_z)
#    print 'fact'
    xpx=xxs.copy()
#    print '...', len(fact)
    xpx[:,1]=xxs[:,1]*p_z   #px
    xpx[:,3]=xxs[:,3]*p_z   #py
    xpx[:,5]=p_z   #pz
#PP    plt.plot(xpx[:,4],xpx[:,5])
#PP    plt.show()
#    print size(xpx)  # Temp diagnostic
    return xpx

def Momentum(E_ref,delta):
# 01/03/16, PP: compute reference momentum given the total beam energy 
    E=E_ref*(1.0+delta); gamma=E/E_ele_eV; beta=sqrt(1.0-1.0/(gamma*gamma))
#    print E_ele_eV*gamma*beta
    return E_ele_eV*gamma*beta
    
def E_reference_xpx (xpx,q):
# 01/03/16, PP: compute reference (average) beam energy from 
# xpx coordinate
    paverage2 = MoM1(xpx[:,5],q,0)**2+MoM1(xpx[:,5],q,0)**2+MoM1(xpx[:,5],q,0)**2
    return sqrt(paverage2 + E_ele_eV**2)
    
def sort1(a) : 
    return argsort(a)

def sort2 (a,k): #changes raws due to indexis 'k'
    return take(a,k,0)

def s_to_cur(A,sigma,q0,v):
    Nsigma=3;
    a=min(A)-Nsigma*sigma;    b=max(A)+Nsigma*sigma
    s=0.25*sigma
    N=int(ceil((b-a)/s)); s=(b-a)/N
    B=zeros((N+1,2)); C=zeros(N+1)
    sss=size(B[:,0])
    #print sss,N,s,a,b
    B[:,0]=arange(0,(N+0.5)*s,s)+a
    cA=(A-a)/s;
    I=floor(cA);
    xiA=1+I-cA
    for k in range(0,size(A)):
        i=int(I[k]);
        C[i]=C[i]+xiA[k];
        C[i+1]=C[i+1]+(1-xiA[k]);
    K=floor(Nsigma*sigma/s+0.5)
    G=exp(-0.5*(arange(-K*s,K*s+s,s)/sigma)**2); G=G/sum(G)
    B[:,1]=convolve(C,G,mode=1)
    koef=q0*v/(s*sum(B[:,1]))
    B[:,1]=koef*B[:,1]
    return B

def wake_conv(H,w,wd):
    def Int1(x): return (1-x)*w(x*cdt)
    def Int2(x): return (1-abs(x))*w((x+i)*cdt)
    N=size(H[:,0]); cta=-H[N-1,0]; ctb=-H[0,0]
    cdt=(ctb-cta)/(N-1)
    W=zeros(N)
    nW=min([100,N-1])
    W[0],err=integrate.quad(Int1,0,1)
    W[0]+=wd/cdt
    for i in range(1,nW):
        W[i],err=integrate.quad(Int2,-1,1)
    if N-1>100:
        W[100:N]=w(arange(100,N)*cdt)
    pH=H[:,1];
    pH=fliplr(pH[NewAxis,:])[0]
    U=cdt*convolve(pH,W,mode=2)
    U=fliplr(U[NewAxis,0:N])[0]
    Out=zeros((N,2));
    Out[:,0]=H[:,0];Out[:,1]=U;
    return Out

def MoM(A,C):
    AB=A.copy()
    csu=sum(C)
    aav=sum(C*AB[:,0])/csu
    AB[:,0]-=aav
    bav=sum(C*AB[:,1])/csu
    AB[:,1]-=bav
    a_a=sum(C*AB[:,0]**2)/csu
    a_b=sum(C*AB[:,0]*AB[:,1])/csu
    b_b=sum(C*AB[:,1]**2)/csu
    detw=sqrt(a_a*b_b-a_b**2)
    return -a_b/detw,a_a/detw,b_b/detw,detw

def T(E1,abm1,E2,abm2):
    w=sqrt(E1/E2)
    alpha1=abm1[0];beta1=abm1[1];alpha2=abm2[0];beta2=abm2[1];
    mue=abm2[2]-abm1[2]
    Out=zeros((2,2))
    Out[0,0]=w*sqrt(beta2/beta1)*(cos(mue)+alpha1*sin(mue))
    Out[0,1]=w*sqrt(beta1*beta2)*sin(mue)
    Out[1,0]=w*((alpha1-alpha2)*cos(mue)-(1+alpha1*alpha2)*sin(mue))/sqrt(beta1*beta2)
    Out[1,1]=w*sqrt(beta1/beta2)*(cos(mue)-alpha2*sin(mue))
    return Out

def ShowTwissPar(xxs,w):
    alpha,beta,gamma,epsilon=MoM(xxs,w)
    print('alpha,beta,gamma,geometric epsilon:')
    print(alpha,beta,gamma,epsilon)
    return alpha,beta,gamma,epsilon

def ShowTwissPar2(xxs,w, E):
    alpha,beta,gamma,epsilon=MoM(xxs,w)
    gamma=E/E_ele_eV; beta=sqrt(1.0-1.0/(gamma*gamma))
    print('alpha,beta,gamma,geometric epsilon, normalized epsilon:')
    print(alpha,beta,gamma,epsilon,beta*gamma*epsilon)
    return alpha,beta,gamma,epsilon

#CRP 11/08/2011   Write Twiss to File
# Four lines:  Alpha, Beta, Gamma, Epsilon
def PrintTwissPar(fileName, xxs,w):
    alpha,beta,gamma,epsilon=MoM(xxs,w)
    print('alpha,beta,gamma,epsilon:')
    print(alpha,beta,gamma,epsilon)
    f=open(fileName,'w');
    f.write(str(alpha))
    f.write(" ") 
    f.write(str(beta)) 
    f.write(" ") 
    f.write(str(gamma)) 
    f.write(" ") 
    f.write(str(epsilon)) 
    f.write(" ")
    f.close()
    return alpha,beta,gamma,epsilon

#CRP 10/30/2012 Note:  Writing to CSRtrack's input file, with the first line as the reference particle.  Whether or not this is wanted depends on what parameters your input deck is set to.
def xpx_to_fmt1_withRef(xpx,q,z0,ct,E_ref):
#   t   r1  r2      r3      r4      r5      r6
#   x1  y1  z1      px1     py1     pz1     q1
#   dx2 dy2 dz2     dpx2    dpy2    dpz2    q2
    N=size(xpx[:,0])
    F=zeros((N+2,7)) 

#Modified, to add the reference particle as the second line.
#    F[1,0]=numpy.average(xpx[:,4]);    F[1,1]=numpy.average(xpx[:,0])
#    F[1,2]=numpy.average(xpx[:,2]);    F[1,3]=numpy.average(xpx[:,5])
#    F[1,4]=numpy.average(xpx[:,1]);    F[1,5]=numpy.average(xpx[:,3])+E_ref/511000
    F[1,0]=0;    F[1,1]=0
    F[1,2]=0;    F[1,3]=E_ref
    F[1,4]=0;    F[1,5]=0
    F[1,6]=0
    F[2:N+2,0]=xpx[:,4];    F[2:N+2,1]=xpx[:,0]
    F[2:N+2,2]=xpx[:,2];    F[2:N+2,3]=xpx[:,5]
    F[2:N+2,4]=xpx[:,1];    F[2:N+2,5]=xpx[:,3]
    F[2:N+2,0:6]-=F[1,0:6]
    F[1,0]=z0
    F[1,6]=q[1]
    F[2:N+2,6]=q
    F[1,6]=0
    F[0,0]=ct
    
    return F


# this is the old, non-reference fmt1 code, to be used when using CSRtracks ability to define reference particles in the input deck.
#def xpx_to_fmt1(xpx,q,z0,ct):
#CRP 10/30/2012 Note:  Writing to CSRtrack's input file, without the first line as the reference particle.  Whether or not this is wanted depends on what parameters your input deck is set to.
def xpx_to_fmt1(xpx,q,z0,ct):
#   t   r1  r2      r3      r4      r5      r6
#   x1  y1  z1      px1     py1     pz1     q1
#   dx2 dy2 dz2     dpx2    dpy2    dpz2    q2
    N=size(xpx[:,0])
    F=zeros((N+1,7))
    F[1:N+1,0]=xpx[:,4];    F[1:N+1,1]=xpx[:,0]
    F[1:N+1,2]=xpx[:,2];    F[1:N+1,3]=xpx[:,5]
    F[1:N+1,4]=xpx[:,1];    F[1:N+1,5]=xpx[:,3]
    F[2:N+1,0:6]-=F[1,0:6]
    F[1,0]=z0
    F[1:N+1,6]=q
    F[0,0]=ct
    return F



#CRP 10/30/2012 Note:  Generic write table.  Opens the filename, dumps the matrix to it, then closes the file.
def writetable(name,A,head=''):
    f=open(name,'w');
    if head!='': f.write(head)
    writeT(f,A);
    f.close()


#CRP 11/15/11: ImpactZ expects the first line to be the number of particles, so its write function is slightly different than just writing an X*Y matrix.
def writetableIMPACT(name,A,N):
    f=open(name,'w');
#AL 06142017
    head="  %d\n" % N
    f.write(head)
    writeT(f,A);
    f.close()


#CRP 10/30/2012 Note: Writes ImpactT without the leading number.  This is necessary if you are processing an output file, such as with down sampling, for use with other scripts or analysis tools.
def writetableIMPACT_noNumber(name,A,N):
    f=open(name,'w');
#    head="%d\n" % N
#    f.write(head)
    writeT(f,A);
    f.close()


#CRP 10/30/2012 Note: This is the function for the basic READ operation.  Reads whatever file.
def readtable(name):
    f=open(name,'r');
    A=readT(f);
    f.close()
    return A

#CRP 10/30/2012 Note: Read a CSRtrack input to xpx
def fmt1_to_xpx(name):
#   t   r1  r2      r3      r4      r5      r6
#   x1  y1  z1      px1     py1     pz1     q1
#   dx2 dy2 dz2     dpx2    dpy2    dpz2    q2
    f=open(name,'r');    F1=readT(f);    f.close()
    N=size(F1[:,0])-1
    xpx=zeros((N,6))
    xpx[:,4]=F1[1:N+1,0];xpx[:,0]=F1[1:N+1,1]
    xpx[:,2]=F1[1:N+1,2];xpx[:,5]=F1[1:N+1,3]
    xpx[:,1]=F1[1:N+1,4];xpx[:,3]=F1[1:N+1,5]
    q=zeros((N))
    q=F1[1:N+1,6]
    xpx[1:N+1,0:6]=xpx[1:N+1,0:6]+xpx[0,0:6]
    zav=sum(xpx[:,4])/N
    xpx[:,4]-=zav
    q_tot=sum(q)
    ct=F1[0,0]
    OUT=[xpx,q,q_tot,zav,ct]
    return OUT



# CRP 10/15/11, this was referenced but never defined.  Created from scratch.
# Compared to CSRtrack documentation, and tested for validity.
def fmt3_to_xxs(name):
# 12 columns
# general
# reference
# particle
# primes
# Charges
# forces
# angular momentum.
#   ct gammar          0          0        0        0        0        sigs        sigr        sigz        0        0
#   hr hprr    vr        vpr        0        deltar        q1        fs1        fh1        fv1        Ls1        Lt1
#   hn hprn    vn        vpn        sn        deltan        qn        fsn        fhn        fvn        Lsn        Ltn
    f=open(name,'r');    F1=readT(f);    f.close()
    N=size(F1[:,0])-1
    xxs=zeros((N,6))
    xxs[:,0]=F1[1:N+1,0]
    xxs[:,1]=F1[1:N+1,1]
    xxs[:,2]=F1[1:N+1,2]
    xxs[:,3]=F1[1:N+1,3]
    xxs[:,4]=F1[1:N+1,4]
    xxs[:,5]=F1[1:N+1,5]
    q=zeros((N))
    q=F1[1:N+1,6]
    xxs[1:N+1,0:6]=xxs[1:N+1,0:6]+xxs[0,0:6]
    zav=sum(xxs[1:,4])/N
    xxs[:,4]-=zav
    q_tot=sum(q)
    ct=F1[0,0]
    E_ref=F1[0,1]*511000.
    OUT=[xxs,q,q_tot,zav,ct,E_ref]
    return OUT

#CRP 10/30/2012 Note: Function for writing to ASTRA.
def xpx_to_astra(xpx,q,z0,ct):
#       0   1   2   3   4   5   6       7       8       9
#       x   y   z   px  py  pz  clock   charge  index   status
#       m           eV/c        ns      nC
    N=size(xpx[:,0])
    print(N)
#pp    AS=zeros((N,10))
    AS=zeros((N,10))
    if z0<0:
       # we dump the time coordinate as z/c and set the longitudinal position to 0
       AS[:,0]=xpx[:,0];AS[:,1]=xpx[:,2]
       AS[:,2]=0.0*xpx[:,4];AS[:,3]=xpx[:,1]
       AS[:,4]=xpx[:,3];AS[:,5]=xpx[:,5]
       AS[1:N,2]-=AS[0,2];AS[1:N,5]-=AS[0,5]
       AS[:,7]=q*1e9
       AS[:,6]=xpx[:,4]/c0*1.e9  # the time is in nsec in ASTRA
       AS[:,8]=1;AS[:,9]=-1
       AS[0,2]=0.0
       print("dumping cathode distribution all z set to 0.")
    else:
       AS[:,0]=xpx[:,0];AS[:,1]=xpx[:,2]
       AS[:,2]=xpx[:,4];AS[:,3]=xpx[:,1]
       AS[:,4]=xpx[:,3];AS[:,5]=xpx[:,5]
       AS[1:N,2]-=AS[0,2];AS[1:N,5]-=AS[0,5]
       AS[:,7]=q*1e9
       AS[:,6]=ct;AS[:,8]=1;AS[:,9]=5
       AS[0,2]=z0
    return AS

#CRP 10/30/2012 Note: Converting to ImpactT.
def xpx_to_impactT(xpx,q,z0,ct):
#PP 11/15/2010
#       0   1      2      3      4      5         
#       x   bgx    y     bgy     z    gbz
#       m   []     m      []     m     []
    N=size(xpx[:,0])
#pp    AS=zeros((N,10))
    AS=zeros((N,6))
#PP 4/13/2017
    if z0<0:
      print("!!! cathode dump as z0<0 !!!")
      AS[:,0]=xpx[:,0];AS[:,1]=xpx[:,1]/E_ele_eV
      AS[:,2]=xpx[:,2];AS[:,3]=xpx[:,3]/E_ele_eV
      AS[:,5]=xpx[:,5]/E_ele_eV
      print('total emission time [sec] =', abs(min(xpx[:,4])-max(xpx[:,4]))/c0)
      print('kinetic energy      [eV]  =', mean(sqrt(AS[:,1]**2+AS[:,3]**2+AS[:,5]**2+1)-1.)*E_ele_eV)
      gamma=1.0+sqrt(AS[:,1]**2+AS[:,3]**2+AS[:,5]**2)      
      betaz=AS[:,5]/gamma   
      AS[:,4]=betaz*(xpx[:,4]-max(xpx[:,4]))-1e-32
#AL 06142017
      AS = insert(AS, 0, [0., 0., 0., 0., AS[0,4],AS[0,5]], axis=0)
      return AS
      
    else:
      AS[:,0]=xpx[:,0];AS[:,1]=xpx[:,1]/E_ele_eV
      AS[:,2]=xpx[:,2];AS[:,3]=xpx[:,3]/E_ele_eV
      AS[:,4]=xpx[:,4];AS[:,5]=xpx[:,5]/E_ele_eV
      print(std(AS[:,0]), std(AS[:,2]), std(AS[:,4]))
      return AS


def BCcoefS(x1,x2,x3,R,gamma): #S-type
    fi=arcsin(x1/R)
    y=R*cos(fi)*tan(fi/2.0)
    return BCcoefM(x1,2*x2+y,2*x3,R,6,gamma)


def BCcoefC(x1,x2,x3,R,gamma): #C-type
    return BCcoefM(x1,x2,x3,R,4,gamma)


def BCcoefM(x1,x2,x3,R,m,gamma):
    t0p5=sqrt(1-(x1/R)**2)
    t1p5=t0p5**3;   t2p5=t0p5**5;   t3p5=t0p5**7
    C=zeros(4);D=zeros(3);
    C[0]=m*R*arcsin(x1/R)+2*x2/t0p5+x3
    C[1]=-m*x1/t0p5+m*R*arcsin(x1/R)-2.0*(x2/t1p5)*(x1/R)**2
    C[2]=m*x1**3.0/(2.0*R**2*t1p5)+3.0*(x2/t1p5)*(x1/R)**2+3.0*(x2/t2p5)*(x1/R)**4
    C[3]=-4.0*(x2/t1p5)*(x1/R)**2-9.0*(x2/t2p5)*(x1/R)**4-5.0*(x2/t3p5)*(x1/R)**6
    C[3]=C[3]-2.0*(m/3.0)*x1**3/(R*R*t1p5)-m/2.0*(x1**5)/((R**4)*t2p5)
    D[0]=1;     D[1]=-1.0/(gamma**2-1); D[2]=3.0/(2.0*(gamma**2-1))+3.0/(2.0*(gamma**2-1)**2)
    return (C[0],C[1]+C[0]*D[1],C[2]+C[1]*D[1]+C[0]*D[2],C[3]+C[2]*D[1]+C[1]*D[2])



####################
#CRP 10/30/2012 Note:  Several sum or weighted-sums, but they are here only for consistency?  This is a relic of the original design that I don't understand.
def Sum0(q):
    return cumsum(q)

def Sum1(x,q):
    y=x*q
    return cumsum(y)

def Sum2(x,y,q):
    z=x*y*q
    return cumsum(z)

##########



#CRP 10/30/2012 Note: ????
def Mh(n,M,N):
    return int(     floor(  min([0.5*M,0.5*n,0.5*(N-n)])+0.501  )   )


######### Slice Analysis Functions ############
#CRP 10/30/2012 Note:  Slice functions are called by slice_analysis, in order.

def Slice_Ana_0(xxs,q,E_ref):
    N=size(xxs[:,0])
    X=zeros((N,5))
    X[:,0]=Sum0(q)                      #   Q
    X[:,1]=xxs[:,4]                     #   z
    X[:,2]=xxs_to_xpx(xxs,E_ref)[:,5]   #   pz
    X[:,3]=Sum1(X[:,2],q)               #   sum(pz*q) -weighted mean value
    X[:,4]=Sum1(xxs[:,4],q)             #   sum (z*q)
    return X

def Slice_Ana_1(a,b,q):
    # calculate weighted moments - covariance matrix
    N=a.shape[0]
    Y=zeros((N,5))
    Y[:,0]=Sum1(a,q)
    Y[:,1]=Sum1(b,q)
    Y[:,2]=Sum2(a,a,q)
    Y[:,3]=Sum2(a,b,q)
    Y[:,4]=Sum2(b,b,q)
    return Y
    
def Slice_Ana_2(X,Y,M):
    N=X.shape[0]
    Z=zeros((N,8))
    N=N-1;
    Pzav=X[N,3]/X[N,0]
    gammabetaav=Pzav/E_ele_eV
    gammaav=sqrt(gammabetaav**2+1)
    for nn in range(N-1):
        mh=Mh(nn+1,M,N)
        n1=nn+1-mh
        n2=nn+1+mh
        dq=X[n2,0]-X[n1,0]
        Z[nn,0]=(X[n2,4]-X[n1,4])/dq
        Z[nn,1]=dq/(X[n2,1]-X[n1,1])
        A=(Y[n2,0]-Y[n1,0])/dq
        B=(Y[n2,1]-Y[n1,1])/dq
        AA=(Y[n2,2]-Y[n1,2])/dq-A*A
        AB=(Y[n2,3]-Y[n1,3])/dq-A*B
        BB=(Y[n2,4]-Y[n1,4])/dq-B*B
        Z[nn,2]=A;  Z[nn,3]=B;  Z[nn,4]=AA;
        Z[nn,5]=AB; Z[nn,6]=BB;
        Z[nn,7]=gammaav*sqrt(abs(AA*BB-AB*AB));
    Z[N,0]=Z[N-1,0]=Z[N-2,0]
    return Z

def Slice_Ana_3(X,q,M):
    N=X.shape[0]
    Z=zeros((N,2))
    z=zeros(N)
    N=N-1;
    for nn in range(N-1):
        mh=Mh(nn+1,M,N)
        n1=nn+1-mh
        n2=nn+1+mh
        dq=X[n2,0]-X[n1,0]
        z[nn+1]=(X[n2,3]-X[n1,3])/dq
    z[0]=z[1]; z[N]=z[N-1]
    Z[:,0]=z
    z=X[:,2]-z
    z=Sum2(z,z,q)
    for nn in range(N-1):
        mh=Mh(nn+1,M,N)
        n1=nn+1-mh
        n2=nn+1+mh
        dq=X[n2,0]-X[n1,0]
        Z[nn,1]=sqrt((z[n2]-z[n1])/dq)
    return Z

########################

#CRP 10/30/2012 Note: Calculates the Twiss Parameters of a bunch
def Twiss(X):
    N=X.shape[0]
    Z=zeros((N,2))
    #h=zeros(N)
    h=abs(X[:,4]*X[:,6]-X[:,5]*X[:,5])
    Z[:,0]=where(h>1e-4*X[:,4]*X[:,6],-X[:,5]*1/sqrt(h),0)
    Z[:,1]=where(h>1e-4*X[:,4]*X[:,6],X[:,4]*1/sqrt(h),0)
    return Z



#CRP 10/30/2012 Note:  Filters a distribution.  Method ???    
def Gauss_Filter(lambdaA,z,sigma):
    Nsigma=3;
    za=min(z)-Nsigma*sigma;    zb=max(z)+Nsigma*sigma
    dz=0.25*sigma
    N=int(ceil((zb-za)/dz)); dz=(zb-za)/N
    H=zeros((N+1,2)); HH=zeros(N+1)
    #print sss,N,s,a,b
    H[:,0]=arange(0,(N+1)*dz,dz)+za
    i2=int(floor((z[0]-za)/dz))
    for i in range(1,size(z)):
        i1=i2
        i2=int(floor((z[i]-za)/dz))
        if i1!=i2:
            X=1-(z[i]-H[i2,0])/dz
            HH[i2]=lambdaA[i]*X+lambdaA[i-1]*(1-X)
    K=int(Nsigma*sigma/dz)
    G=zeros(K+1)
    G[0]=1
    su=1
    for k in range(1,K+1):
        G[k]=exp(-0.5*(k/K)**2)
        su=su+G[k]*2
    G=G/su
    for n in range(K,N-K+1):
        for i in range(-K,K+1):
            H[n-i,1]=H[n-i,1]+HH[n]*G[abs(i)]
    return H


#CRP 10/30/2012 Note:  Copies a file.
def copyfiles(src, dst):
    names = os.listdir(src)
    if not(os.path.exists(dst)): os.mkdir(dst)
    for name in names:
        srcname = os.path.join(src, name)
        dstname = os.path.join(dst, name)
        copy2(srcname, dstname)

def _remove(array, indices):
    def_indices = list(range(array.shape[0]))
    n=size(indices)
    if n>1 :
        for i in range(n):
           def_indices.remove(indices[i])
    else:
        def_indices.remove(indices)
    
    return take(array, def_indices)

def remove(array, indices, axis=0):
    n=size(indices)
    if n>0:
        array = asarray(array)
        if axis == 0:
            array = _remove(array, indices)
        else:
            array = swapaxes(array, 0, axis)
            array = _remove(array, indices)
            array = swapaxes(array, 0, axis)
    return array

#CRP 10/30/2012 Note: The basic write table.  This is only applicable for no specific format, basic matrix.
def writeT(file, a):
    """Write a two-dim. NumPy array a in tabular form."""
    print('writing ...')
    if len(a.shape) > 2:
        raise TypeError("a 2D array is required, shape now is "+str(a.shape))
    else:
        if len(a.shape) == 2:
            for i in range(a.shape[0]):
                for j in range(a.shape[1]):
#AL 06142017
                    s=" %12.5e" % a[i,j]
                    file.write(s)
                file.write("\n")
        else:
            for i in range(a.shape[0]):
                s="%2.10g\n" % a[i]
                file.write(s)
                
#CRP 10/30/2012 Note:   Important part of Read functionality.
def readT(file, commentchar='#'):
    """Load a table with numbers into a two-dim. NumPy array."""
    # based on a version by Mario Pernici <Mario.Pernici@mi.infn.it>
    print('reading...', end=' ')
    location = file.tell()
    import re
    commentchar = re.escape(commentchar)
    while 1:  # might call read several times for a file
        line = file.readline()
        if not line: break  # end of file
        elif line.isspace(): break  # blank line
        elif re.match(commentchar, line): continue # treat next line
        else: break

    shape1 = len(line.split())
    if shape1 == 0: return None
    file.seek(location)

    blankline = re.compile('\n\s*\n',  re.M)
    commentline = re.compile('^%s[^\n]*\n' % commentchar, re.M)
    filestr = file.read()
    # remove lines after a blank line
    m = re.search(blankline, filestr)
    if m:
        filestr = filestr[:m.start()+1]
    # skip lines starting with the comment character
    filestr = re.sub(commentline, '', filestr)
    a = [float(x) for x in filestr.split()]
#pp    data = array(a, Float)
    data = array(a)
    data.shape = (int(len(a)/shape1), shape1)
    return data
        

def writevectors(name,*args):
    nj=len(args);
    ni=size(args[0]);
    A=zeros((ni,nj));
    for j in range(nj):
        A[:,j]=args[j]
    writetable(name,A)


#CRP Added December 9th, 2010
# Technically, the ordering for Elegant file columns is arbitrary, as long as the SDDS header is accurate.  Here, we choose the order to match the xxs format.
def xxs_to_elegant(xxs,E_ref):
#  Is this wrong?
    gammaref=E_ref/E_ele_eV
    #elegantdist=xxs.copy()  #
    firstCol=xxs[:,0]
    NumParticles=size(firstCol)
    elegantdist=zeros((NumParticles,7))
    elegantdist[:,0]=xxs[:,0]
    elegantdist[:,1]=xxs[:,1]  #px
    elegantdist[:,2]=xxs[:,2]
    elegantdist[:,3]=xxs[:,3]   #py
    elegantdist[:,4]=-xxs[:,4]/c0
    elegantdist[:,5]=gammaref*(1 + xxs[:,5])   #pz, or gammaref-xxs[:,5]*gammaref
    for index in range(0,NumParticles,1):
            elegantdist[index,6]=index
    return elegantdist



#CRP 10/30/2012 Note: Depricated
# Due to the need to use sddsprintout to read the Elegant output, this 'Custom' format is essentially arbitrary
#    as long as the conversion back in to xxs is consistent.
def Custom1(name):
#  Line 1:  pCentral Charge s Particles  0   0   0
#  Line 2:  x  xp  y  yp  t  p  dt
    f=open(name,'r');    AS=readT(f);    f.close()
#    pCentral=AS[0,0]
#    Charge=AS[0,1]
#    s=AS[0,2]
#    NumParticles=AS[0,3]
    # extract invalid particles
#    ind_0=[]
#    for i in xrange(AS.shape[0]):
#        if AS[i,9]==0: ind_0.append(i)
#         ind_0.append(i)
#    AS=remove(AS,ind_0)
    n=size(AS[:,0])-1
    print('Custom1')
    print(AS)
#    sddsData=zeros((n,7))
#    sddsData[0,:]=AS[0,:]
#    sddsData[1:,:]=AS[2:,:]
    return AS


#CRP 10/30/2012 Note: Depricated
# CRP 11/10/11  Must match the ordering of whatever sddsprintout script/command is used. 
def Custom1_to_xxs(name):
    sddsData=Custom1(name)
    n=len(sddsData)-1
    print('Custom1_to_xxs')
    print(n)
#    print sddsData
#    sddsData=zeros((n,7))
#    n=size(sddsData[:,0])
    pCentral=sddsData[0,0]
    E_ref=pCentral*E_ele_eV
    print(pCentral)
#    Charge=sddsData[0,1]
    Charge=-3.2e-9
#    Charge[0:n-1]=Charge/n
    s=sddsData[0,2]
    NumParticles=sddsData[0,3]
#    AS[1:n,2]+=AS[0,2]; AS[1:n,5]+=AS[0,5]  #add reference particle parameters
    x=sddsData[1:,0].copy()
    xs=sddsData[1:,1].copy()
    y=sddsData[1:,2].copy()
    ys=sddsData[1:,3].copy()   
    t=sddsData[1:,4].copy()
    z=t*c0
    zav=mean(z)
    z=z-zav
    p=sddsData[1:,5].copy()
    delta=(p-pCentral)/pCentral
#    dt=AS[1:,6].copy()
    AS=sddsData[1:,:]
    AS[:,0]=x
    AS[:,1]=xs
    AS[:,2]=y
    AS[:,3]=ys
####### CRP Oct 5th, 2011:  Changed sign of z?
    AS[:,4]=-z
    AS[:,5]=delta
    q=delta-delta + Charge/n
    print(pCentral)
#    print sddsData
#    print name
    return [AS,q,zav,E_ref,Charge]



#CRP 10/30/2012 Note:  This is the "preferred" Elegant Format.  The old "custom 1" depended on an external script.  This function matches the format in the S2E_tools "sddsprintout" command.
# Due to the need to use sddsprintout to read the Elegant output, this 'SDDSout' format is essentially arbitrary
#    as long as the conversion back in to xxs is consistent.
def SDDSout(name):
#  Line 2: x  xp  y  yp  t  p  
#             m      m      s  gamma
    f=open(name,'r');    AS=readT(f);    f.close()
#    n=size(AS[:,0])-1
    n=size(AS[:,0])
    print('SDDSout')
    print(AS)
    return AS

#CRP 10/30/2012 Note: Converts the arbitrary Elegant format to xxs.
def SDDSout_to_xxs(name):
    sddsData=SDDSout(name)
    n=len(sddsData)
    print('SDDSout_to_xxs')
    print(n)

#    pCentral=sddsData[0,0]
    pCentral=mean(sddsData[:,5])
# Trying out a change:
#    E_ref=pCentral*E_ele_eV
    E_ref=(pCentral)*E_ele_eV
    print(pCentral, E_ref)

#    s=sddsData[0,2]
#    NumParticles=sddsData[0,3]
    NumParticles=n    

#    x=sddsData[1:,0].copy()
#    xs=sddsData[1:,1].copy()
#    y=sddsData[1:,2].copy()
#    ys=sddsData[1:,3].copy()   
#    t=sddsData[1:,4].copy()
#    z=t*c0
#    zav=mean(z)
#    z=z-zav
#    p=sddsData[1:,5].copy()
#    delta=(p-pCentral)/pCentral
#    print zav
    x=sddsData[:,0].copy()
    xs=sddsData[:,1].copy()
    y=sddsData[:,2].copy()
    ys=sddsData[:,3].copy()   
    t=sddsData[:,4].copy()
    z=t*c0
    zav=mean(z)
    z=z-zav
    p=sddsData[:,5].copy()
    delta=(p-pCentral)/pCentral
    print(zav)
#    x=sddsData[:,0].copy()
#    xs=sddsData[:,1].copy()
#    y=sddsData[:,2].copy()
#    ys=sddsData[:,3].copy()   
#    t=sddsData[:,4].copy()
#    z=t*c0
#    zav=mean(z)
#    z=z-zav
#    p=sddsData[:,5].copy()
#    delta=(p-pCentral)/pCentral

#    AS=sddsData[:,:]  # this was the regular
#    AS=sddsData[1:,:]
#    AS[:,0]=x
#    AS[:,1]=xs
#    AS[:,2]=y
#    AS[:,3]=ys
#    AS[:,4]=-z
#    AS[:,5]=delta

    AS=sddsData[:,:]
#    AS=zeros(6,NumParticles)
    AS[:,0]=x
    AS[:,1]=xs
    AS[:,2]=y
    AS[:,3]=ys
    AS[:,4]=-z
    AS[:,5]=delta

    q=zeros((n))
    q+=1.00; q_tot=sum(q) 

#    print pCentral

    return [AS,q,zav,E_ref]

#def(addNumParticles)


# CRP 11/10/11   Elegant requires an SDDS format.  Thankfully, the data is not divided in to pages, but this still requieres a very specific header structure.  First, a list of parameters, then a list of the units for the various columns.
# First few lines are the parameters, in order.
# Then the rest is a columned matrix of the main data, in the order listed.
#CRP 10/30/2012 Note:   Elegant is generally an arbitrary pain to write to, but this is enough for the simulation to start.
def writetableELEGANT(name,A,E_ref):
    gammaref=E_ref/E_ele_eV
    firstCol=A[:,0]
    NumParticles=size(firstCol)
    f=open(name,'w');
    head="SDDS1\n"
    f.write(head)
    head="&parameter name=Step, description=\"Simulation step\", type=long, &end \n"
    f.write(head)
    head="&parameter name=pCentral, symbol=\"p$bcen$n\", units=\"m$be$nc\", description=\"Reference beta*gamma\", type=double, &end\n"
    f.write(head)
    head="&parameter name=Charge, units=C, description=\"Beam charge\", type=double, &end\n"
    f.write(head)
    head="&parameter name=Particles, description=\"Number of particles\", type=long, &end\n"
    f.write(head)
    head="&column name=x, units=m, type=double, &end \n&column name=xp, type=double,   &end \n"
    f.write(head)
    head="&column name=y, units=m, type=double,  &end \n&column name=yp, type=double,  &end \n"
    f.write(head)
    head="&column name=t, units=s, type=double,  &end \n&column name=p, type=double,  &end \n&column name=particleID, type=long,  &end \n"
    f.write(head)
    head="&data mode=ascii, &end \n"
    f.write(head)
    head="1\n" 
    f.write(head)
    s=str(gammaref)
    f.write(s)
    head="\n nan \n"
    f.write(head)
    s=str(NumParticles)
    f.write(s)
    head="\n" 
    f.write(head)
    s=str(NumParticles)
    f.write(s)
    head="\n" 
    f.write(head)
    writeT(f,A);
    f.close()

#CRP 10/30/2012 Note: Finds the Beam Matrix?  I did not make this.
def makeBeamMatrix(xpx):
        beamMatrix=zeros((6,6))
        for n in range (0,6) :
                for m in range (0,6) :  #Range is definied strangely in Python, so this is really 6 by 6).
                        factor1=1;
                        factor2=1;
                        beamMatrix[n,m]=mean(xpx[:,n]*xpx[:,m])*factor1*factor2;
        return [beamMatrix]
    

# CRP 11/15/11: Created and tested.
def emittances(beamMatrix):
        print(beamMatrix)
        ex=numpy.sqrt(beamMatrix[0,0]*beamMatrix[1,1] - beamMatrix[0,1]*beamMatrix[0,1])
        ey=numpy.sqrt(beamMatrix[2,2]*beamMatrix[3,3] - beamMatrix[2,3]*beamMatrix[2,3])
        ez=numpy.sqrt(beamMatrix[4,4]*beamMatrix[5,5] - beamMatrix[4,5]*beamMatrix[4,5])
        return [ex,ey,ez]
    

# CRP 02/01/12:  Split fourier analysis in to this function, to be used by 1D, 2D, and 3D routines located in S2E_tools, as part of the Fourier analysis of Gaussian macroparticles.  "positions" is xxs[:,0], xxs[:,2], xxs[:,4], for x, y, and z analysis.  
def fourier_analysis(positions,q, lowWavelength, highWavelength,NumWavelengths, gaussianSize):
        waveStep=(highWavelength-lowWavelength)/(NumWavelengths-1)
        totalCharge=sum(q)
        wavelengthVector=zeros((NumWavelengths,1))
        frequencyVector=zeros((NumWavelengths,1))
        fourier=zeros((NumWavelengths,1))
        BFF=zeros((NumWavelengths,1))
        wavelengthList=arange(0,NumWavelengths+.01,1)
        cms=2.998e8
        sizeFactor=lowWavelength
        print('STD size (mm):', std(positions)*1000.0)
        print('GaussianSize / STD size:',  gaussianSize/std(positions))
        wavelengthList=sizeFactor**((wavelengthList)/NumWavelengths)
                

        frequencyList=2.998e8/wavelengthList

#        expectedMacroSize=1e-9
        expectedMacroSize= gaussianSize
        numParticles=len(positions)
# Guassian Breakdown
        minD=min(positions)
        maxD=max(positions)
                
        numBins=1000
        bins=arange(minD,maxD,(maxD-minD)/numBins)
        binSize=bins[2]-bins[1]
        I=0.0*bins
        print(binSize)
        histInd=0
        while histInd<numParticles-1:
                I= I + binSize*exp(-(bins-positions[histInd])**2.0/(2.0*gaussianSize**2.0))/sqrt(2.0*pi*gaussianSize**2.0);
                histInd=histInd+1  
#wavelength=c/f
#wavelength=(cms/frequency)
        waveInd=0
        while waveInd < NumWavelengths:
                wavelength=wavelengthList[waveInd]
                frequency=cms/wavelength
                wavelengthVector[waveInd]=wavelength
                frequencyVector[waveInd]=frequency
#                fourier[waveInd]=(1.0/totalCharge)*exp(-expectedMacroSize**2.0*4.0*pi*pi/wavelength**2.0)*(sum(I[:]*sin(2.0*pi*bins[:]/wavelength)) + sum(I[:]*cos(2.0*pi*bins[:]/wavelength) )  )
                fourier[waveInd]=(1.0/totalCharge)*exp(-expectedMacroSize**2.0*4.0*pi*pi/(cms/frequency)**2.0)*(sum(I[:]*sin(2.0*pi*bins[:]/(cms/frequency))) + sum(I[:]*cos(2.0*pi*bins[:]/(cms/frequency)) )  )
                BFF[waveInd]=(1.0/totalCharge**2.0)*exp(-expectedMacroSize**2.0*4.0*pi*pi/wavelength**2.0)*(sum(I[:]**2.0*sin(2.0*pi*bins[:]/wavelength))**2 + sum(I[:]**2*cos(2.0*pi*bins[:]/wavelength) )**2.0  )
#                        print fourier[waveInd], BFF[waveInd]
                waveInd=waveInd+1

        print('fourier done')
        return [wavelengthVector,frequencyVector,fourier,BFF, bins, I]
        
        #conventionalEmittances(ind,1)= sqrt(x2*px2 -xpx^2);
#conventionalEmittances(ind,2) = sqrt(y2*py2 -ypy^2);
#conventionalEmittances(ind,3) = sqrt(z2*pz2 -zpz^2);





        








    
    
    






    
    
    
        
    
    



    
