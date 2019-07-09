import os
import glob
import shutil
import random

## CRP: 11/10/11 Numerous capabilities added:  CSRtrack_Parallel, LPS_straigthen, LPSplot, plot_show
## CRP: 10/05/11 Flipped CSR fmt1 "write" sign for z.
## CRP: 08/18/11 Added shutil for copying (to prevent permission errors)

import sys
import os
from pylab import *

from particle_tools import*
from wakes import *

from numpy import *
import numpy as np
import random
import matplotlib.pyplot as plt



# CRP note:   #OG denotes a function already implemented when PP and I received GlueTrack.

#OG
def RunTask(actions,fig_dir,html,indfig=0,start=0,stop=100000):
    stop=min([stop+1,len(actions)]);    start=max([0,start]);
    print("pp in Run Task")
    for i in range(start,stop):
        action=actions[i];
        name=action['name'];
        tasks=action['tasks'];
        comment=action['comment']
        print(action)
        if name=='processor':
            print('\n');print(comment);        print(name,':',i); 
            html.write('<H1>%s</H1>\n' %comment)
            indfig=Processor(tasks,fig_dir,html,indfig)
            html.write('<H1>done</H1>\n')
            continue
        if name=='astra':
            print('\n');print(comment);        print(name,':',i);     
            html.write('<H1>%s</H1>\n' %comment)
            indfig=ASTRA(tasks,fig_dir,html,indfig)
            html.write('<H1>done</H1>\n')
            continue
## PP 02/24/2011
        if name=='impactT':
            print('\n');print(comment);        print(name,':',i);     
            html.write('<H1>%s</H1>\n' %comment)
            indfig=IMPACTT(tasks,fig_dir,html,indfig)
            html.write('<H1>done</H1>\n')
            continue
## PP 02/25/2011
        if name=='impactZ':
            print('\n');print(comment);        print(name,':',i);     
            html.write('<H1>%s</H1>\n' %comment)
            indfig=IMPACTZ(tasks,fig_dir,html,indfig)
            html.write('<H1>done</H1>\n')
            continue
        if name=='CSRTrack':
            print('\n');print(comment);        print(name,':',i);     
            html.write('<H1>%s</H1>\n' %comment)
            indfig=CSRTrack(tasks,fig_dir,html,indfig)
            html.write('<H1>done</H1>\n')
            continue
## CRP 11/09/2011
# Runs the parallel version of CSRtrack, but requires several extra Properties (see function CSRTrack_Parallel)
        if name=='CSRTrack_Parallel':
            print('\n');print(comment);        print(name,':',i);     
            html.write('<H1>%s</H1>\n' %comment)
            indfig=CSRTrack_Parallel(tasks,fig_dir,html,indfig)
            html.write('<H1>done</H1>\n')
            continue
## CRP 8/18/11  Added Elegant
        if name=='elegant':
            print('\n');print(comment);        print(name,':',i);     
            html.write('<H1>%s</H1>\n' %comment)
            indfig=ELEGANT(tasks,fig_dir,html,indfig)
            html.write('<H1>done</H1>\n')
            continue
# CRP 08/15/11  This was an attempt, but is not the method used.
#        if name=='SDDSprocess':
#            print '\n';print comment;        print name,':',i;     
#            html.write('<H1>%s</H1>\n' %comment)
#            indfig=SDDSprocess(tasks,fig_dir,html,indfig)
#            html.write('<H1>done</H1>\n')
#            continue
        if name=='analysis':
            print('\n'); print(comment);        print(name,':',i);     
            html.write('<H1>%s</H1>\n' %comment) 
#            print "pp go to analysis"
            indfig=Analysis(tasks,fig_dir,html,indfig)
            html.write('<H1>done</H1>\n')
    return indfig


def Processor(tasks,fig_dir,html,indfig):
    for i in range(len(tasks)):
        task=tasks[i];  name=task['name']; properties=task['properties']
        if name=='read':
            print('task: ', name)
#### CRP July 18th, 2012:  I dislike adding the ./, as it SHOULD be able to read files from anywhere.  
###            in_file='./'+properties['in_file']
            in_file=properties['in_file']
            print(in_file)
            format=properties['format']
            extr_part_1=sort_len=False; E_ref=0;
            if 'extr_part_1' in properties: extr_part_1=eval(properties['extr_part_1'])
            if 'sort_len' in properties: sort_len=eval(properties['sort_len'])
            if 'E_ref' in properties:  E_ref=eval(properties['E_ref'])
# CRP modified 10/27/2011 to pass the properties information deeper into GlueTrack
#            Dh=particle_reader(in_file,format,E_ref,extr_part_1,sort_len)
            Dh=particle_reader(in_file,format,E_ref,extr_part_1,sort_len,properties)
            N=Dh[0]; xpx=Dh[1];  xxs=Dh[2]; q=Dh[3]
            E_ref=Dh[4]; q_tot=Dh[5];  zav=Dh[6]; ct=Dh[7]
#PP            plt.plot(Dh[1][:,4], Dh[1][:,5],'.')
#PP            plt.show()
            print('- N=',N,'  E_ref=',E_ref)
            print('q_tot=',q_tot,'  zav=',zav,'  ct=',ct)
            if 'in_file2' in properties:
                in_file2=properties['in_file2']
                Dh2=particle_reader(in_file2,format,E_ref,extr_part_1,sort_len)
                xpx2=Dh2[1];  xxs2=Dh2[2];
            continue
        
        if name=='center':
            print('task: ', name)
            if 'interv' in properties:
                za,zb=eval(properties['interv'])
                w=where(logical_and(xxs[:,4]<zb,xxs[:,4]>za),1,0)
            else:
                w=q
            centerz=centerx=centerxs=centery=centerys=False
            if 'centerz' in properties: centerz=eval(properties['centerz'])
            if 'centerx' in properties: centerx=eval(properties['centerx'])
            if 'centerxs' in properties: centerxs=eval(properties['centerxs'])
            if 'centery' in properties: centery=eval(properties['centerx'])
            if 'centerys' in properties: centerys=eval(properties['centerys'])
            if centerz:
                z_av=MoM1(xxs[:,4],w,0)
                print('z_av(m)=',z_av)
                xxs[:,4]-=z_av  
            if centerx:
                x_av=MoM1(xxs[:,0],w,0)
                print('x_av(m)=',x_av)
                xxs[:,0]-=x_av  
            if centerxs:
                xs_av=MoM1(xxs[:,1],w,0)
                print('xs_av(m)=',xs_av)
                xxs[:,1]-=xs_av 
            if centery:
                y_av=MoM1(xxs[:,2],w,0)
                print('y_av(m)=',y_av)
                xxs[:,2]-=y_av  
            if centerys:
                ys_av=MoM1(xxs[:,3],w,0)
                print('ys_av(m)=',ys_av)
                xxs[:,3]-=ys_av
            continue

# This function shifts the particle energy using simple scaling, and adjust the parameter E.  This is done by 
# Does not work yet???
        if name=='recenter_Energy':
                print('old average energy:', mean((xxs[:,5]+1)*E_ref))
                energies=(xxs[:,5]+1)*E_ref
                shiftedDelta=energies/E_ref-1
                print('old average delta:', mean(shiftedDelta))
                print('old sigma delta:', std(xxs[:,5]))                
                new_energy=eval(properties['new_energyMeV'])*1000000
                E_ref=new_energy
                new_energies=E_ref*(1+shiftedDelta)
                xxs[:,5]=new_energies/E_ref-1
                print('new average energy:', mean((xxs[:,5]+1)*E_ref))
                print('new average delta:', mean(xxs[:,5]))
                print('new sigma delta:', std(xxs[:,5]))
                continue

        if name=='extract_slice_centr':
            print('task: ', name)
            M=eval(properties['M'])
            sigmaI=eval(properties['sigmaI'])
            za,zb=eval(properties['interv'])
            kso=sort1(xxs[:,4]); xxs=sort2(xxs,kso); q=sort2(q,kso)
            iso=sort1(kso)
            qsu=Sum1(q*1e-300+1,q);
            B=s_to_cur(xxs[:,4],sigmaI,1e-9,3e8)
            xsu=Sum1(xxs[:,0],q);   ysu=Sum1(xxs[:,2],q);
            xssu=Sum1(xxs[:,1],q);   yssu=Sum1(xxs[:,3],q);
            def av(Xsu,N,M):
                out=zeros(N,Float)
                for i in range(N):
                    m=min([0.5*M,0.5*i,0.5*(N-1-i)])
                    m=int(floor(m+0.5))
                    if m>0: out[i]=(Xsu[i-m]-Xsu[i+m])/(qsu[i-m]-qsu[i+m])
                return out
            xav=av(xsu,N,M);    yav=av(ysu,N,M)
            xsav=av(xssu,N,M);    ysav=av(yssu,N,M)
            plot(xxs[:,4],xav,'r',xxs[:,4],yav,'b');xlim(za,zb);
            fig_file='fig%g.png' %indfig; indfig=indfig+1
            savefig(fig_dir+fig_file);   hold(False);html.write('<IMG SRC="%s">\n' %fig_file)
            plot(B[:,0],B[:,1]);xlim(za,zb);
            fig_file='fig%g.png' %indfig; indfig=indfig+1
            savefig(fig_dir+fig_file);   hold(False);html.write('<IMG SRC="%s">\n' %fig_file)
            plot(xxs[:,4],xsav,'r',xxs[:,4],ysav,'b');xlim(za,zb);
            fig_file='fig%g.png' %indfig; indfig=indfig+1
            savefig(fig_dir+fig_file);   hold(False);html.write('<IMG SRC="%s">\n' %fig_file)
            xxs[:,0]=xxs[:,0]-xav;  xxs[:,1]=xxs[:,1]-xsav
            xxs=sort2(xxs,iso); q=sort2(q,iso)
            continue

        if name=='define_ref_part':
            print('task: ', name)
            xxsr=zeros(6,Float)
            for k in range(6): xxsr[k]=MoM1(xxs[:,k],q,0)
            print(xxsr)
            xxs=concatenate((xxsr[NewAxis,:],xxs))
            OO_0=zeros(1,Float); q=concatenate((OO_0,q))
            continue


### CRP 11/09/2011
### LPS_straighten        
### Removes LPS curvature, simulating a third-harmonic cavity.
        if name=='LPS_straighten':
#####       numpy.polyfit(x, y, deg, rcond=None, full=False, w=None, cov=False)
            xpx=xxs_to_xpx(xxs,E_ref)
            Coeff=np.polyfit(xpx[:,4],xpx[:,5],2)
            print(Coeff)
#            ZZS=xxs[:,4:6];
#            za,zb=eval(properties['interv'])
#            w=where(logical_and(xxs[:,4]<zb,xxs[:,4]>za),1,0)
#            alpha_z,beta_z,gamma_z,epsilon_z=ShowTwissPar(ZZS,w);
###            xxs[:,5]=xxs[:,5]+Coeff[0]-Coeff[1]*xxs[:,4] -Coeff[2]*xxs[:,4]*xxs[:,4]
            xpx[:,5]=xpx[:,5] - Coeff[0]*xpx[:,4]*xpx[:,4]
#            xpx[:,5]=Coeff[0]+Coeff[1]*xpx[:,4]
            xxs=xpx_to_xxs(xpx,E_ref)
            continue



        if name=='drift':
            print('task: ', name)
            length=eval(properties['length'])
            xxs[:,0]+=length*xxs[:,1]
            xxs[:,2]+=length*xxs[:,3]
            continue

        if name=='drift_to':
            print('task: ', name)
            z0=eval(properties['z0'])
            z_ist=zav; length=z0-z_ist
            xxs[:,0]+=length*xxs[:,1];  xpx[:,0]=xxs[:,0]
            xxs[:,2]+=length*xxs[:,3];  xpx[:,2]=xxs[:,2]
            continue
            

### CRP 11/14/11:  entering 0 for all of the parameters [a,b,m] of xpar OR ypar will not change that dimension at all.
        if name=='matching':
            print('task: ', name)
#            a_x,b_x,m_x=eval(properties['xpar']);
#            a_y,b_y,m_y=eval(properties['ypar']);
            m_x=0; m_y=0
            a_x,b_x,e_x=eval(properties['xpar']);
            a_y,b_y,e_y=eval(properties['ypar']);
            
            za,zb=eval(properties['interv'])
            XXS=xxs[:,0:2];  YYS=xxs[:,2:4];  ZZS=xxs[:,4:6];
            w=where(logical_and(xxs[:,4]<zb,xxs[:,4]>za),1,0)
            print('x, y, z Twiss parameters')
            alpha_x,beta_x,gamma_x,epsilon_x=ShowTwissPar(XXS,w);
            alpha_y,beta_y,gamma_y,epsilon_y=ShowTwissPar(YYS,w);
            alpha_z,beta_z,gamma_z,epsilon_z=ShowTwissPar(ZZS,w);
            if e_x==0:
                e_x=epsilon_x
            if e_y==0:
                e_y=epsilon_y
            if b_x==0:
                e_x=epsilon_x
                a_x=alpha_x
                b_x=beta_x
            if b_y==0:
                e_y=epsilon_y
                a_y=alpha_y
                b_y=beta_y
            Tx=T(E_ref,(alpha_x,beta_x,0),E_ref,(a_x,b_x,m_x));print('Tx=',Tx)
            Ty=T(E_ref,(alpha_y,beta_y,0),E_ref,(a_y,b_y,m_y));print('Ty=',Ty)
#PP            xxs2=matrixmultiply(XXS,transpose(Tx))
#PP            yys2=matrixmultiply(YYS,transpose(Ty))
#PP 11/14/2010 updated to match new numpy1.5 syntax
            xxs2=dot(XXS,transpose(Tx))*sqrt(e_x/epsilon_x)
            yys2=dot(YYS,transpose(Ty))*sqrt(e_y/epsilon_y)
            ShowTwissPar(xxs2,w);ShowTwissPar(yys2,w)
            xxs[:,0]=xxs2[:,0];xxs[:,1]=xxs2[:,1];
            xxs[:,2]=yys2[:,0];xxs[:,3]=yys2[:,1];
            continue

#            ZZS=xxs[:,4:6];
#            za,zb=eval(properties['interv'])
#            w=where(logical_and(xxs[:,4]<zb,xxs[:,4]>za),1,0)
#            alpha_z,beta_z,gamma_z,epsilon_z=ShowTwissPar(ZZS,w);


# CRP 08/15/11
# This function reads a specific Elegant "Cross" file into a file format I call "Custom1".  How data is dumped with any sddsprintout script is essentially arbitrary, and easily change.  Only the first 4 parts of the first line are really of interest: the other three are dumped only to keep the first line 7 values, like the rest of the matrix.
# Line 1: pCentral(gamma) BunchCharge(C)  s(m)  NumParticles    Pass PassLength PassCentralTime
# Line 2+:  x(m)   xp(unitless)   y(m)   yp(unitless)   t(s)   p(gamma)   dt(s)
        if name=='SDDSprocess':
            command=properties['command']
            work_dir='./'+properties['work_dir']
            command=command.replace(',',' ')
            os.chdir(work_dir)
            src='./MiscScripts'
            print(command)
            fullCommand=properties['deck_inp'] + ' ' + properties['dist_inp'] + ' > TempDistributions/temp'
            os.system(fullCommand)
            fin = open( 'TempDistributions/temp', "r" )
            data_list = fin.readlines()
            fin.close()
            print(fullCommand)
            del data_list[1]
            fout = open("TempDistributions/processedSDDSfile", "w")
            fout.writelines(data_list)
            fout.close()

            continue
            
#        if name=='Lmatching':
#            print 'task: ', name
#            m_z=0;
#            chirp,sig_z,e_z=eval(properties['zpar']);
#            b_z=sig_z*sig_z/e_z
#            a_z=-chirp*b_z
#            print "target longitudinal C-S:", a_z, b_z, e_z
#            za,zb=eval(properties['interv'])
#            ZZS=xxs[:,4:6]; 
#            w=where(logical_and(xxs[:,4]<zb,xxs[:,4]>za),1,0)
#            alpha_z,beta_z,gamma_z,epsilon_z=ShowTwissPar(ZZS,w);
#            print "CHIRP [m^-1]=", -alpha_z/beta_z ;
#            if e_z==0:
#                e_z=epsilon_z
#            if b_z==0:
#                e_z=epsilon_z
#                a_z=alpha_z
#                b_z=beta_z
#           Tz=T(E_ref,(alpha_z,beta_z,0),E_ref,(a_z,b_z,m_z));print 'Tz=',Tz
##PP 11/14/2010 updated to match new numpy1.5 syntax
#            zzs2=dot(ZZS,transpose(Tz))*sqrt(e_z/epsilon_z)
#            ShowTwissPar(zzs2,w);
#            xxs[:,4]=zzs2[:,0];xxs[:,5]=zzs2[:,1];
#            continue


# Reorganzied by CRP, 11/07/2011, to allow for e_z=0 without throwing an error.  Previously, epsilon_z was defined after the e_z==0 check.
        if name=='Lmatching':
            print('task: ', name)
            za,zb=eval(properties['interv'])
            ZZS=xxs[:,4:6]; 
            w=where(logical_and(xxs[:,4]<zb,xxs[:,4]>za),1,0)
            alpha_z,beta_z,gamma_z,epsilon_z=ShowTwissPar(ZZS,w);
            m_z=0;
            chirp,sig_z,e_z=eval(properties['zpar']);
            if e_z==0:
                e_z=epsilon_z

# CRP 12/29/11:  Added "default" values.
            if sig_z==0:
                sig_z=sqrt(beta_z*epsilon_z)
#            if chirp==0:
#                chirp=-alpha_z/beta_z

            b_z=sig_z*sig_z/e_z
            a_z=-chirp*b_z
            print("target longitudinal C-S:", a_z, b_z, e_z)          
            print("CHIRP [m^-1]=", -alpha_z/beta_z) ;
            if b_z==0:
                e_z=epsilon_z
                a_z=alpha_z
                b_z=beta_z
            Tz=T(E_ref,(alpha_z,beta_z,0),E_ref,(a_z,b_z,m_z));print('Tz=',Tz)
##PP 11/14/2010 updated to match new numpy1.5 syntax
            zzs2=dot(ZZS,transpose(Tz))*sqrt(e_z/epsilon_z)
            ShowTwissPar(zzs2,w);
            xxs[:,4]=zzs2[:,0];xxs[:,5]=zzs2[:,1];
            continue




#CRP 10/24/11: Function added and tested.  Would like to eventually modify the TargetSize to be specified by the user, but for now, this is good.
        if name=='downSample':
#             Down sample a particle distribution, by pruning unsorted lines from xpx and the corresponding lines from q.  At each step, q needs to be increased.
#def downSampleParticles(xpx, q, N_p,targetSize):
            N_p=size(q)
            print(size(q), size(xxs))
            totalQ=sum(q)
#            targetSize=10000
            targetSize=eval(properties['target_size'])
            adjusted_q=totalQ/targetSize
            particle=0
            validParticles=list(range(2,N_p,1))
            numToCut= N_p-targetSize
            rowsToDelete=random.sample(validParticles,numToCut)
            print(len(rowsToDelete))
            xxs=numpy.delete(xxs, rowsToDelete, 0)
            xpx=numpy.delete(xpx, rowsToDelete, 0)
           
            q=numpy.delete(q, rowsToDelete, 0)

            q[:]=adjusted_q
        
            print(size(q), len(q), size(xxs), size (xpx))

            continue




#CRP 11/14/11: Plots various beam aspects, and SHOWS it, halting the GlueTrack operation until the windows are closed.
        if name=='xpxPlot':
            xpx=xxs_to_xpx(xxs,E_ref)
            #matplotlib
            plt.figure()
            plt.scatter(xpx[1:,4],xpx[1:,5],)
            plt.figure()
            plt.scatter(xpx[1:,0],xpx[1:,2],)
            plt.figure()
            plt.scatter(xpx[1:,0],xpx[1:,2],)
            print(sum(q))
            plt.show()

            continue

        if name=='spotPlot':
#            xpx=xxs_to_xpx(xxs,E_ref)
            #matplotlib
            plt.figure()
            plt.scatter(xpx[1:,0],xpx[1:,2],)
#            plt.figure()
#            plt.scatter(xpx[1:,0],xpx[1:,2],)
#            plt.figure()
#            plt.scatter(xpx[1:,0],xpx[1:,2],)
            print(sum(q))
            plt.show()

            continue

#CRP 11/14/11:  Creates the figure for the an LPS plot, with the appropriate title, but must be coupled with the plot_show command.  The suggested usage is to insert LPSplot throughout the GlueTrack file, and then to use plot_show at the end to display the LPS evolution.
#   I will add more Properties at a later time.
        if name=='LPSplot':
            xpx=xxs_to_xpx(xxs,E_ref)
            #matplotlib
            plt.figure()
            plt.scatter(xpx[1:,4],xpx[1:,5],)
            LPStitle=properties['title']
            plt.title(LPStitle)
            print(sum(q))
##            plt.show()

            continue



#CRP 11/14/11: The show command, to plot any previous figure() created my other GlueTrack Processes (as of this writing, LPS plot, but I will add more later.
        if name=='plot_show':
            plt.show()
            continue

#CRP 08/29/12:  Perform the thin-lens TDC transformation for a given strength T.  This is only a transformation on delta and xprime.  Should I split this off in to a particle_tools functions?  I don't understand why this isn't working.
        if name=='thin_TDC_x':
            kickstrength=eval(properties['Kick_Strength_T']);
            print('task: Thin TDC with kick strength of ' + str(kickstrength))
            xxs[:,1]=xxs[:,1] + kickstrength*(xxs[:,4])
            xpx=xxs_to_xpx(xxs,E_ref)
            print(xpx[:,5])
            print(xxs[:,5])
            print(std(xxs[:,5]), mean(xxs[:,5]))
#            xxs[:,5]=(xxs[:,5]/mean(xxs[:,5]) - 1)  + kickstrength*xxs[:,1]
            xxs[:,5]=xxs[:,5] + kickstrength*(xxs[:,0])
            print(std(xxs[:,5]), mean(xxs[:,5]))
            continue

# "One Cell" TDC, where the one cell means that we use the maximum value for the R65 term.  More cells will reduce it.
        if name=='thick_TDC_oneCell_x':
            kickstrength=eval(properties['Kick_Strength_T']);
            cavitylength=eval(properties['Cavity_Length_m']);
            print('task: Thick TDC with kick strength of ' + str(kickstrength) + ' and Length of ' + str(cavitylength))
            xxs[:,0]=xxs[:,0] +  cavitylength*xxs[:,1] + cavitylength*kickstrength*(xxs[:,4])/2
            xxs[:,1]=xxs[:,1] + kickstrength*(xxs[:,4])

# Don't forget the drift in y!
            xxs[:,2]=xxs[:,2] + cavitylength*(xxs[:,3])
#            xpx=xxs_to_xpx(xxs,E_ref)
#            print xpx[:,5]
#            print xxs[:,5]
#            print std(xxs[:,5]), mean(xxs[:,5])
#            xxs[:,5]=(xxs[:,5]/mean(xxs[:,5]) - 1)  + kickstrength*xxs[:,1]
            xxs[:,5]=xxs[:,5] + kickstrength*(xxs[:,0]) + cavitylength*kickstrength*(xxs[:,1])/2  + kickstrength*kickstrength*cavitylength*(xxs[:,4])/4
            print(std(xxs[:,5]), mean(xxs[:,5]))
            continue

            

#PP 12/15/2010
# slice the beam by introducing a multislit mask
# parameters are
#  direction = 0, 1, 2, 3, 4, 5   for x x' y y' z z'
#  width = a number         the slit(s) width
#  spacing= a number        the slits center-to-center separation
#  position = a number        the center slit position
        if name=='slits':
            print('task: ', name)
            direction=eval(properties['direction']);
            width=eval(properties['width']);
            offset=eval(properties['offset']);
            num=eval(properties['number']);
            spacing=eval(properties['spacing']);
            nn=len(xxs)
            ind=0
            posmin=-num/2.0*spacing+offset;
            posmax= num/2.0*spacing+offset;
            Xslit=xxs;
#PP this double-loop needs to be improved... 
            for j in range(num):
               if num==1:
                  slitpos=offset
               if num>1:
                  posstep = (posmax-posmin)/(num-1)
                  slitpos=posmin+j*posstep;
               print('slitpos:', slitpos)
               print('j:', j)
               for i in range(nn):
                  if abs(Xslit[i,direction]-slitpos)<width/2.0:
#                    print ind, j
                    for k in range(6):
                      xxs[ind,k]=Xslit[i,k]
                    ind=ind+1
# last populated index is      ind=ind-1
            print('Transmission:', ind)
#            xxs[:,0]=Xslit[:,0];xxs[:,1]=Xslit[:,1];
#            xxs[:,2]=Xslit[:,2];xxs[:,3]=Xslit[:,3];
#            xxs[:,4]=Xslit[:,4];xxs[:,5]=Xslit[:,5];
            xxs=delete(xxs, s_[ind:nn+1], axis=0)   
            N=ind-1
            print(len(xxs), nn, N)
            continue

#Slit is not robust enough for a complicated mask.  "masking_statement" allows for an arbitrary logic statement.  This will be difficult to implement and test.  Started September 6th 2012.
        if name=='masking_statement_xxs':
#            xxs=xpx_to_xxs(xpx,E_ref)
            #Logic is a line-by-line reader, an if statement that dumps in to Transmitted, an else that dumps in to Masked, then a copy of them in to xxs and xxp, cutting unused lines.  Important: have copies of both px and xp?  What steps need to be taken to make sure that the length is okay?  Look at random downsample.  Use simple x<0 and y<0 to test (spot size should be very distinct.).  I should use a temp, then outright delete xxs and xpx, then copy over.
#            xstd=std(xxs[0,:])
#            xpstd=std(xxs[1,:])
#            ystd=std(xxs[2,:])
#            ypstd=std(xxs[3,:])
#            zstd=std(xxs[4,:])
#            deltastd=std(xxs[5,:])

# CRP indices were flipped!  Fixed September 13th, 2012
            xstd=std(xxs[:,0])
            xpstd=std(xxs[:,1])
            ystd=std(xxs[:,2])
            ypstd=std(xxs[:,3])
            zstd=std(xxs[:,4])
            deltastd=std(xxs[:,5])
## These values are calculated to allow for quick scaling based on the size of the beam.
            print(xstd, xpstd, ystd, ypstd, zstd, deltastd)

            logicString=properties['logic_statement_xxs']
#            print eval(logicString)
            print(logicString)
#            if ( (xxs[line,0] <.1*xstd) and (xxs[line,0] > -.1*xstd) and (xxs[line,2] < ymax) and  (xxs[line,2] > ymin)  ) or  (xxs[line,0] <1.1*xstd) and (pxxs[line,0] > 0.9*xstd) or (xxs[line,0] > -1.1*xstd) and (xxs[line,0] < -0.9*xstd):
            numTrans=0
            numMasked=0
            line=0

            tempTransmittedParticles=xxs[0,:]
            tempMaskedParticles=xxs[0,:]
            tempq=q[0]


            while line < len(xxs):            
#                    if (xxs[line,0] > 0 ) or (xxs[line,2] > 0):
                    if eval(logicString):
#                exec logicString
                        tempTransmittedParticles=vstack((tempTransmittedParticles, xxs[line,:]))
                        tempq=vstack((tempq, q[line]))
                        numTrans+=1
                    else:
                        tempMaskedParticles=vstack((tempMaskedParticles, xxs[line,:]))
                        numMasked+=1
                    line=line+1

            print(len(xxs), numTrans, numMasked)
                    
            del xxs
            del xpx
            del q


            TransmittedParticles=tempTransmittedParticles[1:,:]
#####            MaskedParticles=tempMaskedParticles[1:,:]  # Commented out, as 0 masking would crash.
            q=tempq[1:]
            
            N_p=len(TransmittedParticles)
            q=tempq
            print(N_p, sum(q), size(TransmittedParticles))


            xxs=TransmittedParticles
            xpx=xxs_to_xpx(xxs,E_ref)



            continue
# I need to figure out what the way to deal with q is.
#while line < len(preMasked):
#        if ( (preMasked['x'][line] <.1*xstd) and (preMasked['x'][line] > -.1*xstd) and (preMasked['y'][line] < ymax) and  (preMasked['y'][line] > ymin)  ) or  (preMasked['x'][line] <1.1*xstd) and (preMasked['x'][line] > 0.9*xstd) or (preMasked['x'][line] > -1.1*xstd) and (preMasked['x'][line] < -0.9*xstd):
#                tempTransmittedParticles=vstack((tempTransmittedParticles, preMasked[:][line]))
#        else:
#                tempMaskedParticles=vstack((tempMaskedParticles, preMasked[:][line]))
#        line=line+1



        if name=='add_wake':
            print('task: ', name)
            #Current distribution
            sigma=eval(properties['sigmaI']);
            za,zb=eval(properties['m_interv'])
            wake=eval(properties['wake'])
            wake_delta=eval(properties['wake_delta'])
            B1=s_to_cur(xxs[:,4],sigma,q_tot,c0)
            C1=B1.copy(); C1[:,1]=C1[:,1]/c0;
            C=wake_conv(C1,wake,wake_delta);
            koef=max(abs(C[:,1]))/max(abs(B1[:,1]));
            plot(-C[:,0],C[:,1],-B1[:,0],B1[:,1]*koef);xlim(za,zb)
            xlabel('z[m]');ylabel('wake]');title('wake')
            fig_file='fig%g.png' %indfig; indfig=indfig+1
            savefig(fig_dir+fig_file);   hold(False);html.write('<IMG SRC="%s">\n' %fig_file)
            int1d=interpolate.interp1d(C[:,0],C[:,1],axis=0)
            y=int1d(xxs[:,4])/E_ref
            plot(xxs[:,4],xxs[:,5],'.r',xxs[:,4],xxs[:,5]+y,'.b');
            xlabel('z[m]');ylabel('E[eV]');title('longitudinal phase space with wake')
            fig_file='fig%g.png' %indfig; indfig=indfig+1
            savefig(fig_dir+fig_file);   hold(False);html.write('<IMG SRC="%s">\n' %fig_file)
            xxs[:,5]=xxs[:,5]+y;
            continue

#PP 01/08/2016 copy beam to make an array 
#
        if name=='multibeam_cart':
            
            print('task: ', name)
            print('shape:',np.shape(xxs))
            pitchx=eval(properties['pitchx']);
            pitchy=eval(properties['pitchy']);
            pitchx_p=eval(properties['pitchx_p']);
            pitchy_p=eval(properties['pitchy_p']);
            numx  =eval(properties['numx']);
            numy  =eval(properties['numy']);
            resamp = 0
            
            if 'resample' in properties:
                 resamp =eval(properties['resample'])

            print('resamp=',resamp)
            print('pitch=',pitchx, pitchy, pitchx_p, pitchy_p)
            posxmin=-(numx-1.)/2.0*pitchx;
            posxmax= (numx-1.)/2.0*pitchx;
            posymin=-(numy-1.)/2.0*pitchy;
            posymax= (numy-1.)/2.0*pitchy;

            posx_pmin=-(numx-1.)/2.0*pitchx_p;
            posx_pmax= (numx-1.)/2.0*pitchx_p;
            posy_pmin=-(numy-1.)/2.0*pitchy_p;
            posy_pmax= (numy-1.)/2.0*pitchy_p;
            
            print('pitch_pos=',posxmin, posxmax, posymin, posymax)
            
#PP this double-loop needs to be improved... 
            for j in range(numy):
               for i in range(numx):
                  if resamp>0:
                     index=random.sample(list(range(0, len(q))), resamp)
                     Xtemp=xxs[index,:]
                     Qtemp=q[index]
                  if resamp==0:   
                     Xtemp=xxs.copy()
                     Qtemp=q.copy()
                     
                  Xpos=posxmin+i*pitchx
                  Ypos=posymin+j*pitchy
                  X_ppos=posx_pmin+i*pitchx_p
                  Y_ppos=posy_pmin+j*pitchy_p
#                  print Xpos, Ypos
                  
                  if pitchx>0:
                     Xtemp[:,0]=Xpos+Xtemp[:,0]
                  if pitchy>0:
                     Xtemp[:,2]=Ypos+Xtemp[:,2]
                  if pitchx_p>0:
                     Xtemp[:,1]=X_ppos+Xtemp[:,1]
#                     Xtemp[:,5]=Xtemp[:,5]-X_ppos
                  if pitchy_p>0:
                     Xtemp[:,3]=Y_ppos+Xtemp[:,3]
#                     Xtemp[:,5]=Xtemp[:,5]-Y_ppos

                  
                  if (i==0) & (j==0):
                     Xnew=Xtemp.copy()
                     Qnew=Qtemp.copy()
#                     print i, j
                  else: 
                     Xnew=np.vstack((Xnew, Xtemp))
                     Qnew=np.hstack((Qnew, Qtemp))
                     
            xxs = Xnew
            q   = Qnew
            print(q)
            N=np.shape(xxs)[0]     
            print(np.shape(xxs))
            
            plt.subplot (1,2,1)
            plt.plot (xxs[:,0], xxs[:,2],'.')
            plt.subplot (1,2,2)
            plt.plot (xxs[:,1], xxs[:,3],'.')
#            plt.plot (xxs[:,4], xxs[:,5],'.')
            plt.show()
            continue


        if name=='multibeam_hex':
            
            print('task: ', name)
            print('shape:',np.shape(xxs))
            pitchx=eval(properties['pitchx']);
            #pitchy=eval(properties['pitchy']);
            pitchy= pitchx/1.7320508076/2.0 ;
            pitchx_p=eval(properties['pitchx_p']);
            pitchy_p=pitchx_p/1.7320508076/2.0;
            numx  =eval(properties['numx']);
            numy  =eval(properties['numy']);
            resamp = 0
            
            if 'resample' in properties:
                 resamp =eval(properties['resample'])

            print('resamp=',resamp)
            posxmin=-(numx-1.)/2.0*pitchx;
            posxmax= (numx-1.)/2.0*pitchx;
            posymin=-(numy-1.)/2.0*pitchy;
            posymax= (numy-1.)/2.0*pitchy;

            posx_pmin=-(numx-1.)/2.0*pitchx_p;
            posx_pmax= (numx-1.)/2.0*pitchx_p;
            posy_pmin=-(numy-1.)/2.0*pitchy_p;
            posy_pmax= (numy-1.)/2.0*pitchy_p;
            
#PP this double-loop needs to be improved... 
            for j in range(numy):
               for i in range(numx):
                  if resamp>0:
                     index=random.sample(list(range(0, len(q))), resamp)
                     Xtemp=xxs[index,:]
                     Qtemp=q[index]
                  if resamp==0:   
                     Xtemp=xxs.copy()
                     Qtemp=q.copy()
                  Xpos=posxmin+i*pitchx
                  Ypos=posymin+j*pitchy
                  X_ppos=posx_pmin+i*pitchx_p
                  Y_ppos=posy_pmin+j*pitchy_p
                  if j%2==1 :
                      Xpos=posxmin+(i-0.5)*pitchx
                      X_ppos=posx_pmin+(i-0.5)*pitchx_p
#                  print Xpos, Ypos
                  Xtemp[:,0]=Xpos+Xtemp[:,0]
                  Xtemp[:,2]=Ypos+Xtemp[:,2]
                  Xtemp[:,1]=X_ppos+Xtemp[:,1]
                  Xtemp[:,3]=Y_ppos+Xtemp[:,3]

                  if (i==0) & (j==0):
                     Xnew=Xtemp.copy()
                     Qnew=Qtemp.copy()
#                     print i, j
                  else: 
                     Xnew=np.vstack((Xnew, Xtemp))
                     Qnew=np.hstack((Qnew, Qtemp))
                     
            xxs = Xnew
            q   = Qnew
            print(q)
            N=np.shape(xxs)[0]     
            print(np.shape(xxs))
            continue

        if name=='multibeam_hex_long_y':
            
            print('task: ', name)
            print('shape:',np.shape(xxs))
            pitchx=eval(properties['pitchx']);
            #pitchy=eval(properties['pitchy']);
            pitchy= 1.7320508076*pitchx/2.0 ;
            numx  =eval(properties['numx']);
            numy  =eval(properties['numy']);
            resamp = 0
            
            if 'resample' in properties:
                 resamp =eval(properties['resample'])

            print('resamp=',resamp)
            posxmin=-(numx-1.)/2.0*pitchx;
            posxmax= (numx-1.)/2.0*pitchx;
            posymin=-(numy-1.)/2.0*pitchy;
            posymax= (numy-1.)/2.0*pitchy;
            
#PP this double-loop needs to be improved... 
            for j in range(numy):
               for i in range(numx):
                  if resamp>0:
                     index=random.sample(list(range(0, len(q))), resamp)
                     Xtemp=xxs[index,:]
                     Qtemp=q[index]
                  if resamp==0:   
                     Xtemp=xxs.copy()
                     Qtemp=q.copy()
                  Xpos=posxmin+i*pitchx
                  Ypos=posymin+j*pitchy
                  if j%2==1 :
                      Xpos=posxmin+(i-0.5)*pitchx
#                  print Xpos, Ypos
                  Xtemp[:,0]=Xpos+Xtemp[:,0]
                  Xtemp[:,2]=Ypos+Xtemp[:,2]

                  if (i==0) & (j==0):
                     Xnew=Xtemp.copy()
                     Qnew=Qtemp.copy()
#                     print i, j
                  else: 
                     Xnew=np.vstack((Xnew, Xtemp))
                     Qnew=np.hstack((Qnew, Qtemp))
                     
            xxs = Xnew
            q   = Qnew
            print(q)
            N=np.shape(xxs)[0]     
            print(np.shape(xxs))
            continue



        if name=='write':
            print('task: ', name)
            format=properties['format'];out_file='./'+properties['out_file']
            xpx=xxs_to_xpx(xxs,E_ref)

#            plt.plot (xpx[:,1], xpx[:,3],'.')
#            plt.show()
            
            print('shape before write=', np.shape(xpx))
            if format=='astra':
                z0=eval(properties['z0']);    ct=0;
                F=xpx_to_astra(xpx,q,z0,ct)
                path,file=os.path.split(out_file)
                if not(os.path.exists(path)):os.makedirs(path)
                writetable(out_file,F)

##CRP 11/14/11: CSRtrack can be used in a variety of ways.  The reference particle, by command in the csrtrk.in input deck, can be specified inside the input deck, or read from of the second line in the particle file.  

# This is the original fmt1, reference information specified in the csrtrk.in file.
            if format=='fmt1':
                F=xpx_to_fmt1_noRef(xpx,q,-0.1,-0.1)
#                F=xpx_to_fmt1(xpx,q,-0.1,-0.1,E_ref)
                path,file=os.path.split(out_file)
                if not(os.path.exists(path)):os.makedirs(path)
                writetable(out_file,F)

# This is the modified fmt1, reference information specified in the particle file.
            if format=='fmt1refPart':                
                F=xpx_to_fmt1_withRef(xpx,q,-0.1,-0.1,E_ref)
                path,file=os.path.split(out_file)
                if not(os.path.exists(path)):os.makedirs(path)
                writetable(out_file,F)

#CRP added December 9th,2010
            if format=='elegant':
                z0=eval(properties['z0']);    ct=0;
                F=xxs_to_elegant(xxs,E_ref)
                path,file=os.path.split(out_file)
                if not(os.path.exists(path)):os.makedirs(path)
                writetableELEGANT(out_file,F,E_ref)

#CRP:  impactT is the Impact Z format of output -2, and input 35.   
# NOTE THE CAPITALIZATION
            if format=='impactT':
                z0=eval(properties['z0']);    ct=0;
                F=xpx_to_impactT(xpx,q,z0,ct)
                N=size(F[:,0])
                path,file=os.path.split(out_file)
                if not(os.path.exists(path)):os.makedirs(path)
                writetableIMPACT(out_file,F,N)

            if format=='impactT_noNumber':
                z0=eval(properties['z0']);    ct=0;
                F=xpx_to_impactT(xpx,q,z0,ct)
                N=size(F[:,0])
                path,file=os.path.split(out_file)
                if not(os.path.exists(path)):os.makedirs(path)
                writetableIMPACT_noNumber(out_file,F,N)
#PP original code the writetable was common to all code
#PP this needs to be changed has impactT require the # of macropart on 1st line...
#            path,file=os.path.split(out_file)
#            if not(os.path.exists(path)):os.makedirs(path)
#            writetable(out_file,F)
            continue

        if name=='mix':
            print('task: ', name)
            z1=eval(properties['z1']);z2=eval(properties['z2'])
            zs1=eval(properties['zs1']);
            mix=zeros(N,Float)
            for n in range(N):
                if (xxs[n,5]>zs1)and(xxs[n,4]<z2):
                    mix[n]=1.0;
                    if xxs[n,4]>=z1:
                        z=(z2-xxs[n,4])/(z2-z1)
                        mix[n]=0.5*(1-cos(pi*z))
            mix.shape=(N,1)
            de=xxs[:,5].copy()
            xxs[:,0:6]=xxs[:,0:6]*(1-mix)+xxs2[:,0:6]*mix
            plot(xxs[:,4],de,'r.')
            hold(True);   plot(xxs[:,4],xxs[:,5],'b.')
            if 'axlim' in properties: axis(eval(properties['axlim']))
            fig_file='fig%g.png' %indfig; indfig=indfig+1
            savefig(fig_dir+fig_file);   hold(False);html.write('<IMG SRC="%s">\n' %fig_file)

        if name=='showEmittances':
            print('task: ', name)
            beamMatrix=makeBeamMatrix(xpx)
            emittanceValues=emittances(beamMatrix)
            print(emittanceValues)

    return indfig
# when action=analysis 
# CRP 11/14/11:  Please be careful with this.  action=analysis cannot be combined in parts with action=processor, so the analysis READ function has to  be used to pass particle data between GlueTrack parts.
def Analysis(tasks,fig_dir,html,indfig):
    print("ANALYSIS")
    for i in range(len(tasks)):
        task=tasks[i];  name=task['name']; properties=task['properties']
        if name=='read':
            print('task: ', name)
            in_file=properties['in_file']
            format=properties['format']
            if 'extr_part_1' in properties: 
                extr_part_1=eval(properties['extr_part_1'])
            else:
                extr_part_1=False
            if 'sort_len' in properties: 
                sort_len=eval(properties['sort_len'])
            else:
                sort_len=False
            E_ref=eval(properties['E_ref'])
#            Dh=particle_reader(in_file,format,E_ref,extr_part_1,sort_len)
# Modified to pass properties in to the reader to match improvement to the Processor read function
            Dh=particle_reader(in_file,format,E_ref,extr_part_1,sort_len, properties)
            N=Dh[0]; xpx=Dh[1];  xxs=Dh[2]; q=Dh[3]
            E_ref=Dh[4]; q_tot=Dh[5];  zav=Dh[6]; ct=Dh[7]
            print('N=',N,'  E_ref=',E_ref)
            print('q_tot=',q_tot,'  zav=',zav,'  ct=',ct)
            continue

        if name=='slice_analysis':
            print('task: ', name)
            M_av=eval(properties['M_av'])
            # Q, z, pz, E(pz*q), E(z*q)
            SA0=Slice_Ana_0(xxs,q,E_ref)
            # E(x*q), E(xs*q),E(x*x*q),E(x*xs*q), E(xs*xs*q)
            SAh=Slice_Ana_1(xxs[:,0],xxs[:,1],q) 
            SAv=Slice_Ana_1(xxs[:,2],xxs[:,3],q)
            # averaged parameters
            SA=Slice_Ana_2(SA0,SAh,M_av)
            Z=SA[:,0];  lambdaA=SA[:,1]; epsnx=SA[:,7]
            xa=SA[:,2]; xsa=SA[:,3];    Twx=Twiss(SA)
            SA=Slice_Ana_2(SA0,SAv,M_av); epsny=SA[:,7]
            ya=SA[:,2]; ysa=SA[:,3];    Twy=Twiss(SA)
            SA=Slice_Ana_3(SA0,q,M_av)
            p_av=SA[:,0];  p_rms=SA[:,1];
            print('N=',N,'  E_ref=',E_ref)
            print('q_tot=',q_tot)
            continue




#CRP 12/29/11: Function copied and pasted to analysis.
        if name=='downSample':
#             Down sample a particle distribution, by pruning unsorted lines from xpx and the corresponding lines from q.  At each step, q needs to be increased.
#def downSampleParticles(xpx, q, N_p,targetSize):
            N_p=size(q)
            print(size(q), size(xxs))
            totalQ=sum(q)
#            targetSize=10000
            targetSize=eval(properties['target_size'])
            adjusted_q=totalQ/targetSize
            particle=0
            validParticles=list(range(2,N_p,1))
            numToCut= N_p-targetSize
            rowsToDelete=random.sample(validParticles,numToCut)
            print(len(rowsToDelete))
            xxs=numpy.delete(xxs, rowsToDelete, 0)
            xpx=numpy.delete(xpx, rowsToDelete, 0)
           
            q=numpy.delete(q, rowsToDelete, 0)

            q[:]=adjusted_q
        
            print(size(q), len(q), size(xxs), size (xpx))

            continue




#########return alpha,beta,gamma,epsilon
# CRP 11/07/2011:  Added to record data to user-specified file name
        if name=='slice_analysis_to_file':
            print('task: ', name)
            M_av=eval(properties['M_av'])
            output_root = properties['output_str']
#            full_filename=output_root + '_' AScanNum + '_' BScanNum
            full_filename=output_root
            alphaInit=eval(properties['initial_alpha'])
            betaInit=eval(properties['initial_beta'])
            # Q, z, pz, E(pz*q), E(z*q)
            SA0=Slice_Ana_0(xxs,q,E_ref)
            # E(x*q), E(xs*q),E(x*x*q),E(x*xs*q), E(xs*xs*q)
            SAh=Slice_Ana_1(xxs[:,0],xxs[:,1],q) 
            SAv=Slice_Ana_1(xxs[:,2],xxs[:,3],q)
# CRP 12/17/11
            SAz=Slice_Ana_1(xxs[:,4],xxs[:,5],q)
            # averaged parameters
            SA=Slice_Ana_2(SA0,SAh,M_av)
            Z=SA[:,0];  lambdaA=SA[:,1]; epsnx=SA[:,7]
            xa=SA[:,2]; xsa=SA[:,3];    Twx=Twiss(SA)
            SA=Slice_Ana_2(SA0,SAv,M_av); epsny=SA[:,7]
            ya=SA[:,2]; ysa=SA[:,3];    Twy=Twiss(SA)
# CRP 12/17/11
            SA=Slice_Ana_2(SA0,SAz,M_av); epsnz=SA[:,7]
            za=SA[:,2]; zsa=SA[:,3];    Twz=Twiss(SA)
            SA=Slice_Ana_3(SA0,q,M_av)
            p_av=SA[:,0];  p_rms=SA[:,1];
# It is calculating y, ysa, just not recording them.  I will change this.
#            ZZS=xxs[:,4:6]; 
#            w=where(logical_and(xxs[:,4]<zb,xxs[:,4]>za),1,0)
#            alpha_z,beta_z,gamma_z,epsilon_z=ShowTwissPar(ZZS,w);
            za,zb=eval(properties['interv'])
            XXS=xxs[:,0:2];  YYS=xxs[:,2:4];  ZZS=xxs[:,4:6];
            w=where(logical_and(xxs[:,4]<zb,xxs[:,4]>za),1,0)
            PrintTwissPar(full_filename + '_x', XXS,w)
            PrintTwissPar(full_filename + '_y', YYS,w)
            PrintTwissPar(full_filename + '_z', ZZS,w)
            f=open(full_filename + '_x_initial','w');
            f.write(str(alphaInit))
            f.write("\n") 
            f.write(str(betaInit)) 
            f.close()


            energy=E_ref
            momentumSpread=std(xpx[:,5])
            f=open(full_filename + '_energy','w');
            f.write('%s %s \n' % (energy, momentumSpread))
            f.close()

#            f=open(full_filename + '_y_initial','w');
#            f.write(str(alphaInit))
#            f.write("\n") 
#            f.write(str(betaInit)) 
#            f.close()

            print('N=',N,'  E_ref=',E_ref)
            print('q_tot=',q_tot)
            continue

        if name=='m_current':
            print('task: ', name)
            sigmaI=eval(properties['sigmaI'])
            gf=Gauss_Filter(lambdaA,Z,sigmaI)
            plot(Z*1e3,lambdaA*c0*1e-3,gf[:,0]*1e3,c0*gf[:,1]*1e-3)
            if 'xlim' in properties:  xlim(eval(properties['xlim']))
            if 'ylim' in properties:  ylim(eval(properties['ylim']))
            if 'savefile' in properties:
                savefile=(properties['savefile'])
                writevectors(savefile,gf[:,0]*1e3,c0*gf[:,1]*1e-3)
            xlabel('z[mm]');ylabel('I[kA]');title('current')
            fig_file='fig%g.png' %indfig; indfig=indfig+1
            savefig(fig_dir+fig_file); hold(False);  html.write('<IMG SRC="%s">\n' % fig_file)
            continue



        if name=='m_emmitance':
            print('task: ', name)
            m1=max(epsnx*1e6);m2=max(epsny*1e6);m=max([m1,m2]); koef=m/max(lambdaA)
            if 'savefile' in properties:
                savefile=(properties['savefile'])
                writevectors(savefile,Z*1e3,epsnx*1e6,Z*1e3,epsny*1e6,Z*1e3,koef*lambdaA)
            plot(Z*1e3,epsnx*1e6,Z*1e3,epsny*1e6,Z*1e3,koef*lambdaA)
            if 'xlim' in properties:  xlim(eval(properties['xlim']))
            if 'ylim' in properties:  ylim(eval(properties['ylim']))
            xlabel('z[mm]');title('norm. emit. hor&vert')
            fig_file='fig%g.png' %indfig; indfig=indfig+1
            savefig(fig_dir+fig_file); hold(False);  html.write('<IMG SRC="%s">\n' % fig_file)
            continue

        if name=='m_momentum':
            print('task: ', name)
            if 'savefile' in properties:
                savefile=(properties['savefile'])
                writevectors(savefile,xxs[:,4]*1e3,xpx[:,5]*1e-6,p_av*1e-6)
            plot(xxs[:,4]*1e3,xpx[:,5]*1e-6,'.b',markersize=1)
            hold(True);
            plot(xxs[:,4]*1e3,p_av*1e-6,'r')
            if 'xlim' in properties:  xlim(eval(properties['xlim']))
            if 'ylim' in properties:  ylim(eval(properties['ylim']))
            xlabel('z[mm]');ylabel('MeV/c');title('aver. of momentum')
            fig_file='fig%g.png' %indfig; indfig=indfig+1
            savefig(fig_dir+fig_file); hold(False);  html.write('<IMG SRC="%s">\n' % fig_file)
            continue

        if name=='m_rms_momentum':
            print('task: ', name)
            m=max(p_rms*1e-6);koef=m/max(lambdaA)
            if 'savefile' in properties:
                savefile=(properties['savefile'])
                writevectors(savefile,Z*1e3,p_rms*1e-6,lambdaA*koef)
            plot(Z*1e3,p_rms*1e-6,Z*1e3,lambdaA*koef)
            if 'xlim' in properties:  xlim(eval(properties['xlim']))
            if 'ylim' in properties:  ylim(eval(properties['ylim']))
            xlabel('z[mm]');ylabel('MeV/c');title('rms of momentum')
            fig_file='fig%g.png' %indfig; indfig=indfig+1
            savefig(fig_dir+fig_file); hold(False);  html.write('<IMG SRC="%s">\n' % fig_file)
            continue


#CRP 11/22/11:  Created "combo" slice analysis, that combines the write functions in to a single file.  Not really sure what the Koefficients are for, but they are saved here:  Emittance and P use the same one, while Twiss X and Twiss Y's are different. Also, added a modifed "slice current" or "slice charge" column, to better identify the bunch center.  
#SortedData=np.sort(currentParticleData, order='h')
# I want to sort, then associate particles with the nearest z-bin.
        if name=='m_combo_write':
            print('task: ', name)
            m1=max(epsnx*1e6);m2=max(epsny*1e6);mEmit=max([m1,m2]); koef=mEmit/max(lambdaA)
            m1=max(Twx[:,0]);m2=max(Twx[:,1]);m=max([m1,m2]); koefx=m/max(lambdaA)
            m1=max(Twy[:,0]);m2=max(Twy[:,1]);m=max([m1,m2]); koefy=m/max(lambdaA)
            m=max(p_rms*1e-6);koef=m/max(lambdaA)

            sigmaI=eval(properties['sigmaI'])
            print('gauss filter')
            gf=Gauss_Filter(lambdaA,Z,sigmaI)
            print('gauss filter completed')
# gf[:,0]*1e3,c0*gf[:,1]*1e-3    This is what gets recorded.  I need to adapt it to the other format.  Just copy it to surrounding Z values, I guess?


            Zrms=std(Z)
            BinSize=Zrms*0o1
            NumPoints=len(Z)
            cms=2.998e8
            BunchCharge=sum(q)
            ParticleDensity=zeros(NumPoints)
            ChargePerMacro=q[1]
            N=0
            M=0
            gLength=len(gf)
            adjustedCurrents=zeros(NumPoints)

# This is necessary to adapt the current calculation to the same Z locations
# Original design was wrong.
            print('current matching')
#            while N<NumPoints:
#                    nearestIndex=(numpy.abs(gf[:,0]-Z[N])).argmin()
#                    adjustedCurrents[N]=gf[nearestIndex,1]
#####                    M=0
#                    N=N+1


# This version is much faster and is of the order N, instead of N^2 or worse.
            while N<NumPoints:
                    while Z[N]>gf[M,0]:
                        if M<gLength:
                            M=M+1
                    if Z[N]<gf[M,0]:
                        adjustedCurrents[N]=gf[M,1]
                    N=N+1


            print('current matching completed')


            


#sliceNeg1=vstack((sliceNeg1, currentParticleData[:][line]))

            print(sum(q))
            if 'savefile' in properties:
                savefile=(properties['savefile'])
                writevectors(savefile,Z*1e3,epsnx*1e6,epsny*1e6,epsnz*1e6,koef*lambdaA,p_av*1e-6,p_rms*1e-6,Twx[:,0],Twx[:,1],lambdaA*koefx,Twy[:,0],Twy[:,1],Twz[:,0],Twz[:,1], lambdaA*koefy, adjustedCurrents*c0,xa*1e3,xsa*1e3,ya*1e3,ysa*1e3,za*1e3,zsa*1e3)
            #Z, emitX, emitY, general Koef, RMS momentum, alphax, betax, koefx, alphay, betay, koefy

            continue



        if name=='m_hor_offset':
            print('task: ', name)
            m=max(xa*1e3);koef=m/max(lambdaA)
            if 'savefile' in properties:
                savefile=(properties['savefile'])
                writevectors(savefile,Z*1e3,xa*1e3,lambdaA*koef)
            plot(Z*1e3,xa*1e3,Z*1e3,lambdaA*koef)
            if 'xlim' in properties:  xlim(eval(properties['xlim']))
            if 'ylim' in properties:  ylim(eval(properties['ylim']))
            xlabel('z[mm]');ylabel('mm');title('horizontal offset')
            fig_file='fig%g.png' %indfig; indfig=indfig+1
            savefig(fig_dir+fig_file); hold(False);  html.write('<IMG SRC="%s">\n' % fig_file)
            continue
    
        if name=='m_hor_slope':
            print('task: ', name)
            m=max(xsa*1e3);koef=m/max(lambdaA)
            if 'savefile' in properties:
                savefile=(properties['savefile'])
                writevectors(savefile,Z*1e3,xsa*1e3,lambdaA*koef)
            plot(Z*1e3,xsa*1e3,Z*1e3,lambdaA*koef)
            if 'xlim' in properties:  xlim(eval(properties['xlim']))
            if 'ylim' in properties:  ylim(eval(properties['ylim']))
            xlabel('z[mm]');ylabel('mrad');title('horizontal slope')
            fig_file='fig%g.png' %indfig; indfig=indfig+1
            savefig(fig_dir+fig_file); hold(False);  html.write('<IMG SRC="%s">\n' % fig_file)
            continue

        if name=='m_Twissx':
            print('task: ', name)
            m1=max(Twx[:,0]);m2=max(Twx[:,1]);m=max([m1,m2]); koef=m/max(lambdaA)
            if 'savefile' in properties:
                savefile=(properties['savefile'])
                writevectors(savefile,Z*1e3,Twx[:,0],Twx[:,1],lambdaA*koef)
            plot(Z*1e3,Twx[:,0],Z*1e3,Twx[:,1],Z*1e3,lambdaA*koef)
            if 'xlim' in properties:  xlim(eval(properties['xlim']))
            if 'ylim' in properties:  ylim(eval(properties['ylim']))
            xlabel('z[mm]');ylabel('');title('alpha, beta')
            fig_file='fig%g.png' %indfig; indfig=indfig+1
            savefig(fig_dir+fig_file); hold(False);  html.write('<IMG SRC="%s">\n' % fig_file)
            continue

        if name=='m_Twissy':
            print('task: ', name)
            m1=max(Twy[:,0]);m2=max(Twy[:,1]);m=max([m1,m2]); koef=m/max(lambdaA)
            if 'savefile' in properties:
                savefile=(properties['savefile'])
                writevectors(savefile,Z*1e3,Twy[:,0],Twy[:,1],lambdaA*koef)
            plot(Z*1e3,Twy[:,0],Z*1e3,Twy[:,1],Z*1e3,lambdaA*koef)
            if 'xlim' in properties:  xlim(eval(properties['xlim']))
            if 'ylim' in properties:  ylim(eval(properties['ylim']))
            xlabel('z[mm]');ylabel('');title('alpha, beta')
            fig_file='fig%g.png' %indfig; indfig=indfig+1
            savefig(fig_dir+fig_file); hold(False);  html.write('<IMG SRC="%s">\n' % fig_file)
            continue

        if name=='m_topview':
            print('task: ', name)
            if 'savefile' in properties:
                savefile=(properties['savefile'])
                writevectors(savefile,xxs[:,4]*1e3,xpx[:,0]*1e3)
            plot(xxs[:,4]*1e3,xpx[:,0]*1e3,'.r',markersize=1);
            if 'xlim' in properties:  xlim(eval(properties['xlim']))
            if 'ylim' in properties:  ylim(eval(properties['ylim']))
            xlabel('z[mm]');ylabel('');title('top view')
            fig_file='fig%g.png' %indfig; indfig=indfig+1
            savefig(fig_dir+fig_file); hold(False);  html.write('<IMG SRC="%s">\n' % fig_file)



# CRP: 12/14/11
# Takes loaded xpx format with a user-specified sigma_long, and both plots the cummulative current profile, as well as a fourier frequency analysis
        if name=='fourier_plot':
                lowWavelength=.1e-9
                highWavelength=10e-6
                NumWavelengths=2000
                waveStep=(highWavelength-lowWavelength)/(NumWavelengths-1)
                totalCharge=sum(q)
                wavelengthVector=zeros((NumWavelengths,1))
                fourier=zeros((NumWavelengths,1))
                BFF=zeros((NumWavelengths,1))

                wavelengthList=np.arange(0,NumWavelengths+.01,1)
                wavelengthList=.1e-8**(wavelengthList/NumWavelengths)
                
#                expectedMacroSize=1e-9
                expectedMacroSize=totalCharge


                waveInd=0
                while waveInd < NumWavelengths:
                #        wavelength=lowWavelength+waveInd*waveStep
                        wavelength=wavelengthList[waveInd]
                        wavelengthVector[waveInd]=wavelength

#                        print len(fourier), size(q), len(q), len(xxs), size(xxs)
                        value0=exp(-expectedMacroSize**2.0*4.0*pi**2/wavelength**2.0)
                        value4=-expectedMacroSize**2.0*4.0*pi**2/wavelength**2.0
                        value1=sum(q[:]*cos(2.0*pi*xxs[:,4]/wavelength) )
                        value2=sum(q[:]*sin(2.0*pi*xxs[:,4]/wavelength) )
#                        print totalCharge, expectedMacroSize, wavelength, value0, value1, value2, value4
                        fourier[waveInd]=(1.0/totalCharge)*exp(-expectedMacroSize**2.0*4.0*pi*pi/wavelength**2.0)*(sum(q[:]*sin(2.0*pi*xxs[:,4]/wavelength)) + sum(q[:]*cos(2.0*pi*xxs[:,4]/wavelength) )  )
                        scalingFactor=1/totalCharge**2
                        BFF[waveInd]=scalingFactor*np.exp(-expectedMacroSize**2*4*np.pi*np.pi/wavelength**2.0)*(sum(q[:]**2.0*np.sin(2.0*np.pi*xxs[:,4]/wavelength))**2 + sum(q[:]**2.0*np.cos(2.0*np.pi*xxs[:,4]/wavelength) )**2.0  )
#                        print fourier[waveInd], BFF[waveInd]
                        waveInd=waveInd+1

                print('fourier done')
                xbins=np.arange(-2,2,.01)*std(xxs[:,4])
                binSize=xbins[2]-xbins[1]

                Ix=0*xbins

                histInd=0

                numParticles=len(xxs)
                print(len(xbins), len(Ix), numParticles)

#                while histInd<numParticles-1:
#                        Ix= Ix + xxs[histInd,4]*binSize*exp(-(xbins-xxs[histInd,4])**2/(2*expectedMacroSize**2))/sqrt(2*pi*expectedMacroSize**2);
#                print histInd
#                histInd=histInd+1  

                plt.figure()
                semilogx(wavelengthVector,fourier)
                plt.title('fourier')


                plt.figure()
                loglog(wavelengthVector,BFF)
                plt.title('BFF')

#                figure()
#                plot(xbins,Ix)
#                title('current distribution')

                plt.show()



                continue



# CRP: 12/29/11
# Takes loaded xpx format with a user-specified sigma_long, and both plots the cummulative current profile, as well as a fourier frequency analysis
# THIS IS WHAT IS ACTUALLY USED FOR CSRTRACK!!!
        if name=='gaussian_fourier_plot':
                lowWavelength=.1e-9
                highWavelength=10e-6
                NumWavelengths=2000
                gLength=eval(properties['sigma_long'])

                waveStep=(highWavelength-lowWavelength)/(NumWavelengths-1)
                totalCharge=sum(q)
                wavelengthVector=zeros((NumWavelengths,1))
                fourier=zeros((NumWavelengths,1))
                BFF=zeros((NumWavelengths,1))

                wavelengthList=np.arange(0,NumWavelengths+.01,1)
                wavelengthList=.1e-8**(wavelengthList/NumWavelengths)
                
#                expectedMacroSize=1e-9
                expectedMacroSize=gLength
                numParticles=len(xxs)
# Guassian Breakdown
                minZ=min(xxs[:,4])
                maxZ=max(xxs[:,4])
                
                numBins=2000
                zbins=arange(minZ,maxZ,(maxZ-minZ)/numBins)
                binSize=zbins[2]-zbins[1]
                Iz=0.0*zbins
                print(binSize)

                histInd=0
                while histInd<numParticles-1:
                        Iz= Iz + binSize*exp(-(zbins-xxs[histInd,4])**2/(2*gLength**2))/sqrt(2*pi*gLength**2);
#                        print histInd
                        histInd=histInd+1  

#                plt.figure()
#                plot(zbins, Iz)
#                plt.show()
#                scalingFactor=totalCharge**2

                waveInd=0
                while waveInd < NumWavelengths:
                        wavelength=wavelengthList[waveInd]
                        wavelengthVector[waveInd]=wavelength
                        fourier[waveInd]=(1.0/totalCharge)*exp(-expectedMacroSize**2.0*4.0*pi*pi/wavelength**2.0)*(sum(Iz[:]*sin(2.0*pi*zbins[:]/wavelength)) + sum(Iz[:]*cos(2.0*pi*zbins[:]/wavelength) )  )
                        BFF[waveInd]=(1.0/totalCharge**2)*np.exp(-expectedMacroSize**2*4*np.pi*np.pi/wavelength**2.0)*(sum(Iz[:]**2*np.sin(2.0*np.pi*zbins[:]/wavelength))**2 + sum(Iz[:]**2*np.cos(2.0*np.pi*zbins[:]/wavelength) )**2.0  )
#                        print fourier[waveInd], BFF[waveInd]
                        waveInd=waveInd+1

                print('fourier done')
                xbins=np.arange(-2,2,.01)*std(xxs[:,4])
                binSize=xbins[2]-xbins[1]

                Ix=0*xbins

                histInd=0

                numParticles=len(xxs)
                print(len(xbins), len(Ix), numParticles)


                plt.figure()
                plot(zbins, Iz)
                plt.title('current')

                plt.figure()
                semilogx(wavelengthVector,fourier)
                plt.title('fourier')


                plt.figure()
                loglog(wavelengthVector,BFF)
                plt.title('BFF')
                xlabel('wavelength (m)')


                cms=2.998e8


                plt.show()



                continue



# CRP: 12/29/11
# Takes loaded xpx format with a user-specified sigma_long, and both plots the cummulative current profile, as well as a fourier frequency analysis
# THIS IS WHAT IS ACTUALLY USED FOR CSRTRACK!!!
        if name=='gaussian_3D_fourier_plot':
                lowWavelength=.1e-9
                highWavelength=1e-8
                NumWavelengths=20000
                gLength=eval(properties['sigma_long'])
                gHor=eval(properties['sigma_hor'])
                gVer=eval(properties['sigma_ver'])

# Need to 

                lambdasX, fX, fourierCurveX, BFFCurveX, xbins, Ix=fourier_analysis(xxs[:,0],q, lowWavelength, highWavelength,NumWavelengths, gHor)
                lambdasY, fY, fourierCurveY, BFFCurveY, ybins, Iy=fourier_analysis(xxs[:,2],q, lowWavelength, highWavelength,NumWavelengths, gVer)
                lambdasZ, fZ, fourierCurveZ, BFFCurveZ, zbins, Iz=fourier_analysis(xxs[:,4],q, lowWavelength, highWavelength,NumWavelengths, gLength)

                plt.figure()
                subplot(221)
                semilogx(fX,fourierCurveX)
                subplot(222)
                semilogx(fY,fourierCurveY)
                subplot(223)
                semilogx(fZ,fourierCurveZ)
                title('Fourier')

                plt.figure()
                subplot(221)
                loglog(lambdasX,BFFCurveX)
                subplot(222)
                loglog(lambdasY,BFFCurveY)
                subplot(223)
                loglog(lambdasZ,BFFCurveZ)
                title('BFF')

                plt.figure()
                subplot(221)
                plot(xbins,Ix)
                subplot(222)
                plot(ybins,Iy)
                subplot(223)
                plot(zbins,Iz)
                title('Imposed Distributions')


                show()
                continue



#CRP: 02/01/12
# Break the Fourier content in to a function in particle_tools, to make it work with both 1D and 3D.

# CRP: 11/21/11
# I should make an m_ALL file that writes all of these to one file, or at least the ones that have the same length (column Z)


        if name=='write':
            print('task: ', name)
            format=properties['format'];out_file='./'+properties['out_file']
            if format=='astra':
                xpx=xxs_to_xpx(xxs,E_ref)
                z0=eval(properties['z0']);    ct=0;
                F=xpx_to_astra(xpx,q,z0,ct)
            if format=='fmt1':
                xpx=xxs_to_xpx(xxs,E_ref)
#                F=xpx_to_fmt1(xpx,q,-0.1,-0.1)
                F=xpx_to_fmt1(xpx,q,-0.1,-0.1,E_ref)
            if format=='genesis':
                sigmaI=eval(properties['sigmaI'])
                betax,betay=eval(properties['beta'])
                nstep=eval(properties['nstep'])
                qf=Gauss_Filter(lambdaA,Z,sigmaI)
                ZF=qf[:,0]
                Nz=size(ZF)
                if 'interv' in properties:
                    start,stop=eval(properties['interv'])
                    z0=ZF[0]; nstart=nstop=0;
                    while  z0<start or nstart+10>Nz-1:
                        nstart=nstart+nstep;  z0=ZF[nstart]
                    nstop=nstart
                    while  z0<stop or nstop+nstep>Nz-1:
                        nstop=nstop+nstep;  z0=ZF[nstop]
                else:
                    nstart=1; nstop=Nz-1
                Nz=size(list(range(nstart,nstop,nstep)));
                F=zeros((Nz,8),Float)
                F[:,0]=ZF[nstart:nstop:nstep]-ZF[nstart];     #z
                point=Gauss_Filter(p_av,Z,sigmaI);
                F[:,1]=point[nstart:nstop:nstep,1]/E_ele_eV;   #gamma
                point=Gauss_Filter(p_rms,Z,sigmaI);
                F[:,2]=point[nstart:nstop:nstep,1]/E_ele_eV;  #delgam
                point=Gauss_Filter(epsnx,Z,sigmaI);
                F[:,3]=abs(point[nstart:nstop:nstep,1]);      #emitx
                point=Gauss_Filter(epsny,Z,sigmaI);
                F[:,4]=abs(point[nstart:nstop:nstep,1]);      #emity
                F[:,5]=sqrt(F[:,3]*betax/F[:,1]);           #rxbeam
                F[:,6]=sqrt(F[:,4]*betay/F[:,1]);           #rybeam
                F[:,7]=qf[nstart:nstop:nstep,1]*c0       #I
                z=F[:,0];y=F[:,1]; curr=F[:,7]
                m2=max(y);m1=min(y)
                koef=(m2-m1)/max(abs(curr));
                title('GAMMA0');plot(z,m1+curr*koef,z,y);
                fig_file='fig%g.png' %indfig; indfig=indfig+1
                savefig(fig_dir+fig_file); hold(False);  html.write('<IMG SRC="%s">\n' % fig_file)
                y=F[:,2];        koef=max(abs(y))/max(abs(curr));
                title('DELGAMMA');plot(z,curr*koef,z,y);
                fig_file='fig%g.png' %indfig; indfig=indfig+1
                savefig(fig_dir+fig_file); hold(False);  html.write('<IMG SRC="%s">\n' % fig_file)
                y=F[:,3];koef=max(abs(y))/max(abs(curr));
                title('EMITX');plot(z,curr*koef,z,y);
                fig_file='fig%g.png' %indfig; indfig=indfig+1
                savefig(fig_dir+fig_file); hold(False);  html.write('<IMG SRC="%s">\n' % fig_file)
                y=F[:,4];koef=max(abs(y))/max(abs(curr));
                title('EMITY');plot(z,curr*koef,z,y);
                fig_file='fig%g.png' %indfig; indfig=indfig+1
                savefig(fig_dir+fig_file); hold(False);  html.write('<IMG SRC="%s">\n' % fig_file)
                y=F[:,5];koef=max(abs(y))/max(abs(curr));
                title('RX');plot(z,curr*koef,z,y);
                fig_file='fig%g.png' %indfig; indfig=indfig+1
                savefig(fig_dir+fig_file); hold(False);  html.write('<IMG SRC="%s">\n' % fig_file)
                y=F[:,6];koef=max(abs(y))/max(abs(curr));
                title('RY');plot(z,curr*koef,z,y);
                fig_file='fig%g.png' %indfig; indfig=indfig+1
                savefig(fig_dir+fig_file); hold(False);  html.write('<IMG SRC="%s">\n' % fig_file)
                
                head='#\n? VERSION = 1.0\n? SIZE = %s\n? COLUMNS ZPOS GAMMA0 DELGAM EMITX EMITY RXBEAM RYBEAM CURPEAK\n' %Nz
            writetable(out_file,F,head)
            continue
    return indfig

def WriteTask(out_file,actions):
    f=open(out_file,'w')
    for i in range(len(actions)):
        action=actions[i];
        name=action['name'];     tasks=action['tasks']
        comment=action['comment']+'\n'
        f.write(comment)
        str='action='+name+'\n'; f.write(str)
        for j in range(len(tasks)):
            task=tasks[j];  name=task['name']; properties=task['properties']
            str='\ttask='+name+'\n'; f.write(str)
            for mkey in list(properties.keys()):
                str='\t\t'+mkey+'=%s\n'; f.write(str %properties[mkey])
        f.write('\n')
    f.close()

def ReadTask(in_file):
    f=open(in_file,'r')
    actions=[];actionname='';
    tasks=[];taskname='';
    properties={};

    for line in f:
        line=line.replace('\n','');line=line.replace('\t',' ')
        pos=line.find('!');
        if pos>=0: line=line[0:pos]
        if len(line)==0: continue
        if line[0]=='#':
            if actionname!='':
                if taskname!='':
                    task={'name':taskname,'properties':properties}
                    tasks.append(task)
                    properties={};taskname=''
                action={'name':actionname,'tasks':tasks,'comment':comment}
                actions.append(action)
                tasks=[]
            comment=line; continue
        line=line.replace(' ','')
        if line.find('action')>=0:
#pp            cols=line.split('=')
            colt=line.split('\r')
            cols=colt[0].split('=')
            actionname=cols[1]
#            print actionname
            continue
        if line.find('task')>=0:
            if taskname!='':
                task={'name':taskname,'properties':properties}
                tasks.append(task)
                properties={}
#pp            cols=line.split('=')
            colt=line.split('\r')
            cols=colt[0].split('=')
            taskname=cols[1]; continue
#pp            print taskname; continue
        if line.find('=')>=0:
#pp            cols=line.split('=')
            colt=line.split('\r')
            cols=colt[0].split('=')
            property=cols[0]
            value=cols[1];
            properties[property]=value
            

    if actionname!='':
        if taskname!='':
            task={'name':taskname,'properties':properties}
            tasks.append(task)
        action={'name':actionname,'tasks':tasks,'comment':comment}
        actions.append(action)
        print("pp in ReadTask", action)
    f.close()
    return actions

def ASTRA(tasks,fig_dir,html,indfig):
    for i in range(len(tasks)):
        task=tasks[i];  name=task['name']; properties=task['properties']
        if name=='run':
            print('task: ', name)
            if 'command' in properties: 
                command=properties['command']
                work_dir='./'+properties['work_dir']
                command=command.replace(',',' ')
                os.chdir(work_dir)
                src='../Codes/ASTRA'
                copyfiles(src,'./')
                names = os.listdir(src)
                print(command)
                os.system(command)
                for name in names:
                    os.remove(name)
                os.chdir('..')
                continue
        if name=='os':
            print('task: ', name)
            if 'command' in properties: 
                command=properties['command']
                command=command.replace(',',' ')
                print(command)
                os.system(command)
                continue
        if name=='adjust_acceleration':
            print('task: ', name)
            in_file=properties['in_file']
            out_file=properties['out_file']
            V0,phi0,V1,phi1,omega=eval(properties['V0_phi0_V1_phi1_omega'])
            omega=2*pi*omega; phi0=phi0*pi/180; phi1=phi1*pi/180;
            AdjustAcceleration(in_file,out_file,V0,phi0,V1,phi1,omega)
    return indfig


### This is so messed up (June 28th, 2012
def ELEGANT(tasks,fig_dir,html,indfig):
    for i in range(len(tasks)):
        task=tasks[i];  name=task['name']; properties=task['properties']
        if name=='run':
            print('task: ', name)
            if 'command' in properties: 
                os.system ('pwd')

                work_dir=properties['work_dir']
                deck_root=properties['deck_inp']
                deck_inp=properties['deck_inp']+'.ele'

                dist_inp=properties['dist_inp']
                lattice_inp=properties['lattice_inp']

                deck_std=properties['work_dir']+'/elegantDeck.ele'
                dist_std=properties['work_dir']+'/' + properties['dist_inp_rename']

                lattice_std=properties['work_dir']+'/'

                command_line=properties['command'] + ' elegantDeck.ele'

                print(command_line)
                print('cp %s %s' % (deck_inp, deck_std))
                print('cp %s %s' % (dist_inp, dist_std))
                print('cp %s %s' % (lattice_inp, lattice_std))
                os.system ('pwd')
                os.system ('cp %s %s' % (deck_inp, deck_std))
                os.system ('cp %s %s' % (dist_inp, dist_std))
                os.system ('cp %s %s' % (lattice_inp, lattice_std))

                os.chdir(work_dir)  
                os.system (command_line)
#                for fname in glob.glob('fort.*'):
#                   print fname
#                   fnameout='.' + deck_roo+fname[4:]
#                   print fnameout
#
#                   os.rename(fname, fnameout)
#
#                   
#                os.system ('mv partcl.data lastDist.data')
                os.chdir('..')                
                continue
#        if name=='os':
#            print 'task: ', name
#            if 'command' in properties: 
#                command=properties['command']
#                command=command.replace(',',' ')
#                print command
#                os.system(command)
#                continue
#    return indfig
        if name=='adjust_acceleration':
            print('task: ', name)
            in_file=properties['in_file']
            out_file=properties['out_file']
            V0,phi0,V1,phi1,omega=eval(properties['V0_phi0_V1_phi1_omega'])
            omega=2*pi*omega; phi0=phi0*pi/180; phi1=phi1*pi/180;
            AdjustAcceleration(in_file,out_file,V0,phi0,V1,phi1,omega)
    return indfig



def IMPACTT(tasks,fig_dir,html,indfig):
    for i in range(len(tasks)):
        task=tasks[i];  name=task['name']; properties=task['properties']
        if name=='run':
            print('task: ', name)
            if 'command' in properties: 
                command=properties['command']
                work_dir='./'+properties['work_dir']
                command=command.replace(',',' ')
                os.chdir(work_dir)
                src='/opt/nicadd/contrib/piot/Impact_T/impactT_v1.6_32'
                copyfiles(src,'./')
                names = os.listdir(src)
                print(command)
                os.system(command)
                for name in names:
                    os.remove(name)
                os.chdir('..')
                continue
            if 'input' in properties: 
                input=properties['input']
                work_dir='./'+properties['work_dir']
                command=command.replace(',',' ')
                os.chdir(work_dir)
                src='/opt/nicadd/contrib/piot/Impact_T/impactT_v1.6_32'
                copyfiles(src,'./')
                names = os.listdir(src)
                print(command)
                os.system(command)
                for name in names:
                    os.remove(name)
                os.chdir('..')
                continue
        if name=='os':
            print('task: ', name)
            if 'command' in properties: 
                command=properties['command']
                command=command.replace(',',' ')
                print(command)
                os.system(command)
                continue
    return indfig


# September 7th, 2012:  CRP:  Now that I need to start working with Impact-T, I've found the original function to be completely inadequate.  For now, I am copying the Impact-Z functionality, of copying file names, instead of using Impact-T's mandatory defaults.
def IMPACTT_Copy(tasks,fig_dir,html,indfig):
    for i in range(len(tasks)):
        task=tasks[i];  name=task['name']; properties=task['properties']
        if name=='run':
            print('task: ', name)
            if 'command' in properties: 
                os.system ('pwd')
                command=properties['command']
                work_dir=properties['work_dir']
                deck_roo=properties['deck_inp']
                deck_inp=properties['deck_inp']+'.in'

# CRP Change:   Want to grab out of ANY directory, not one from ImpactZ.
#                dist_inp=work_dir+properties['dist_inp']
                dist_inp=properties['dist_inp']
        #PP standard input file required by Impact-T
                deck_std=properties['work_dir']+'/ImpactT.in'
                dist_std=properties['work_dir']+'/partcl.data'
        #CRP I need to make it so GlueTrack copies to FILENAME.18, etc., rather than the fort.18 format it uses by default.
#                os.system ('pwd')
#                os.chdir(work_dir)
#                os.system ('pwd')
                print(command)
                print(deck_inp)
                print(dist_inp)
                os.system ('pwd')
                os.system ('cp %s %s' % (deck_inp, deck_std))
                os.system ('cp %s %s' % (dist_inp, dist_std))

#Needs to be before OS system?
                os.chdir(work_dir) 
                os.system ('pwd')


# should it copy ALL rfdatas, and ALL T7 datas?  This is not at all clear how I should do this.
#### Testing purposes
#                os.system ('ls *.T7')



                os.system (command)
#                os.chdir('..')        
#                for fname in glob.glob(work_dir + 'fort.*'):
                for fname in glob.glob('fort.*'):
                   print(fname)
# CRP adding the dot?  no, the 4: should handle that
#
                   fnameout='.' + deck_roo+fname[4:]
                   print(fnameout)
# This line is giving errors, not reading the things as strings.

#This is the original
#                   os.system ('cp %s %s' % (fname, fnameout))                   
#                   os.system ('cp %s %s' % (fname, fnameout))
                   os.rename(fname, fnameout)
#                   os.rename(fname, '.' + deck_roo+fname[4:])
                   
                os.system ('mv partcl.data lastDist.data')
                os.chdir('..')                
                continue


            if 'input' in properties: 
                input=properties['input']
                work_dir='./'+properties['work_dir']
                command=command.replace(',',' ')
                os.chdir(work_dir)
                src='/opt/nicadd/contrib/piot/Impact_T/impactT_v1.6_32'
                copyfiles(src,'./')
                names = os.listdir(src)
                print(command)
                os.system(command)

# CRP Note 01/18/12:  Is it really necessary to remove the old executables?
                for name in names:
                    os.remove(name)
                os.chdir('..')
                continue
        if name=='os':
            print('task: ', name)
            if 'command' in properties: 
                command=properties['command']
                command=command.replace(',',' ')
                print(command)
                os.system(command)
                continue
    return indfig







#PP 02/25/2011:
def IMPACTZ(tasks,fig_dir,html,indfig):
    for i in range(len(tasks)):
        task=tasks[i];  name=task['name']; properties=task['properties']
        if name=='run':
            print('task: ', name)
            if 'command' in properties: 
                os.system ('pwd')
                command=properties['command']
                work_dir=properties['work_dir']
                deck_roo=properties['deck_inp']
                deck_inp=properties['deck_inp']+'.in'

# CRP Change:   Want to grab out of ANY directory, not one from ImpactZ.
#                dist_inp=work_dir+properties['dist_inp']
                dist_inp=properties['dist_inp']
        #PP standard input file required by Impact-Z
                deck_std=properties['work_dir']+'/test.in'
                dist_std=properties['work_dir']+'/partcl.data'
        #CRP I need to make it so GlueTrack copies to FILENAME.18, etc., rather than the fort.18 format it uses by default.
#                os.system ('pwd')
#                os.chdir(work_dir)
#                os.system ('pwd')
                print(command)
                print(deck_inp)
                print(dist_inp)
                os.system ('pwd')
                os.system ('cp %s %s' % (deck_inp, deck_std))
                os.system ('cp %s %s' % (dist_inp, dist_std))

#Needs to be before OS system?
                os.chdir(work_dir) 
                os.system ('pwd')

#### Testing purposes
                os.system ('ls *.T7')



                os.system (command)
#                os.chdir('..')        
#                for fname in glob.glob(work_dir + 'fort.*'):
                for fname in glob.glob('fort.*'):
                   print(fname)
# CRP adding the dot?  no, the 4: should handle that
#
                   fnameout='.' + deck_roo+fname[4:]
                   print(fnameout)
# This line is giving errors, not reading the things as strings.

#This is the original
#                   os.system ('cp %s %s' % (fname, fnameout))                   
#                   os.system ('cp %s %s' % (fname, fnameout))
                   os.rename(fname, fnameout)
#                   os.rename(fname, '.' + deck_roo+fname[4:])
                   
                os.system ('mv partcl.data lastDist.data')
                os.chdir('..')                
                continue


            if 'input' in properties: 
                input=properties['input']
                work_dir='./'+properties['work_dir']
                command=command.replace(',',' ')
                os.chdir(work_dir)
                src='/opt/nicadd/contrib/piot/Impact_Z/impactZ_32builtbeta'
                copyfiles(src,'./')
                names = os.listdir(src)
                print(command)
                os.system(command)

# CRP Note 01/18/12:  Is it really necessary to remove the old executables?
                for name in names:
                    os.remove(name)
                os.chdir('..')
                continue
        if name=='os':
            print('task: ', name)
            if 'command' in properties: 
                command=properties['command']
                command=command.replace(',',' ')
                print(command)
                os.system(command)
                continue
    return indfig




### Modified CRP  08/18/11
# Need to simplify and get its copying working, due to CSRtrack's rather finicky input requirements.
# Needs to auto copy input file to csrtrk.in
# Need to somehow save the output files.
# I WANT to replace the input distribution name in the .in file, as there is no command line for this.
# Updated CRP 02/21/12 to use new executable
def CSRTrack(tasks,fig_dir,html,indfig):
    for i in range(len(tasks)):
        task=tasks[i];  name=task['name']; properties=task['properties']
        if name=='run':
            print('task: ', name)
            if 'command' in properties: 
                command=properties['command']
                print(command)
                os.system('pwd')
                work_dir='./'+properties['work_dir']
                dist_inp='./'+properties['dist_inp']
                deck_inp='./'+properties['deck_inp']
                dist_std=work_dir+'/in/genericGlueTrackFileName.fmt1'
                command=command.replace(',',' ')
                os.system('pwd')
                os.system ('cp %s %s' % (dist_inp, dist_std))
                os.chdir(work_dir)
#                src='/opt/nicadd/contrib/piot/CSRtrack_10/CsrTrack10_64'
                src='/opt/nicadd/contrib/piot/CSRtrack_10/CsrTrack10_64'
                shutil.copy2(src, '.')
#                src='/opt/nicadd/contrib/piot/CSRtrack_10/CsrTrack10'
                src='/opt/nicadd/contrib/piot/CSRtrack_10/CsrTrack10'
                shutil.copy2(src, '.')
                src='/pdata/prokop/GlueTrack/CurrentExecutables/csrtrack_x64_openmpi_1.4.3'
                shutil.copy2(src, '.')
                src='/pdata/prokop/GlueTrack/CurrentExecutables/CSRtrack_ThisNode'
                shutil.copy2(src, '.')
                src='/pdata/prokop/GlueTrack/CurrentExecutables/CSRtrack_OneNode'
                shutil.copy2(src, '.')
                src='/pdata/prokop/GlueTrack/CurrentExecutables/OneProcessor_phpc1.LINUX'
                shutil.copy2(src, '.')
                os.system('pwd')
                os.system('ls')
#                copyfiles(src,'.')
#                copyfiles(properties['deck_inp'],'/csrtrk.in')
                shutil.copy2(properties['deck_inp'],'csrtrk.in')
#                names = os.listdir(src)
                if not(os.path.exists('out')):
                    os.mkdir('out')
                print(command)
                os.system(command)
#                os.system('rm CsrTrack10_64')
#                os.system('rm CsrTrack10')
                os.system('rm csrtrk.in')
                os.system('rm in/genericGlueTrackFileName.fmt1')
#                for name in names:
#                    os.remove(name)
                os.chdir('..')
        if name=='os':
            print('task: ', name)
            if 'command' in properties: 
                command=properties['command']
                command=command.replace(',',' ')
                print(command)
                os.system(command)
    return indfig




# CRP 11/09/2011
# Parallel CSRtrack is ultimately more useful than serial CSRtrack.  However, the command prompts are quite a bit different, so I will add a copied CSRTrack_Parallel command to handle this, with a separate command line and user specification of the machines.Linux file to be used to specify the node list.
######/usr/mpi/gcc/openmpi-1.4.2/bin/orterun --mca btl_tcp_if_exclude lo,eth0 -np 40 -hostfile machines_phpc0.LINUX /opt/nicadd/contrib/piot/CSRtrack_10/CsrTrack10_64
def CSRTrack_Parallel(tasks,fig_dir,html,indfig):
    for i in range(len(tasks)):
        task=tasks[i];  name=task['name']; properties=task['properties']
        if name=='run':
            print('task: ', name)
            if 'command' in properties: 
                command=properties['command']
                print(command)
                os.system('pwd')
                work_dir='./'+properties['work_dir']
                dist_inp='./'+properties['dist_inp']
                deck_inp='./'+properties['deck_inp']
                node_list='./'+properties['node_list']
                num_nodes=properties['num_nodes']
                dist_std=work_dir+'/in/genericGlueTrackFileName.fmt1'
                command=command.replace(',',' ')
                os.system('pwd')
                os.system ('cp %s %s' % (dist_inp, dist_std))
                os.chdir(work_dir)
                src='/opt/nicadd/contrib/piot/CSRtrack_10/CsrTrack10_64'
                shutil.copy2(src, '.')
                src='/opt/nicadd/contrib/piot/CSRtrack_10/CsrTrack10'
                shutil.copy2(src, '.')
                os.system('pwd')
                os.system('ls')
#                copyfiles(src,'.')
#                copyfiles(properties['deck_inp'],'/csrtrk.in')
                shutil.copy2(properties['deck_inp'],'csrtrk.in')
#                names = os.listdir(src)
                if not(os.path.exists('out')):
                    os.mkdir('out')
# Create Command based on previous material
####/usr/mpi/gcc/openmpi-1.4.2/bin/orterun --mca btl_tcp_if_exclude lo,eth0 -np 40 -hostfile machines_phpc0.LINUX /opt/nicadd/contrib/piot/CSRtrack_10/CsrTrack10_64
#                command='/usr/mpi/gcc/openmpi-1.4.2/bin/orterun --mca btl_tcp_if_exclude lo,eth0 -np ' + num_nodes + ' -hostfile ' + node_list + ' /opt/nicadd/contrib/piot/CSRtrack_10/CsrTrack10_64'
                command='./parallelCSR_node40launcher ' + num_nodes + ' ' + node_list
                print(command)
#                os.system(command)
                os.system(command)
                os.system('rm CsrTrack10_64')
                os.system('rm CsrTrack10')
                os.system('rm csrtrk.in')
                os.system('rm in/genericGlueTrackFileName.fmt1')
#                for name in names:
#                    os.remove(name)
                os.chdir('..')
        if name=='os':
            print('task: ', name)
            if 'command' in properties: 
                command=properties['command']
                command=command.replace(',',' ')
                print(command)
                os.system(command)
    return indfig




# CRP 08/18/11
# Uses a new script "DumpElegantParticleToText", which converts an Elegant output deck (any sdds dumped at a cross, or the .out file itself), and converts it in to a custom format that I call "Custom 1"
def SDDSprocess(tasks,fig_dir,html,indfig):
    for i in range(len(tasks)):
        task=tasks[i];  name=task['name']; properties=task['properties']
        if name=='run':
            print('task: ', name)
            if 'command' in properties: 
                command=properties['command']
                work_dir='./'+properties['work_dir']
                command=command.replace(',',' ')
                os.chdir(work_dir)
#                src='../Codes/CSRtrack'
                src='./MiscScripts'
#                copyfiles(src,'./')
#                names = os.listdir(src)
                if not(os.path.exists('TempDistributions')):
                    os.mkdir('TempDistributions')
                print(command)
                os.system(command)
#                os.system('./MiscScripts/DumpElegantParticleToText InitialDistributions/072611_90deg_dist.xUBC1 > TempDistributions/processedSDDSfile')
#                for name in names:
#                    os.remove(name)
                os.chdir('..')
        if name=='os':
            print('task: ', name)
            if 'command' in properties: 
                command=properties['command']
                command=command.replace(',',' ')
                print(command)
                os.system(command)
    return indfig


def FindRestartParameters(html):
    f=html
    indfig=0;start=0;stop=100000
    for line in f:
        if line.find('done')>=0:    start=start+1;
        if line.find('fig')>=0:     indfig=indfig+1;
    if start!=0: html.write('\n\n<H1>!!! RESTART </H1>\n')
    return  start,stop,indfig
    
#def plotXXS






#def plotXPX


                   
            
            
    
     









