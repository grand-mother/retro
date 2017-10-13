import sys
import numpy as np
pi = 3.1415927
import pylab as pl
import matplotlib.pyplot as plt
import linecache
from scipy.fftpack import rfft, irfft, rfftfreq
from scipy.interpolate import interp1d


#===========================================================================================================
def effective_zenith(x_ant, y_ant, z_ant, zen, azim, alpha, injectionheight=3000.) :
#===========================================================================================================	
	zen = zen * np.pi / 180
	azim = azim * np.pi / 180
	alpha = alpha * np.pi / 180
	
	print 'effective zenith angle used for calculating the antenna response'
	###get the effectiv zenith
	## find approx. Xmax position
	# shower direction: where it goes to
	v = np.array([np.cos(azim)*np.sin(zen), np.sin(azim)*np.sin(zen), np.cos(zen)]) # or *-1: change the direction
		
	# vector for zenith=90
	b = np.array([np.cos(azim), np.sin(azim) ,0.])
			# projection v onto b
	v_p= np.linalg.norm(v)* np.cos(0.5*pi-zen) *b/np.linalg.norm(b)
	v_p= v_p/np.linalg.norm(v_p)

	hor_dist=8000. # approx. horizontal distance from injection point to xmax
	Inj= np.array([0.,0.,injectionheight])
	Xmax= Inj + np.array([0,0,1.]) * hor_dist/ np.tan(zen) + hor_dist * v_p
	#print(v, Xmax)
	#print ' for x,y, :', v_p*8000.,   'for z:',  8000. / np.tan(zen),  Xmax
		
	# antenna position
	ant= np.array([x_ant, y_ant, z_ant])
	# unit vetor between xmax and antenna
	u_xmax=  Xmax - ant
	u_xmax= u_xmax/np.linalg.norm(u_xmax)
	#print u_xmax
		
	# antenna vector
	# at the moment slope always facing the shower => azim_ant = 180. -azim, mountain slope alpha works as zenith
	u_ant= np.array([ np.cos(pi-azim)* np.sin(alpha),  np.sin(pi-azim)* np.sin(alpha), np.cos(alpha)])
	#print(u_ant)
		
	# get effective zenith, the angle between antenna vector and vector between xmax and antenna
	cos_zen_eff= np.dot(u_xmax, u_ant)
	zen_eff=  np.arccos(cos_zen_eff)
		
	print 'effective zenith: ', zen_eff*180./pi, 'deg'

	zen_eff = zen_eff*180./pi

	return [zen_eff , Xmax]


### efield in V/m,azimuth, zenith and alpha in deg, time in s # # # EW=1: EW file, EW=0: NS file
#===========================================================================================================
def get_voltage(time1=None,Ex=None, Ey=None, Ez=None, zenith=None, azimuth =None, EW=1):
#===========================================================================================================
    effective=1 # if you wanna use the effective zenith angle between antenna and shower set value to 1


    azim=azimuth
    
        
    delt = time1[1]-time1[0];
    Fs = 1/delt
    
    time_off=time1[0] # time offset, to get absolute time
    time1 = (time1-time1[0]) #resetted to zero
    print 'Efield signal length: ', time1[-1]*1e9 , 'ns.'


##### the angle conversion    
    zen = zenith*pi/180;
    azim = azim*pi/180;
	
    ## azimuth: internal conversion from GRAND conventions to antenna response conventions since GRAND North=0deg, antenna EAST=0deg
    azim= azim + 0.5*pi # az_ant=az_GRAND +90deg
    # TODO zenith correction: is NEC running fo upgoing or downgoing showers , at the moment that is like GRAND conv if this is true
    # zenith= 0,..,90deg
    print "azimuth in antenna frame ", azim* 180./pi
    
##### efield in antenna frame
    szen = np.sin(zen);
    czen = np.cos(zen);
    saz = np.sin(azim);
    caz = np.cos(azim);

    #amplitudet = czen*(caz*Ex+saz*Ey)-szen*Ez
    #amplitudep = -saz*Ex+caz*Ey
    amplitudet = czen*(caz*Ey-saz*Ex)-szen*Ez
    amplitudep = -saz*Ey-caz*Ex
    #Anne :
    #Above we have to rotate the vector by 90deg around the z axis to go into the NEC sytsem 
    #(y_NEC = x_Grand & x_NEC=-y_GRAND, since delta az( NEC- GRAND) =90deg => 90deg rotation around zaxis) 
    #and then calculate the efield vector in the antenna system.

    ## for the antenna response do back to degree
    azim = azim* 180./pi
    zen = zen* 180./pi
    
    
    #if (azim > 180.):
        #azim = azim- 360.
    #if (azim < 0):
        #azim = abs(azim) # symmetric
        
    #if (zen > 90.):
        #zen= 180. - zen # symmetric
        #print zen



##################################
    #fileleff='1ant10mavgnd_leff.npy' #not loaded,
    if EW==1:
        fileleff=antenna_path+'butwidebasedoublehflat_4p5mhalf_leff.npy' #last one, not loaded
    if EW==0: #NS component
        fileleff=antenna_path+'butwidebasedoublehflat_4p5mhalf_leff.npy' #last one, not loaded --- has to be changed to NS file
        
        # still our fake NS component
        azim= azim + 90. #0.5*pi # az=az +90deg, to get antenna response for a rotated shower which would fake the NS antenna component
        if azim >360.:
            azim= azim- 360.


    azstep=1 #step in azimuth in npy file
    freqscale=1 #freq*2 if h/2 and sizeant/2
    outputpower=0 #if wanted output is power
    loaded=1 #if antenna is loaded or not in npy file

    # imp RLC R = 300;C = 6.5e-12;L = 2e-6;
    #RL =np.array([2.946189699625225,   2.693640403293285,   2.368014134942750,   2.048063136874684,   1.763350109975618,   1.520403822949526,   1.316647804763135, 1.146675541086378,   1.004790940858134,   0.885905649139099,   0.785760827398088,   0.700893633058194,   0.628523271578350,   0.566425575932399,  0.512820527310258, 0.466279040988616])
    #RL=RL*100
    #XL =np.array([  -0.398164981757440,  -0.908417628420895,  -1.223336201352510,  -1.396290370230734,  -1.476701296666265,  -1.499861221583199,  -1.488751816960041,  -1.457793478766247,  -1.415898297104643,  -1.368535029964790,  -1.319038439288611,  -1.269420739673441,  -1.220871947347262,  -1.174069331311762,  -1.129370040642129, -1.086931910885215])
    #XL=XL*100

    #impRLC R = 300;C = 6.5e-12;L = 1e-6; 20 300 MHz
    RLp=np.array([0.536733768083299,   0.840010121593293,   1.200896057862110,   1.600090229038176,   2.006667705049151,   2.381373444983652,   2.685327039754095,2.890920915220645,   2.989008053352027,   2.988573924755322,   2.910298546933116,   2.778572562448605,   2.615623513015058,   2.438774302675680,2.260069404472798,   2.087095391770503,   1.924138291818527,   1.773244204213567,   1.635038639773870 ,  1.509302421351773,   1.395352352338537,1.292281470704730,   1.199104504595046 ,  1.114841938895556,   1.038565549999183,   0.969420412709374 ,  0.906632958509609 ,  0.849511074453546, 0.797439919119660,   0.749875669027360,   0.706338496130072 ,  0.666405514438127,   0.629704091546769,   0.595905715891539 ,  0.564720490529909, 0.535892256306097,   0.509194310914832,   0.484425672932856 ,  0.461407833493318,   0.439981938133439,   0.420006344534532 ,  0.401354506652400, 0.383913141085314,   0.367580636870446 ,  0.352265674931463 ,  0.337886027975219 ,  0.324367515703664,   0.311643093771139 ,  0.299652058008190, 0.288339348095085 ,  0.277654937150177,   0.267553295648130,   0.257992919745857,   0.248935915510479 ,  0.240347631749694,   0.232196335171802, 0.22445292247744])
    RLp=RLp*100
    XLp=np.array([1.149833973427903 ,  1.347001618559050 ,  1.469876468210024,   1.496656923296413,   1.411838459831799,   1.213746617081856 ,  0.919238711558534, 0.561550538777911 ,  0.181259529550333 , -0.184790883266832 , -0.510938360781749,  -0.784380139061162 , -1.002688474655982 , -1.169915727151229, -1.293172262455534 , -1.380332931202411 , -1.438786540600538,  -1.474902574702375,  -1.493909155794964 , -1.499971154708314,  -1.496345170687206,  -1.485548051254960 , -1.469510769217091 , -1.449707994034062,  -1.427262501557595 , -1.403027192021063 , -1.377648559710691,  -1.351615388245273,  -1.325295941574338 , -1.298966315222550 , -1.272832045980507,  -1.247044599699961 , -1.221713972961578 , -1.196918345352935 , -1.172711490165158, -1.149128477825457 , -1.126190075642858 , -1.103906149182129 , -1.082278296775352,  -1.061301893203182,  -1.040967676805738,  -1.021262982755671,  -1.002172701363368 , -0.983680022166382,  -0.965767010753354 , -0.948415054722766 , -0.931605207084645 , -0.915318449184856 , -0.899535890421292,  -0.884238918293781 , -0.869409309431790 , -0.855029309984293 , -0.841081691988703 , -0.827549790949401 , -0.814417528766047 , -0.801669425292115, -0.789290601124629])
    XLp=XLp*100
    fr=np.arange(20,301,5)

    freq,realimp,reactance,theta,phi,lefftheta,leffphi,phasetheta,phasephi=np.load(fileleff)
    RL=interp1d(fr, RLp, bounds_error=False, fill_value=0.0)(freq[:,0])
    XL=interp1d(fr, XLp, bounds_error=False, fill_value=0.0)(freq[:,0])

    nfreq=len(freq[:,0])
    f=np.zeros(nfreq)
    RA=np.zeros(nfreq)
    XA=np.zeros(nfreq)
    ltr=np.zeros(nfreq)
    lta=np.zeros(nfreq)
    lpr=np.zeros(nfreq)
    lpa=np.zeros(nfreq)
    if azstep==5:
        roundazimuth=round(azim/10)*10+round((azim-10*round(azim/10))/5)*5
    elif azstep==1:
        roundazimuth=round(azim)
    else:
        print('Error on azimuth step!')
        return(0)
    if roundazimuth>=91 and roundazimuth<=180:
        roundazimuth=180-roundazimuth
    if roundazimuth>=181 and roundazimuth<=270:
        roundazimuth=roundazimuth-180
    if roundazimuth>=271 and roundazimuth<=360:
        roundazimuth=360-roundazimuth
    for i in range(nfreq):
        f[i]=freq[i,0]*freqscale
        indtheta=np.nonzero(theta[i,:]==round(zen))[0]
        indphi=np.nonzero(phi[i,:]==roundazimuth)[0]
        indcom=np.intersect1d(indtheta,indphi)
        ltr[i]=lefftheta[i,indcom]
       	lta[i]=phasetheta[i,indcom]*np.pi/180
       	lpr[i]=leffphi[i,indcom]
       	lpa[i]=phasephi[i,indcom]*np.pi/180
        if loaded==0:
            RA[i]=realimp[i,0]
            XA[i]=reactance[i,0]
            Rlefft=ltr[i]*np.cos(lta[i])
            Xlefft=ltr[i]*np.sin(lta[i])
            Rleffp=lpr[i]*np.cos(lpa[i])
            Xleffp=lpr[i]*np.sin(lpa[i])
            Rleqt=((Rlefft*RL[i]-Xlefft*XL[i])*(RA[i]+RL[i]) + (Rlefft*XL[i]+Xlefft*RL[i])*(XA[i]+XL[i])) / ((RA[i]+RL[i])**2+(XA[i]+XL[i])**2)
            Xleqt=((Rlefft*RL[i]+Xlefft*XL[i])*(XA[i]+XL[i]) + (Rlefft*XL[i]+Xlefft*RL[i])*(RA[i]+RL[i])) / ((RA[i]+RL[i])**2+(XA[i]+XL[i])**2)
            ltr[i]=np.sqrt(Rleqt**2+Xleqt**2)
            print(Rleqt,Xleqt,ltr[i])
            lta[i]=np.arccos(Rleqt/ltr[i])
            Rleqp=((Rleffp*RL[i]-Xleffp*XL[i])*(RA[i]+RL[i]) + (Rleffp*XL[i]+Xleffp*RL[i])*(XA[i]+XL[i])) / ((RA[i]+RL[i])**2+(XA[i]+XL[i])**2)
            Xleqp=((Rleffp*RL[i]+Xleffp*XL[i])*(XA[i]+XL[i]) + (Rleffp*XL[i]+Xleffp*RL[i])*(RA[i]+RL[i])) / ((RA[i]+RL[i])**2+(XA[i]+XL[i])**2)
            lpr[i]=np.sqrt(Rleqp**2+Xleqp**2)
            print(Rleqp,lpr[i])
            lpa[i]=np.arccos(Rleqp/lpr[i])
    if loaded==0:#phases are not unwrap! so:        
        for i in range(1,nfreq):
            while lpa[i]-lpa[i-1]<-180*np.pi/180:
                lpa[i]=lpa[i]+360*np.pi/180
            while lpa[i]-lpa[i-1]>180*np.pi/180:
                lpa[i]=lpa[i]-360*np.pi/180
            while lta[i]-lta[i-1]<-180*np.pi/180:
                lta[i]=lta[i]+360*np.pi/180
            while lta[i]-lta[i-1]>180*np.pi/180:
                lta[i]=lta[i]-360*np.pi/180

    #print(round(zenith),roundazimuth,f,ltr,lta,lpr,lpa,Fs)

    fmin=f[0]
    fmax=f[-1]
    f=f*1e6

    nf  = int(2**np.floor(np.log(len(amplitudet))/np.log(2)))
    while Fs/nf > fmin*1e6:   # <== Make sure that the DFT resolution is at least fmin.
        nf *= 2
    F = rfftfreq(nf)*Fs


    modulust = interp1d(f, ltr, bounds_error=False, fill_value=0.0)(F)
    phaset   = interp1d(f, lta, bounds_error=False, fill_value=0.0)(F)
    modulusp = interp1d(f, lpr, bounds_error=False, fill_value=0.0)(F)
    phasep   = interp1d(f, lpa, bounds_error=False, fill_value=0.0)(F)
    if outputpower:
        RLinter  = interp1d(f, RL, bounds_error=False, fill_value=0.0)(F)

    #print(F,modulust)

    phaset -= phaset[0] # Switch the phase origin to be consistent with a real signal.
    phasep -= phasep[0] # Switch the phase origin to be consistent with a real signal.

    #if we want P=V2/RL -> incorrect
    if outputpower:
        modulusp[RLinter!=0]=modulusp[RLinter!=0]/np.sqrt(RLinter[RLinter!=0])
        modulust[RLinter!=0]=modulust[RLinter!=0]/np.sqrt(RLinter[RLinter!=0])

    #B and D are V in freq domain, they are complex
    A = rfft(amplitudet, nf)
    ct = np.cos(phaset)
    st = np.sin(phaset)
    B = np.zeros(A.shape)
    B[1:-1:2] = modulust[1:-1:2]*(A[1:-1:2]*ct[1:-1:2]-A[2:-1:2]*st[2:-1:2])
    B[2:-1:2] = modulust[2:-1:2]*(A[1:-1:2]*st[1:-1:2]+A[2:-1:2]*ct[2:-1:2])
    B[0]  = A[0]*modulust[0]
    B[-1] = A[-1]*modulust[-1]

    C = rfft(amplitudep, nf)
    cp = np.cos(phasep)
    sp = np.sin(phasep)
    D = np.zeros(C.shape)
    D[1:-1:2] = modulusp[1:-1:2]*(C[1:-1:2]*cp[1:-1:2]-C[2:-1:2]*sp[2:-1:2])
    D[2:-1:2] = modulusp[2:-1:2]*(C[1:-1:2]*sp[1:-1:2]+C[2:-1:2]*cp[2:-1:2])
    D[0]  = C[0]*modulusp[0]
    D[-1] = C[-1]*modulusp[-1]



#we should apply 1/sqrt(RL) to the real part
#and then put vt and vp squared, after ifft

    vt=irfft(B)
    vp=irfft(D)

    if outputpower:
        vt=vt**2
        vp=vp**2
        
    voltage = vp+ vt
    

    
    timet     = np.arange(0, len(vt))/Fs
    timep     = np.arange(0, len(vp))/Fs
    
    print '\n'
    print 'peak to peak voltage amplitude = ', max(voltage) - min(voltage)

    return(voltage, timet+time_off)



#===========================================================================================================
# Compute the time dependent voltage
#===========================================================================================================
## example: python computevoltage.py 85 45 0 ./ 0 a0.trace 
## if effective zenith wanted -- inside function, set effective to 1 plus hand over antenna postion x,y,z in m:  python computevoltage.py 85 45 0 ./ 0 a0.trace 0 34000 3000
#'''
# Read zen and az value from input
zenith_sim=float(sys.argv[1])
azimuth_sim=float(sys.argv[2])

## include a mountain slope - correction of zenith angle
alpha_sim=float(sys.argv[3])

# which efield trace do you wanna read in. to be consistent the script works with the antenna ID
path=sys.argv[4] #folder containing the traces and where the output should go to
l=sys.argv[5] # antenna ID
efieldtxt=path+'a'+str(l)+'.trace'

print 'Wave direction in GRAND conv.: zenith = ', zenith_sim, ' deg, azimuth = ', azimuth_sim, 'deg.', 'efield file: ', efieldtxt , ' mountain slope: ', alpha_sim


# Model the input signal.
#time1_sim, Ex_sim, Ey_sim,Ez_sim = np.loadtxt(efieldtxt,delimiter='  ',usecols=(1,2,3,4),unpack=True)
time1_sim, Ey_sim, Ex_sim,Ez_sim = np.loadtxt(efieldtxt,delimiter='  ',usecols=(1,2,3,4),unpack=True)
# NOTE: adapt to your time from whatever to s
time1_sim= time1_sim*1e-6 # time has to be handed in s

print 'Now computing antenna response...'

effective = float(sys.argv[6])

# Compute effeective zenith
if effective==1:
	
	#antenna position
	try :
		
		x_sim = float(sys.argv[7])
		y_sim = float(sys.argv[8])
		z_sim = float(sys.argv[9])
	
	except :
	
		try :
			
			numberline = int(l) + 1
			line = linecache.getline(path + 'antpos.dat', numberline)
			[x_sim, y_sim, z_sim] = map(float, line.split())
		
			print("antenna position at", x_sim, y_sim, z_sim)
	
		except :
		
			print("No antenna position file find, please enter antenna position as arguments like x y z")
			sys.exit()

	[zenith_eff, Xmax] = effective_zenith(x_sim, y_sim, z_sim, zenith_sim, azimuth_sim, alpha_sim)
		
else: # in case effective zenith not wanted, one has to account for mountain slope.
## zenith: correct for mountain slope
	zenith_eff = zenith_sim +alpha # zen_ant= 180deg-(zen_aires-alpha_slope) = zen_GRAND+alpha_slope
	
if zenith_eff > 180 :
	print ' --- GRAND zenith larger than 90deg, antenna response not valid for that angle --- '
	exit()


# Compute the output voltage for EW component
voltage_EW, timeEW  = get_voltage( time1=time1_sim,Ex=Ex_sim, Ey=Ey_sim, Ez=Ez_sim, zenith=zenith_eff, azimuth=azimuth_sim, EW=1)

# Compute the output voltage for NS component -- TODO add here the correct antenna file at some point
#voltage_NS, timeNS = get_voltage(Ex=Ey_sim, Ey=-Ex_sim, Ez=Ez_sim, zenith=zenith_sim, azimuth=azimuth_sim-90, alpha=alpha_sim, time1=time1_sim, EW=0) 
voltage_NS, timeNS = get_voltage( time1=time1_sim,Ex=Ex_sim, Ey=Ey_sim, Ez=Ez_sim, zenith=zenith_eff, azimuth=azimuth_sim, EW=0)


pl.savetxt('out_'+str(l)+'.txt', (timeEW, voltage_EW, voltage_NS))#, voltage_NS))




###plots


plt.figure(1,  facecolor='w', edgecolor='k')
plt.subplot(211)
plt.plot(time1_sim*1e9,Ex_sim, label="Ex = EW")
plt.plot(time1_sim*1e9,Ey_sim, label="Ey = NS")
plt.plot(time1_sim*1e9,Ez_sim, label="Ez = UP")
plt.xlabel('Time [nsec]')
plt.ylabel('Electric field [muV/m]')
plt.legend(loc='best')
##plt.show()
plt.subplot(212)
#plt.ylabel('Etheta [muV]')
#plt.show()
plt.plot(timeEW*1e9,voltage_EW, label="EW")
plt.plot(timeNS*1e9,voltage_NS, label="NS")
plt.xlabel('Time [nsec]')
plt.ylabel('Voltage [muV]')
plt.legend(loc='best')

plt.show()
#'''