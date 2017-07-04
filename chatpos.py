import astropy.coordinates as coord
from astroquery.irsa import Irsa
import astropy.units as u
from astropy.coordinates import SkyCoord
from pylab import *
Irsa.ROW_LIMIT = 10000


def get_center(ra,dec):
	"""
	   This function returns the optimal position in the CHAT CCD for
	   a given traget based on the comparison stars to make good differential photometry
	"""
	ccdx = 2048
	ccdy = 2064
	d1 = 20
	d2 = 300
	scale = 0.62
	delt = 50
	width = scale*3*(ccdx - 2*d2)/60.
	perc = 2.

	#ra, dec = 154.271213, -25.276251
	#ra, dec = 180.113964, -47.243851
	#ra, dec = 217.877016, -59.801907
	#ra, dec = 288.951172, -56.354244

	infocus = False

	table1 = Irsa.query_region(coord.SkyCoord(ra, dec, unit=(u.deg,u.deg)), catalog="fp_psc", spatial="Cone",radius=10.*u.arcmin/60.)

	tlen = len(table1['j_m'])
	if tlen==1:
		target = table1[0]
	else:
		ras,decs = table1['ra'],table1['dec']
		dists = np.sqrt((ras-ra)**2 + (decs-dec)**2)
		I = np.argmin(dists)
		target = table1[I]

	I = np.where((table1['j_m'] - target['j_m'] < 2)&(table1['j_m'] - target['j_m'] != 0 ))[0]
	if len(I)>0:
		infocus =True

	otable = Irsa.query_region(coord.SkyCoord(ra, dec, unit=(u.deg,u.deg)), catalog="fp_psc", spatial="Cone",radius=22 * u.arcmin)
	I = np.argsort(otable['j_m'])
	#print target['j_m'],lim1,lim2
	cond = True
	ni=0
	while cond:
		lim1 = 2.5*np.log10(0.5)
		lim2 = 2.5*np.log10(1.3)+ni
		print target['j_m']+lim2
		I = np.where((otable['j_m']>target['j_m']+lim1) & (otable['j_m']<target['j_m']+lim2))[0]
		table = otable[I]
		I = np.argsort(table['j_m'])
		table = table[I]
		if len(table)>20 or lim2>10:
			cond = False
		else:
			ni+=0.1
	#print table
	#print table
	w2 = (table['j_h']-target['j_h'])/target['j_h']
	w3 = (table['j_k']-target['j_k'])/target['j_k']

	I = np.where(np.absolute(w2)<0.001)[0]
	if len(I)>0:
		w2[I] = 0.001
	I = np.where(np.absolute(w3)<0.001)[0]
	if len(I)>0:
		w3[I] = 0.001

	dist2 = np.sqrt((w2)**2 + (w3)**2)
	I = np.argsort(dist2)
	subtable =  table[I]
	weights  = dist2[I]

	good_ids=[]
	for i in range(len(I)):
		temp_table= Irsa.query_region(coord.SkyCoord(subtable[i]['ra'], subtable[i]['dec'], unit=(u.deg,u.deg)), catalog="fp_psc", spatial="Box",width=12.*u.arcmin/60.)
		vecc = np.array(temp_table['j_m'] - subtable[i]['j_m'])
		I = np.where((vecc<3)&(vecc!=0))[0]
		if  len(I)==0:
			good_ids.append(i)
	good_ids = np.array(good_ids)
	subtable = subtable[good_ids]	

	rangex1 = np.arange(d2,int(0.5*ccdx-d1),delt)
	rangex2 = np.arange(int(0.5*ccdx+d1),ccdx-d2,delt)
	rangex  = np.hstack((rangex1,rangex2))

	rangey = np.arange(d2,ccdy-d2,delt)
	counter = []
	nstars = []
	refx, refy =[],[]
	for x in rangex:
		for y in rangey:
			counts = 0.
			number = 0
			for i in range(len(subtable['j_m'])):
				comp = subtable[i]
				weight = weights[i]
				pra  = x  - 3600.*(target['ra'] - comp['ra'])/scale
				pdec = y  - 3600.*(target['dec'] - comp['dec'])/scale

				if pra > d2 and pra < ccdx - d2 and pdec > d2 and pdec < ccdy - d2 and (target['ra']!=comp['ra'] and target['dec']!=comp['dec']) :
					counts += 1./weight
					number += 1
			#print pra,pdec,counts, number
			counter.append(counts)
			nstars.append(number)
			refx.append(x)
			refy.append(y)
	counter,nstars = np.array(counter),np.array(nstars)
	refx,refy = np.array(refx),np.array(refy)
	Im = np.argmax(counter)
	#print Im
	I = np.where(counter == np.max(counter))[0]
	refx,refy,counter,nstars = refx[I],refy[I],counter[I],nstars[I]
	if len(I)>0:
		posx = np.median(refx)
		posy = np.median(refy)
	else:
		posx = refx
		posy = refy


	shiftx = (posx - 0.5*ccdx)*scale/3600.
	shifty = (posy - 0.5*ccdy)*scale/3600.
	#print posx, posy
	#print shiftx,shifty
	cenra,cendec = ra - shiftx, dec - shifty

	#print cenra,cendec
	I = []
	for i in range(len(subtable['j_m'])):
		comp = subtable[i]
		pra =  posx - 3600.*(target['ra'] - comp['ra'])/scale
		pdec = posy - 3600.*(target['dec'] - comp['dec'])/scale
		if pra > d2 and pra < ccdx - d2 and pdec > d2 and pdec < ccdy - d2 and weight != 0 and (target['ra']!=comp['ra'] and target['dec']!=comp['dec']):
			c = SkyCoord(ra=float(comp['ra'])*u.degree, dec=float(comp['dec'])*u.degree, frame='icrs')
			I.append(i)
			#print c.ra.hms, c.dec.dms, comp['j_m']
	I = np.array(I)

	tarx = 0.5*ccdx + (ra - cenra) / (scale/3600.) 
	tary = 0.5*ccdy + (dec - cendec) / (scale/3600.)
	#return cenra,cendec,infocus,subtable[I]
	return tarx, tary 
