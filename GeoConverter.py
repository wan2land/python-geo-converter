"""
porting from Java Source 
by wani (me@wani.kr)

origin source :
	from http://www.androidpub.com/1043970
"""

import math

class GeoPoint :
	def __init__(self, x=0, y=0, z = 0) :
		self.x = x
		self.y = y
		self.z = z

	def getX(self) :
		return self.x

	def getY(self) :
		return self.y

	def getZ(self) :
		return self.z
	def setX(self, x) :
		self.x = x
	def setY(self, y) :
		self.y = y
	def setZ(self, z) :
		self.z = z

def degree2radian( degree ) :
	return degree * math.pi / 180.0

def radian2degree( radian ) :
	return radian * 180.0 / math.pi

def e0fn( x ) :
	return 1.0 - 0.25 * x * (1.0 + x / 16.0 * (3.0 + 1.25 * x))

def e1fn( x ) :
	return 0.375 * x * (1.0 + 0.25 * x * (1.0 + 0.46875 * x))

def e2fn( x ) :
	return 0.05859375 * x * x * (1.0 + 0.75 * x)

def e3fn( x ) :
	return x * x * x * (35.0 / 3072.0)

def mlfn( e0,  e1,  e2,  e3,  phi) :
	return e0 * phi - e1 * math.sin(2.0 * phi) + e2 * math.sin(4.0 * phi) - e3 * math.sin(6.0 * phi)

def asinz( value ) :
	if value > 0 :
		value = min(1, value)
	else :
		value = max(-1, value)
	return math.asin(value)

GEO = 0
KATEC = 1
TM = 2
GRS80 = 3

EPSLN = 0.0000000001

m_arMajor = [6378137.0, 6377397.155, 6377397.155]
m_arMinor = [6356752.3142, 6356078.9633422494, 6356078.9633422494]

m_arScaleFactor = [1, 0.9996, 1.0, 1.0] # KATEC -> 0.9999?
m_arLonCenter = [0.0, 2.22529479629277, 2.21661859489671]
m_arLatCenter = [0.0, 0.663225115757845, 0.663225115757845]
m_arFalseNorthing = [0.0, 600000.0, 500000.0]
m_arFalseEasting = [0.0, 400000.0, 200000.0]

datum_params = [-146.43, 507.89, 681.46]

m_Ind = [0,0,0]
m_Es = [0,0,0]
m_Esp = [0,0,0]
src_m = [0,0,0]
dst_m = [0,0,0]

tmp = m_arMinor[GEO] / m_arMajor[GEO]
m_Es[GEO] = 1.0 - tmp * tmp
m_Esp[GEO] = m_Es[GEO] / (1.0 - m_Es[GEO])

if (m_Es[GEO] < 0.00001) :
	m_Ind[GEO] = 1.0
else :
	m_Ind[GEO] = 0.0


tmp = m_arMinor[KATEC] / m_arMajor[KATEC]
m_Es[KATEC] = 1.0 - tmp * tmp
m_Esp[KATEC] = m_Es[KATEC] / (1.0 - m_Es[KATEC])

if (m_Es[KATEC] < 0.00001) :
	m_Ind[KATEC] = 1.0
else :
	m_Ind[KATEC] = 0.0


tmp = m_arMinor[TM] / m_arMajor[TM]
m_Es[TM] = 1.0 - tmp * tmp
m_Esp[TM] = m_Es[TM] / (1.0 - m_Es[TM])

if (m_Es[TM] < 0.00001) :
	m_Ind[TM] = 1.0
else :
	m_Ind[TM] = 0.0


src_m[GEO] = m_arMajor[GEO] * mlfn(e0fn(m_Es[GEO]), e1fn(m_Es[GEO]), e2fn(m_Es[GEO]), e3fn(m_Es[GEO]), m_arLatCenter[GEO])
dst_m[GEO] = m_arMajor[GEO] * mlfn(e0fn(m_Es[GEO]), e1fn(m_Es[GEO]), e2fn(m_Es[GEO]), e3fn(m_Es[GEO]), m_arLatCenter[GEO])
src_m[KATEC] = m_arMajor[KATEC] * mlfn(e0fn(m_Es[KATEC]), e1fn(m_Es[KATEC]), e2fn(m_Es[KATEC]), e3fn(m_Es[KATEC]), m_arLatCenter[KATEC])
dst_m[KATEC] = m_arMajor[KATEC] * mlfn(e0fn(m_Es[KATEC]), e1fn(m_Es[KATEC]), e2fn(m_Es[KATEC]), e3fn(m_Es[KATEC]), m_arLatCenter[KATEC])
src_m[TM] = m_arMajor[TM] * mlfn(e0fn(m_Es[TM]), e1fn(m_Es[TM]), e2fn(m_Es[TM]), e3fn(m_Es[TM]), m_arLatCenter[TM])
dst_m[TM] = m_arMajor[TM] * mlfn(e0fn(m_Es[TM]), e1fn(m_Es[TM]), e2fn(m_Es[TM]), e3fn(m_Es[TM]), m_arLatCenter[TM])

def convert( srctype, dsttype, in_pt) :
	tmpPt = GeoPoint()
	out_pt = GeoPoint()

	if srctype == GEO :
		tmpPt.setX( degree2radian(in_pt.getX()) )
		tmpPt.setY( degree2radian(in_pt.getY()) )

	else :
		tm2geo(srctype, in_pt, tmpPt)

	if (dsttype == GEO) :
		out_pt.setX( radian2degree(tmpPt.getX() ) )
		out_pt.setY( radian2degree(tmpPt.getY() ) )

	else :
		geo2tm(dsttype, tmpPt, out_pt)
	
	return out_pt

def geo2tm( dsttype, in_pt, out_pt ) :
	
	transform(GEO, dsttype, in_pt)
	delta_lon = in_pt.getX() - m_arLonCenter[dsttype]
	sin_phi = math.sin(in_pt.getY())
	cos_phi = math.cos(in_pt.getY())

	if (m_Ind[dsttype] != 0) :
		b = cos_phi * math.sin(delta_lon)

		if ((math.abs(math.abs(b) - 1.0)) < EPSLN) :
			pass
			#Log.d("infinite error")
			#System.out.println("infinite error")

	else :
		b = 0
		x = 0.5 * m_arMajor[dsttype] * m_arScaleFactor[dsttype] * math.log((1.0 + b) / (1.0 - b))
		con = math.acos(cos_phi * math.cos(delta_lon) / math.sqrt(1.0 - b * b))

		if (in_pt.getY() < 0) :
			con = con * -1
			y = m_arMajor[dsttype] * m_arScaleFactor[dsttype] * (con - m_arLatCenter[dsttype])
	
	al = cos_phi * delta_lon
	als = al * al
	c = m_Esp[dsttype] * cos_phi * cos_phi
	tq = math.tan(in_pt.getY())
	t = tq * tq
	con = 1.0 - m_Es[dsttype] * sin_phi * sin_phi
	n = m_arMajor[dsttype] / math.sqrt(con)
	ml = m_arMajor[dsttype] * mlfn(e0fn(m_Es[dsttype]), e1fn(m_Es[dsttype]), e2fn(m_Es[dsttype]), e3fn(m_Es[dsttype]), in_pt.getY())

	out_pt.setX( m_arScaleFactor[dsttype] * n * al * (1.0 + als / 6.0 * (1.0 - t + c + als / 20.0 * (5.0 - 18.0 * t + t * t + 72.0 * c - 58.0 * m_Esp[dsttype]))) + m_arFalseEasting[dsttype] )
	out_pt.setY( m_arScaleFactor[dsttype] * (ml - dst_m[dsttype] + n * tq * (als * (0.5 + als / 24.0 * (5.0 - t + 9.0 * c + 4.0 * c * c + als / 30.0 * (61.0 - 58.0 * t + t * t + 600.0 * c - 330.0 * m_Esp[dsttype]))))) + m_arFalseNorthing[dsttype] )


def tm2geo( srctype,  in_pt, out_pt) :
	tmpPt = GeoPoint(in_pt.getX(), in_pt.getY())
	max_iter = 6
	if (m_Ind[srctype] != 0) :
		f = math.exp(in_pt.getX() / (m_arMajor[srctype] * m_arScaleFactor[srctype]))
		g = 0.5 * (f - 1.0 / f)
		temp = m_arLatCenter[srctype] + tmpPt.getY() / (m_arMajor[srctype] * m_arScaleFactor[srctype])
		h = math.cos(temp)
		con = math.sqrt((1.0 - h * h) / (1.0 + g * g))
		out_pt.setY( asinz(con) )

		if (temp < 0) :
			out_pt.setY( out_pt.getY * -1 )

		if ((g == 0) and (h == 0)) :
			out_pt.setX( m_arLonCenter[srctype] )

		else :
			out_pt.setX( math.atan(g / h) + m_arLonCenter[srctype] )
		
	tmpPt.setX( tmpPt.getX() - m_arFalseEasting[srctype] )
	tmpPt.setY( tmpPt.getY() - m_arFalseNorthing[srctype] )
	
	con = (src_m[srctype] + tmpPt.getY() / m_arScaleFactor[srctype]) / m_arMajor[srctype]
	phi = con

	i = 0

	while (1) :
		delta_Phi = ((con + e1fn(m_Es[srctype]) * math.sin(2.0 * phi) - e2fn(m_Es[srctype]) * math.sin(4.0 * phi) + e3fn(m_Es[srctype]) * math.sin(6.0 * phi)) / e0fn(m_Es[srctype])) - phi

		phi = phi + delta_Phi

		if (abs(delta_Phi) <= EPSLN) :
			break

		if (i >= max_iter) :
			break
		
		i += 1
	

	if (abs(phi) < (math.pi / 2)) :
		sin_phi = math.sin(phi)
		cos_phi = math.cos(phi)
		tan_phi = math.tan(phi)
		c = m_Esp[srctype] * cos_phi * cos_phi
		cs = c * c
		t = tan_phi * tan_phi
		ts = t * t
		cont = 1.0 - m_Es[srctype] * sin_phi * sin_phi
		n = m_arMajor[srctype] / math.sqrt(cont)
		r = n * (1.0 - m_Es[srctype]) / cont
		d = tmpPt.getX() / (n * m_arScaleFactor[srctype])
		ds = d * d
		out_pt.setY( phi - (n * tan_phi * ds / r) * (0.5 - ds / 24.0 * (5.0 + 3.0 * t + 10.0 * c - 4.0 * cs - 9.0 * m_Esp[srctype] - ds / 30.0 * (61.0 + 90.0 * t + 298.0 * c + 45.0 * ts - 252.0 * m_Esp[srctype] - 3.0 * cs))) )
		out_pt.setX( m_arLonCenter[srctype] + (d * (1.0 - ds / 6.0 * (1.0 + 2.0 * t + c - ds / 20.0 * (5.0 - 2.0 * c + 28.0 * t - 3.0 * cs + 8.0 * m_Esp[srctype] + 24.0 * ts))) / cos_phi) )
	else :
		out_pt.setY( math.pi * 0.5 * math.sin(tmpPt.getY()) )
		out_pt.setX( m_arLonCenter[srctype] )
	
	transform(srctype, GEO, out_pt)


def getDistancebyGeo( pt1,  pt2) :
	lat1 = D2R(pt1.getY())
	lon1 = D2R(pt1.getX())
	lat2 = D2R(pt2.getY())
	lon2 = D2R(pt2.getX())

	longitude = lon2 - lon1
	latitude = lat2 - lat1

	a = math.pow(math.sin(latitude / 2.0), 2) + math.cos(lat1) * math.cos(lat2) * math.pow(math.sin(longitude / 2.0), 2)
	return 6376.5 * 2.0 * math.atan2(math.sqrt(a), math.sqrt(1.0 - a))

def getDistancebyKatec( pt1, pt2 ) :
	pt1 = convert(KATEC, GEO, pt1)
	pt2 = convert(KATEC, GEO, pt2)

	return getDistancebyGeo(pt1, pt2)


def getDistancebyTm( pt1,  pt2) :
	pt1 = convert(TM, GEO, pt1)
	pt2 = convert(TM, GEO, pt2)

	return getDistancebyGeo(pt1, pt2)


def getTimebySec(distance) :
	return math.round(3600 * distance / 4)


def getTimebyMin(distance) :
	return math.ceil(getTimebySec(distance) / 60)


"""
Author:       Richard Greenwood rich@greenwoodmap.com
License:      LGPL as per: http://www.gnu.org/copyleft/lesser.html
"""

"""
 * convert between geodetic coordinates (longitude, latitude, height)
 * and gecentric coordinates (X, Y, Z)
 * ported from Proj 4.9.9 geocent.c
"""

HALF_PI = 0.5 * math.pi
COS_67P5  = 0.38268343236508977  #/* cosine of 67.5 degrees */
AD_C      = 1.0026000 

def transform( srctype, dsttype, point) :
	if (srctype == dsttype) :
		return
	
	if (srctype != 0 or dsttype != 0) :
		# Convert to geocentric coordinates.
		geodetic_to_geocentric(srctype, point)
		
		# Convert between datums
		if (srctype != 0) :
			geocentric_to_wgs84(point)
		
		
		if (dsttype != 0) :
			geocentric_from_wgs84(point)
		
		
		# Convert back to geodetic coordinates
		geocentric_to_geodetic(dsttype, point)

def geodetic_to_geocentric ( type, p ) :

	Longitude = p.getX()
	Latitude = p.getY()
	Height = p.getZ()

	if (Latitude < -HALF_PI and Latitude > -1.001 * HALF_PI ) :
		Latitude = -HALF_PI 

	elif (Latitude > HALF_PI and Latitude < 1.001 * HALF_PI ) :

		Latitude = HALF_PI
	elif ((Latitude < -HALF_PI) or (Latitude > HALF_PI)) : # Latitude out of range
		return true


	#/* no errors */
	if (Longitude > math.pi) :
		Longitude -= (2*math.PI)

	Sin_Lat = math.sin(Latitude)
	Cos_Lat = math.cos(Latitude)
	Sin2_Lat = Sin_Lat * Sin_Lat
	Rn = m_arMajor[type] / (math.sqrt(1.0e0 - m_Es[type] * Sin2_Lat))
	X = (Rn + Height) * Cos_Lat * math.cos(Longitude)
	Y = (Rn + Height) * Cos_Lat * math.sin(Longitude)
	Z = ((Rn * (1 - m_Es[type])) + Height) * Sin_Lat

	p.setX(X)
	p.setY(Y)
	p.setZ(Z)

	return False

def geocentric_to_geodetic ( type, p ) :

	X = p.getX()
	Y = p.getY()
	Z = p.getZ()
	Latitude = 0.

	At_Pole = False
	if (X != 0.0) :
		Longitude = math.atan2(Y,X)

	else :
		if (Y > 0) :
			Longitude = HALF_PI
		
		elif (Y < 0) :
			Longitude = -HALF_PI
		
		else :
			At_Pole = true
			Longitude = 0.0
			if (Z > 0.0) : #/* north pole */
				Latitude = HALF_PI
			
			elif (Z < 0.0) :  #/* south pole */
				Latitude = -HALF_PI
			
			else :  #/* center of earth */
				Latitude = HALF_PI
				Height = -m_arMinor[type]
				return
			
		
	W2 = X*X + Y*Y
	W = math.sqrt(W2)
	T0 = Z * AD_C
	S0 = math.sqrt(T0 * T0 + W2)
	Sin_B0 = T0 / S0
	Cos_B0 = W / S0
	Sin3_B0 = Sin_B0 * Sin_B0 * Sin_B0
	T1 = Z + m_arMinor[type] * m_Esp[type] * Sin3_B0
	Sum = W - m_arMajor[type] * m_Es[type] * Cos_B0 * Cos_B0 * Cos_B0
	S1 = math.sqrt(T1*T1 + Sum * Sum)
	Sin_p1 = T1 / S1
	Cos_p1 = Sum / S1
	Rn = m_arMajor[type] / math.sqrt(1.0 - m_Es[type] * Sin_p1 * Sin_p1)
	if (Cos_p1 >= COS_67P5) :
		Height = W / Cos_p1 - Rn
	
	elif (Cos_p1 <= -COS_67P5) :
		Height = W / -Cos_p1 - Rn
	
	else :
		Height = Z / Sin_p1 + Rn * (m_Es[type] - 1.0)
	
	if (At_Pole == False) :
		Latitude = math.atan(Sin_p1 / Cos_p1)
	
	p.setX(Longitude)
	p.setY(Latitude)
	p.setZ(Height)
	return

def geocentric_to_wgs84( p ) :
	p.setX( p.getX() + datum_params[0] )
	p.setY( p.getY() + datum_params[1] )
	p.setZ( p.getZ() + datum_params[2] )

def geocentric_from_wgs84( p ) :
	p.setX( p.getX() - datum_params[0] )
	p.setY( p.getY() - datum_params[1] )
	p.setZ( p.getZ() - datum_params[2] )


#GeoPoint output
if __name__ == "__main__" :
	print "example :D"

	pt = GeoPoint(205989.36192, 449778.885301)
	print pt.getX(), pt.getY()

	output = convert(TM, GEO, pt)
	print output.getX(), output.getY()

	output = convert(GEO, KATEC, output)
	print output.getX(), output.getY()

