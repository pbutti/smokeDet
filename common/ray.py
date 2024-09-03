import numpy as np
import math

class ray:
    
    #This is defined wrt (0,0)
    def __init__(self,o,a,l,sigma_angle=0.,sigma_x=0.,sigma_y=0.):
        
        #This is the original information
        self.origin = o
        self.angle  = a
        
        #This is the ray reference point at creation (equals the origin)
        self.refPoint = self.origin
        
        #y-y0 = m*(x-x0), with m = tan(angle)
        #k    = -m*x0 +y0
        #Ay = Bx + C
        self.tanA = math.tan(a)
        self.m = math.tan(a)
        self.k = -self.m*self.origin[0] + self.origin[1]        
        self.length = l
        self.sigma_angle = sigma_angle
        self.sigma_x = sigma_x
        self.sigma_y = sigma_y
        
        self.d0 = self.getD0fromPoint(self.origin) #this should be always 0.0
        self.d0err = self.getInitialD0Error() #initial error coming from the uncertainty of the origin+angle

        #2x2 diagonal matrix
        self.covariance = np.array([[self.d0err*self.d0err,0.0],[0.0,self.sigma_angle*self.sigma_angle]])
        
        self.end = [self.origin[0] + math.cos(self.angle)*self.length, self.origin[1] + math.sin(self.angle)*self.length]
        self.plotOrigin = [self.origin[0], self.end[0]]
        self.plotEnd    = [self.origin[1], self.end[1]]
        
    def printTrack(self):
        print("Track Properties:")
        print("d0=",self.d0," phi0=",self.angle)
        print(self.covariance)
        
    #I'm not fully sure this is correct, but should give a conservative error for d0
    #There should be a correlation between the d0 and angle that doesn't appear from this formula. 
    #I need to understand if this is correct.
    
    def getInitialD0Error(self):
        sinphi = math.sin(self.angle)
        cosphi = math.cos(self.angle)
        #d0new  = (self.d0 + dx * sinphi - dy * cosphi)
        #Let's compute d0error by the error propagation.
        dd0dx2   = sinphi*sinphi
        dd0dy2   = (-cosphi)*(-cosphi)
        dd0dphi2 = (self.sigma_x*cosphi + self.sigma_y*sinphi)*(self.sigma_x*cosphi + self.sigma_y*sinphi)
        sx2 = self.sigma_x*self.sigma_x
        sy2 = self.sigma_y*self.sigma_y
        sa2 = self.sigma_angle*self.sigma_angle
        #print(sx2,sy2,sa2)
        #print(dd0dx2*sx2,dd0dy2*sy2,dd0dphi2*sa2)
        
        d0err = math.sqrt(dd0dx2*sx2 + dd0dy2*sy2 + dd0dphi2*sa2)
        return d0err
    
        
    #Signed distance, angle between [0,2pi]
    def transformToPoint(self,point):
        dx = self.refPoint[0] - point[0]
        dy = self.refPoint[1] - point[1]
        angle=self.angle
        
        #while (angle > math.pi / 2):
        #    angle-=math.pi
        #while (angle < -math.pi / 2):
        #    angle+=math.pi
        
        sinphi = math.sin(angle)
        cosphi = math.cos(angle)
        d0new  = (self.d0 + dx * sinphi - dy * cosphi)
        self.d0 = d0new
        
        #Jacobian of the transformation
        J = np.array([[1.,dx*cosphi + dy*sinphi],[0.0,1.]])
        JT = J.transpose()
        
        self.covariance = J @ (self.covariance @ JT)
        self.d0err = math.sqrt(self.covariance[0,0])
        #Change the reference point
        self.refPoint = point
        
    
    #d = |k + mx0 -y0| / sqrt(1+m2)
    #Signed distance without transforming
    def getD0fromPoint(self,point):
        #return (-1)*(self.m*point[0] - point[1] + self.k) / math.sqrt(1 + self.m*self.m)
        return (self.m*point[0] - point[1] + self.k) / math.sqrt(1 + self.m*self.m)




## Helper functions for testing purposes
    
def generateOrigins(point,size):
    x1 = np.random.uniform(point[0],point[0]+size[0])
    y1 = np.random.uniform(point[1],point[1]+size[1])
    return [x1,y1]

def generateAngles(angle,size):
    angle = np.random.uniform(angle,angle+size)
    return angle

#y - y1 = (y2-y1) / (x2-x1) *(x-x1)
#m = (y2-y1) / (x2-x1)
#k = -(y2-y1) / (x2-x1)*(x1) + y1
#This generates some noisy data towards a particular point given an origin 
def generatePointingRay(point, origin,sigmaAngle=0.05,sigmaXY=[4,4],sigmaL=2,distr="uniform"):
    scaleFactor=2
    angle  = math.atan2((point[1] - origin[1]), (point[0] - origin[0]))
    dx = (point[0] - origin[0])
    dy = (point[1] - origin[1])
    
    length = math.sqrt(dx*dx + dy*dy)
    
    if distr=="uniform":
        rand_angle  = generateAngles(angle,sigmaAngle)
        #rand_origin = generateOrigins(origin,sigmaXY)
        rand_length = generateAngles(length,sigmaL)*scaleFactor

    return [origin,rand_angle,rand_length]
