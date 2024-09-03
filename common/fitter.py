import numpy as np
import numpy.linalg as la
import math

#Compute the vertex location by angular Chi2
#Chi2 is of the form Sum_i|alpha_i - atan2(Vy-Oy, Vx-Ox)|^2

class Angular2DFitter():
    def __init__(self, tracks,seed=np.array([0,0])):
        self.tracks = tracks
        self.V = seed
        self.chi2 = 0
        
    def ConvertToZero2Pi(self,angle):
        #This function convert an angle between -pi,pi to 0,2pi range
        print("Original Angle:",angle)
        if (angle < 0):
            angle = 2*math.pi + angle
        print("Return Angle:", angle)
        return angle
        
    #delta_alpha = alpha_i - atan2(Vo_Y - O_Y, Vo_X - O_X)
    def computeDeltaAlpha(self,track,debug=False):
        #t_angle = self.ConvertToZero2Pi(track.angle)
        #v_angle = self.ConvertToZero2Pi(math.atan2(self.V[1] - track.origin[1],self.V[0] - track.origin[0]))
        
        t_angle = track.angle
        v_angle = math.atan2(self.V[1] - track.origin[1],self.V[0] - track.origin[0])
        
        if (debug):
            print("Track Angle=",t_angle)
            print("Angle Error=",track.sigma_angle)
            print("Vertex Angle = ", v_angle)
        
        delta_alpha= t_angle - v_angle
    
        if (debug):
            print("Delta Angle Raw=",delta_alpha )
                        
        if (delta_alpha > math.pi):
            delta_alpha -= 2*math.pi
        if (delta_alpha < -math.pi):
            delta_alpha += 2*math.pi
            
        if (debug):
            print("Delta Angle Corr=",delta_alpha)
            print("--------")
        
        return delta_alpha
                 
        
        
    def computeD(self, track,debug=False):
        Ox = track.origin[0]
        Oy = track.origin[1]
        Vx = self.V[0]
        Vy = self.V[1]
        if (debug):
            print("Ox=",Ox,"Oy=",Oy)
            print("Vx=",Vx,"Vy=",Vy)
                 
        rho2 = (Vx - Ox)**2 + (Vy - Oy)**2
        if (debug):
            print("rho2", rho2)
        dfdVx = -(Vy - Oy) / rho2
        dfdVy = (Vx - Ox) / rho2
        D = np.array([dfdVx, dfdVy])
        if (debug):
            print(D)
        return D
    
    def computeTi(self, track):
        Dtrk = self.computeD(track,False)
        Wi = 1./ (track.sigma_angle**2)
        Ti = Dtrk*Wi*self.computeDeltaAlpha(track)
        return Ti
    
    def computewi(self,track):
        Dtrk = self.computeD(track)
        #In this case Wi is the inverse of sigma_angle2
        Wi = 1./ (track.sigma_angle**2)
        wi = np.outer(Dtrk,Wi*Dtrk)
        return wi
    
    #This returns an tuple made of (vtxLocation, vtxErrorMatrix)
    def computeVtxLocation(self):
        Swi  = np.array([[0.0,0.0],[0.0,0.0]])
        T    = np.array([0.0,0.0])
        for trk in self.tracks:
            T += self.computeTi(trk)
            Swi += self.computewi(trk)
            
        #print("Swi Matrix:", Swi)
        A=la.pinv(Swi)
        V = A @ T
        return (V,A)
