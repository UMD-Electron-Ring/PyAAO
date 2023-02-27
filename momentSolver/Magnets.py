class QuadProfile:

    def __init__(self, func, sf=1.0, hardedge=True):
        self.func = func
        self.sf = sf
        self.hardedge = hardedge

    def __repr__(self):
        if self.hardedge:
            return "QuadProfile, dbdx[T/m]: " + str(self.GetValue(0))
        return "QuadProfile" 

    def __str__(self):
        return self.__repr()

    def GetValue(self, x):
        return (self.func(x) * self.sf)

    def SetScaleFactor(self, sf):
        self.sf = sf

class SolenoidProfile:

    def __init__(self, func, sf=1.0, hardedge=True):
        self.func = func
        self.sf = sf
        self.hardedge = hardedge        

    def __repr__(self):
        if self.hardedge:
            return "SolenoidProfile, B[T]: " + str(self.GetValue(0))
        return "SolenoidProfile"

    def __str__(self):
        return self.__repr()      

    def GetValue(self, x):
        return ( self.func(x, self.sf) )

    def SetScaleFactor(self, sf):
        self.sf = sf