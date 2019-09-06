import numpy as np
from numpy import pi
import matplotlib.pyplot as plt
from matplotlib import rcParams

class GCSlabModel(object):
    
    def __init__(self):
        self._h = []
        self._n = []
        self._a_src = []
        self._Nlayers = 0
        self._alphainv = 0.5 * np.array([[1.0, -1.0], [1.0, 1.0]], dtype=np.complex128)
        
    def add_slab(self, h, n, source_amplitude=None, y_src=None):
        if(source_amplitude is not None):
            if(y_src is None):
                if(np.isinf(h)):
                    y_src = 0
                else:
                    y_src = h/2
                
            self._h.append(y_src)
            self._n.append(n)
            self._a_src.append(source_amplitude)
            self._h.append(h-y_src)
            self._Nlayers += 1 
        else:
            self._h.append(h)
            
        self._n.append(n)
        self._a_src.append(0)
        self._Nlayers += 1 
        
    def directionality(self, wavelength, angle):
        """
        Parameters
        ----------
        wavelength : float
            The vacuum wavelength of the excitation.
        angle : float
            The excitation angle as measured in the top cladding.
        """
        S = np.zeros((2,1), dtype=np.complex128)
        B = np.eye(2,2, dtype=np.complex128)
        
        for i in range(0, self._Nlayers-1):
            m = self._Nlayers - i - 2
            n = self._Nlayers - i - 1
            Bmn = self.__Bmn(m, n, wavelength, angle)
            
            B = Bmn @ B
            S = Bmn @ S + self.__Sm(m, n, wavelength, angle)
            
        k0 = 2*pi / wavelength
        k1 = k0 * self._n[0]
        kN = k0 * self._n[-1]
        
        
        EN = -S[1] / B[1,1]
        E1 = S[0] + B[0,1]*EN
        P1 = k1 * np.abs(E1)**2
        PN = kN * np.abs(EN)**2
        eta =  P1 / (P1 + PN)
        
        return eta[0]
    
    def __Bmn(self, m, n, wavelength, angle):
        Bmn = np.zeros((2,2), dtype=np.complex128)
        
        n1 = self._n[0]
        nm = self._n[m]
        nn = self._n[n]
        angle_m = np.arcsin(n1/nm * np.sin(angle))
        angle_n = np.arcsin(n1/nn * np.sin(angle))

        k0 = 2*pi / wavelength
        kmz = k0 * nm * np.cos(angle_m)
        knz = k0 * nn * np.cos(angle_n)
        
        hn = self._h[n]
        if(np.isinf(hn)): hn = 0

        pnm = knz/kmz
        Bmn[0,0] = (1 + pnm) * np.exp(1j*knz*hn)
        Bmn[0,1] = (1 - pnm) * np.exp(-1j*knz*hn)
        Bmn[1,0] = (1 - pnm) * np.exp(1j*knz*hn)
        Bmn[1,1] = (1 + pnm) * np.exp(-1j*knz*hn)
        
        return 0.5*Bmn
    
    def __Sm(self, m, n, wavelength, angle):
        S_amp = self._a_src[m]
        Sm = S_amp*np.array([[0],[1.0]], dtype=np.complex128)
        Sm = self._alphainv @ Sm
        
        return Sm
        
    def visualize_stack(self):
        """Plot the layer stack

        Generate a plot of the layer stack in which darkness of color corresponds to refractive index
        (darker is higher).  The source region in the stack is overlayed with a hatched pattern.
        """

        try:
                import matplotlib.pyplot as plt
        except Exception as e:
                print(e)
                print('Matplotlib must be installed in order to visualize the grating stack.')
                return

        if(self._Nlayers < 2):
                print('Error: Need at least two layers to visualize')
                return
        elif(self._Nlayers == 2):
                h_avg = 1.0
                x = [0, 2.0]
        else:
                h_avg = np.mean(self._h[1:-1])
                x = [0, (h_avg + np.sum(self._h[1:-1]))]

        nmax = np.max(self._n)*2.0

        f = plt.figure()
        ax = f.add_subplot(111, facecolor='w')

        # Plot all but the first and last layer
        y0 = 0.0
        for i in range(1,self._Nlayers-1):
                y1 = [y0-self._h[i], y0-self._h[i]]
                y2 = [y0, y0]
                if(self._a_src[i] != 0):
                        ax.fill_between(x, y1, y2, facecolor='blue',
                                        alpha=self._n[i]/nmax,
                                        hatch='///', edgecolor='0.5')
                        
                        ax.plot(x, [y1, y1],
                                'r--', linewidth=1.0)
                elif(self._a_src[i-1] != 0):
                        ax.fill_between(x, y1, y2, facecolor='blue',
                                        alpha=self._n[i]/nmax,
                                        hatch='///', edgecolor='0.5')
                else:
                        ax.fill_between(x, y1, y2, facecolor='blue',
                                        alpha=self._n[i]/nmax)
                y0 -= self._h[i]

        # Plot the first and last layer
        y1 = [h_avg, h_avg]
        y2 = [0, 0]
        ax.fill_between(x, y1, y2, facecolor='blue', alpha=self._n[0]/nmax)

        y1 = [y0-h_avg, y0-h_avg]
        y2 = [y0, y0]
        ax.fill_between(x, y1, y2, facecolor='blue', alpha=self._n[-1]/nmax)
        
        plt.gca().set_aspect('equal', adjustable='box')
        ax.set_xlim(x)
        ax.set_ylim([y0-h_avg, h_avg])

        plt.tight_layout()
        plt.show()

