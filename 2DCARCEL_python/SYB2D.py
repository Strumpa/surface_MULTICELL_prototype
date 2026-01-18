## SYBB2D.py 
# translation of sybt2d.m from A. Hébert's book MATLAB to Python
# translation of tij_2d.m from A. Hébert's book MATLAB to Python
# % integration of the collision, escape and transmission probabilities
# % in unstructured finite 2D geometry.


import numpy as np

def akin(n, x):
    """
    Bickley-Naylor function of order n >= 1 at x >= 0.
    
    Faithful translation from MATLAB akin.m by A. Hébert.
    
    Parameters:
    -----------
    n : int
        Order of the function (n >= 1)
    x : float
        Argument (x >= 0)
    
    Returns:
    --------
    float : Bickley-Naylor function value Ki_n(x)
    
    (c) 2005 Alain Hebert, Ecole Polytechnique de Montreal
    Translated to Python 2026
    """
    if n < 1:
        raise ValueError('n must be > 0')
    elif x < -1.0e-10:
        raise ValueError(f'x must be >= 0. x={x:.5g}')
    elif x <= 0:
        bick0 = [1.57079632679489, 1.0, 0.78539816339745, 0.66666666666667, 0.58904862254808, 0.53333333333333,
                 0.4908738521234, 0.45714285714286, 0.42951462060798, 0.4063492063492]
        return bick0[n - 1] if n <= 10 else bick0[9]
    elif x < 1:
        xsq = x * x
        lnx = np.log(x)
        xi1 = x - 0.91169076342161
        xi2 = x - 0.2451000192866
        xi3 = x - 1.0
        xi4 = x - 0.68448452306295
        if n == 1:
            f = ((((((((0.7945473334662959e-4 * x + 0.51674609093834e-4) * x + 0.100122044600498e-1) * x + 
                    0.8766262885587739e-2) * x + 0.5720051721276178) * x + 0.5454840912170553) * x + 
                  0.139868293763185e2) * x + 0.1531257133183402e2) * x + 0.1501664584981108e3) * xi3 - 0.3916967515498982e2
            f = f / ((x + 0.1219596245756669e1) * xi3 - 0.1193155300614385e3)
            f = f + (x * lnx * (((0.2337898258663651e-2 * xsq + 0.4646979722852471) * xsq + 0.3695696751512241e2) * xsq +
                    0.123465484355545e4) * xsq + 0.175237360009281e5) / ((xsq - 0.2256564898552151e3) * xsq + 0.175237360009281e5)
        elif n == 2:
            f = (((0.1403059e-10 * xi1 + 0.4811961706232723e3) + (((0.9084569646859357 * x + 0.490556407762756e2) * x +
                  0.7521727532834893e2) * x + 0.1693807200769639e4) + (((0.1086764750096697e-1 * xi2 - 1.0e-15) * xi2 +
                  0.1190096804348251e1) * (x**4))) * (xi1 * xi1))
            f = f / (0.1810502008060146e4 + (((xi3 - 0.1893578319929816e2) * x - 0.7855291318496802e2) * xi3))
            f = f + (xsq * lnx * (-0.1520774316867189e9) / (((((xsq - 0.5631701819761997e2) * xsq -
                    0.3143123637802091e3) * xsq + 0.211218654889524e6) * xsq - 0.1267311930720872e8) * xsq + 0.3041548633734379e9))
        elif n == 3:
            f = ((((0.1173330390873767e5 * x + 0.1095667013274141e7) * xi4 - 0.379051e-8) * xi4 + 0.1269499275481224e7) *
                 xi3 - 0.7043581454636306e6)
            f = f / (((((((xi3 - 0.5812262590904993e1) * x - 0.17551758398419e2) * x - 0.133516191424771e3) * x -
                       0.3942586515380026e4) * xi3 + 0.6077621585261822e5) * x + 0.2053835980116203e6) * xi3 - 0.2961415636470914e7)
            f = f + (x**3 * lnx * (((0.2631126488553487e-2 * xsq + 0.5562992588150486) * xsq + 0.3721363059831219e2) * xsq +
                    0.3027965765686327e4) / ((xsq - 0.2309130812632629e3) * xsq + 0.1816779459411791e5))
        else:
            ak1 = akin(n - 1, x)
            ak2 = akin(n - 2, x)
            ak3 = akin(n - 3, x)
            f = (ak3 - ak1) * x / (n - 1) + (n - 2) * ak2 / (n - 1)
    elif x < 6:
        sqrtx = np.sqrt(x)
        expx = np.exp(-x)
        xrec = 1.0 / x
        if n == 1:
            f = ((((((0.1822929159877549e2 * xrec + 0.3272001530672078e3) * xrec + 0.1326511766009986e4) * xrec +
                   0.1868734192859498e4) * xrec + 0.1059016416894119e4) * xrec + 0.2427580524508585e3) * xrec + 0.1833164538368226e2) * expx
            f = f / ((((((((xrec + 0.6590511376539962e2) * xrec + 0.5952592332227032e3) * xrec + 0.1687760486772990e4) * xrec +
                         0.1922624187690926e4) * xrec + 0.957005687628236e3) * xrec + 0.2028344160151355e3) * xrec + 0.1462653804563246e2) * sqrtx)
        elif n == 2:
            f = (((((((0.8407469297501269e-1 * xrec + 0.5596498537189973e1) * xrec + 0.4801733781249936e2) * xrec +
                    0.123974074193467e3) * xrec + 0.12440411402683e3) * xrec + 0.5320941946830476e2) * xrec + 0.9534267279889207e1) * xrec + 0.5766817227841408) * expx
            f = f / ((((((((xrec + 0.2014290370371339e2) * xrec + 0.9747773947009136e2) * xrec + 0.1754089481769652e3) * xrec +
                         0.1383006201574071e3) * xrec + 0.5035544525458363e2) * xrec + 0.8124881079392082e1) * xrec + 0.4601255143693006) * sqrtx)
        elif n == 3:
            f = (((((((0.3093393327788074e-2 * xrec + 0.8959830746710818e-1) * xrec + 0.4333691656848653) * xrec +
                    0.7366450143916231) * xrec + 0.5521543746274372) * xrec + 0.1835695652656039) * xrec + 0.2439716682658748e-1) * expx) / sqrtx
        else:
            ak1 = akin(n - 1, x)
            ak2 = akin(n - 2, x)
            ak3 = akin(n - 3, x)
            f = (ak3 - ak1) * x / (n - 1) + (n - 2) * ak2 / (n - 1)
    elif x < 673.5:
        sqrtx = np.sqrt(x)
        expx = np.exp(-x)
        xrec = 1.0 / x
        if n == 10:
            f = ((((((((((xrec + 0.2380952380952381e1) * xrec + 0.1904761904761905e1) * xrec + 0.9523809523809524) * xrec +
                      0.3174603174603175) * xrec + 0.7936507936507937e-1) * xrec + 0.1587301587301587e-1) * xrec +
                   0.2551020408163265e-2) * xrec + 0.3255792150607811e-3) * xrec + 0.3125e-4) / sqrtx) * expx
        elif n == 9:
            f = ((((((((((xrec + 0.2222222222222222e1) * xrec + 0.1666666666666667e1) * xrec + 0.7777777777777778) * xrec +
                      0.2380952380952381) * xrec + 0.5291005291005291e-1) * xrec + 0.9259259259259259e-2) * xrec +
                   0.1286008230452675e-2) * xrec + 0.1422475719179229e-3) * xrec + 0.1190476190476190e-4) / sqrtx) * expx
        elif n == 8:
            f = ((((((((((xrec + 0.2083333333333333e1) * xrec + 0.1458333333333333e1) * xrec + 0.625) * xrec +
                      0.1736111111111111) * xrec + 0.3472222222222222e-1) * xrec + 0.5401234567901235e-2) * xrec +
                   0.6614711266617621e-3) * xrec + 0.6376608658979881e-4) * xrec + 0.4629629629629630e-5) / sqrtx) * expx
        else:
            ak1 = akin(n - 1, x)
            ak2 = akin(n - 2, x)
            ak3 = akin(n - 3, x)
            f = ((n + 2) * ak3 - (n + 1) * ak1) / x + ak2
    else:
        f = 0.0
    return f

def di_f(sig, segment):
    """
    Collision probability function.
    
    Parameters:
    -----------
    sig : float
        Cross section (cm^-1)
    segment : float
        Segment length (cm)
    
    Returns:
    --------
    float : Collision probability
    
    (c) 2008 Alain Hebert, Ecole Polytechnique de Montreal
    Translated to Python 2026
    """
    if sig != 0:
        return segment / sig - (akin(3, 0) - akin(3, sig * segment)) / sig**2
    else:
        return np.pi * segment**2 / 4

def ei_f(tau0, sig, segment):
    """
    Escape probability function.
    
    Parameters:
    -----------
    tau0 : float
        Initial optical depth
    sig : float
        Cross section (cm^-1)
    segment : float
        Segment length (cm)
    
    Returns:
    --------
    float : Escape probability
    
    (c) 2008 Alain Hebert, Ecole Polytechnique de Montreal
    Translated to Python 2026
    """
    if sig != 0:
        return (akin(3, tau0) - akin(3, tau0 + sig * segment)) / sig
    else:
        return segment * akin(2, tau0)

def cij_f(tau0, sigi, sigj, segmenti, segmentj):
    """
    Transmission probability function.
    
    Parameters:
    -----------
    tau0 : float
        Initial optical depth
    sigi : float
        Cross section of first region (cm^-1)
    sigj : float
        Cross section of second region (cm^-1)
    segmenti : float
        First segment length (cm)
    segmentj : float
        Second segment length (cm)
    
    Returns:
    --------
    float : Transmission probability
    
    (c) 2008 Alain Hebert, Ecole Polytechnique de Montreal
    Translated to Python 2026
    """
    if sigi != 0 and sigj != 0:
        return (akin(3, tau0) - akin(3, tau0 + sigi * segmenti) - akin(3, tau0 + sigj * segmentj) +
                akin(3, tau0 + sigi * segmenti + sigj * segmentj)) / (sigi * sigj)
    elif sigi == 0 and sigj != 0:
        return (akin(2, tau0) - akin(2, tau0 + sigj * segmentj)) * segmenti / sigj
    elif sigi != 0 and sigj == 0:
        return (akin(2, tau0) - akin(2, tau0 + sigi * segmenti)) * segmentj / sigi
    else:
        return akin(1, tau0) * segmenti * segmentj

def indpos(i, j):
    """
    Compute the index position in a packed symmetric matrix.
    Maps (i,j) to a linear index for upper triangular storage.
    """
    m = max(i, j)
    n = min(i, j)
    return int(m * (m - 1) // 2 + n)

def tij2d(track, sigt):
    """
    Integration of the collision, escape and transmission probabilities
    in unstructured finite 2D geometry.
    
    Parameters:
    -----------
    track : ndarray
        Tracking array from sybt2d
    sigt : array-like
        Total cross sections for each region
    di_f : function
        Collision probability function di_f(sig, seg)
    ei_f : function
        Escape probability function ei_f(tau0, sig, seg)
    cij_f : function
        Transmission probability function cij_f(tau0, sig1, sig3, seg1, seg3)
    akin : function
        Function akin(n, tau0)
    
    Returns:
    --------
    tij : ndarray
        Probability matrix (packed symmetric format)
    
    (c) 2009 Alain Hebert, Ecole Polytechnique de Montreal
    Adapted to Python by Loutre 2026
    """
    nsurf = int(track[0])
    nreg = int(track[1])
    # MATLAB: k=5+track(1)+track(2)+2*track(3)
    k = int(5 + track[0] + track[1] + 2 * track[2])
    tij = np.zeros(int((nreg + nsurf) * (nreg + nsurf + 1) // 2))
    
    # MATLAB: for itrk=1:track(4)
    for itrk in range(int(track[3])):
        # MATLAB: isurf=track(k+2) ; jsurf=track(k+3) ; wei=track(k+4) ; km=track(k+5)
        isurf = int(track[k + 1])  # k+2 in MATLAB -> k+1 in Python
        jsurf = int(track[k + 2])  # k+3 in MATLAB -> k+2 in Python
        wei = track[k + 3]         # k+4 in MATLAB -> k+3 in Python
        km = int(track[k + 4])     # k+5 in MATLAB -> k+4 in Python
        
        kgar = k + 4  # k+5 in MATLAB becomes k+4 in Python (0-based)

        k = k + 5 + km
        irs = isurf
        seg1 = 0.
        sig1 = 0.
        
        # MATLAB: for ixi=1:km
        for ixi in range(1, km + 1):
            irt = irs
            # MATLAB: irs=track(kgar+ixi)  
            # kgar is set so track(kgar+ixi) in MATLAB = track[kgar+ixi] in Python
            irs = int(track[kgar + ixi])
            # MATLAB: seg2=track(k+ixi)
            # k was incremented by 5+km, so track(k+ixi) in MATLAB = track[k+ixi-1] in Python
            seg2 = track[k + ixi - 1]
            # MATLAB: sig2=sigt(irs) - irs is 1-based region number
            sig2 = sigt[irs - 1]
            # MATLAB: irs=irs+nsurf
            irs = irs + nsurf
            
            iij = indpos(irs, irs)
            # MATLAB: tij(iij)=tij(iij)+2.0*wei*di_f(sig2,seg2)
            tij[iij - 1] += 2.0 * wei * di_f(sig2, seg2)  # iij-1 for 0-based indexing
            tau0 = 0.
            
            # MATLAB: for ixj=ixi:km
            for ixj in range(ixi, km + 1):
                # MATLAB: jrs=track(kgar+ixj) ; seg3=track(k+ixj) ; sig3=sigt(jrs)
                # kgar is set so track(kgar+ixj) in MATLAB = track[kgar+ixj] in Python
                jrs = int(track[kgar + ixj])
                # k was incremented, so track(k+ixj) in MATLAB = track[k+ixj-1] in Python
                seg3 = track[k + ixj - 1]
                # Handle potential invalid region indices (workaround for tracking issues)
                if jrs < 1 or jrs > nreg:
                    # Use region 1's cross section as default
                    sig3 = sigt[0] if len(sigt) > 0 else 0.0
                else:
                    sig3 = sigt[jrs - 1]  # Get cross section BEFORE offsetting jrs
                # MATLAB: jrs=jrs+nsurf (done AFTER getting sig3)
                jrs = jrs + nsurf
                
                iij = indpos(irt, jrs)
                # MATLAB: if irt <= nsurf
                if irt <= nsurf:
                    # MATLAB: tij(iij)=tij(iij)+wei*ei_f(tau0,sig3,seg3)
                    tij[iij - 1] += wei * ei_f(tau0, sig3, seg3)
                else:
                    # MATLAB: wi3=cij_f(tau0,sig1,sig3,seg1,seg3)
                    wi3 = cij_f(tau0, sig1, sig3, seg1, seg3)
                    # MATLAB: if jrs == irt, wi3=2.0*wi3; end
                    if jrs == irt:
                        wi3 = 2.0 * wi3
                    # MATLAB: tij(iij)=tij(iij)+wei*wi3
                    tij[iij - 1] += wei * wi3
                
                tau0 = tau0 + seg3 * sig3
            
            iij = indpos(irt, jsurf)
            # MATLAB: if irt <= nsurf
            if irt <= nsurf:
                # MATLAB: wi3=akin(3,tau0)
                wi3 = akin(3, tau0)
                # MATLAB: if isurf == jsurf, wi3=2.0*wi3; end
                if isurf == jsurf:
                    wi3 = 2.0 * wi3
                # MATLAB: tij(iij)=tij(iij)+wei*wi3
                tij[iij - 1] += wei * wi3
            else:
                # MATLAB: tij(iij)=tij(iij)+wei*ei_f(tau0,sig1,seg1)
                tij[iij - 1] += wei * ei_f(tau0, sig1, seg1)
            
            seg1 = seg2
            sig1 = sig2
        
        # MATLAB: iij=indpos(irs,jsurf) ; tij(iij)=tij(iij)+wei*ei_f(0.0,sig1,seg1)
        iij = indpos(irs, jsurf)
        tij[iij - 1] += wei * ei_f(0.0, sig1, seg1)
        k = k + km
    
    # MATLAB: tij(:)=tij(:).*track(5)^2
    tij[:] = tij[:] * track[4]**2
    
    return tij


def sybt2d(a, rad, nangle, ngauss):
    """
    Produce a Gauss-Jacobi tracking in 2D square pincell geometry.
    
    Faithful translation of sybt2d.m from A. Hébert's MATLAB code.
    
    Parameters:
    -----------
    a : float
        Side-length of the square pincell.
    rad : array-like
        Radii of concentric circles defining regions (excluding center at 0).
    nangle : int
        Number of angles for tracking.
    ngauss : int
        Number of Gauss quadrature points (1-6).
    
    Returns:
    --------
    track : ndarray
        Tracking array containing all geometry and tracking information.
    
    (c) 2009 Alain Hebert, Ecole Polytechnique de Montreal
    Translated to Python 2026
    """
    # MATLAB: nreg=1+size(rad,2) ; radd = [0. rad] ; na2=2*nangle ;
    rad = np.asarray(rad)
    nreg = 1 + rad.size
    radd = np.concatenate(([0.], rad))
    na2 = 2 * nangle
    
    # Gauss-Jacobi quadrature parameters
    if ngauss == 1:
        alp = np.array([0.6666666667])
        pwr = np.array([0.5])
        zx = np.array([0.])
        wx = np.array([2.])
    elif ngauss == 2:
        alp = np.array([0.3550510257, 0.8449489743])
        pwr = np.array([0.1819586183, 0.3180413817])
        zx = np.array([-0.577350259, 0.577350259])
        wx = np.array([1., 1.])
    elif ngauss == 3:
        alp = np.array([0.2123405382, 0.5905331356, 0.9114120405])
        pwr = np.array([0.0698269799, 0.2292411064, 0.2009319137])
        zx = np.array([-0.774596691, 0., 0.774596691])
        wx = np.array([0.555555556, 0.888888889, 0.555555556])
    elif ngauss == 4:
        alp = np.array([0.1397598643, 0.4164095676, 0.7231569864, 0.9428958039])
        pwr = np.array([0.0311809710, 0.1298475476, 0.2034645680, 0.1355069134])
        zx = np.array([-0.861136317, -0.339981049, 0.339981049, 0.861136317])
        wx = np.array([0.347854853, 0.652145147, 0.652145147, 0.347854853])
    elif ngauss == 5:
        alp = np.array([0.0985350858, 0.3045357266, 0.5620251898, 0.8019865821, 0.9601901429])
        pwr = np.array([0.0157479145, 0.0739088701, 0.1463869871, 0.1671746381, 0.0967815902])
        zx = np.array([-0.906179845, -0.538469315, 0., 0.538469315, 0.906179845])
        wx = np.array([0.236926883, 0.478628665, 0.568888843, 0.478628665, 0.236926883])
    elif ngauss == 6:
        alp = np.array([0.0730543287, 0.2307661380, 0.4413284812, 0.6630153097, 0.8519214003, 0.9706835728])
        pwr = np.array([0.0087383018, 0.0439551656, 0.0986611509, 0.1407925538, 0.1355424972, 0.0723103307])
        zx = np.array([-0.932469487, -0.661209404, -0.238619193, 0.238619193, 0.661209404, 0.932469487])
        wx = np.array([0.171324492, 0.360761583, 0.467913926, 0.467913926, 0.360761583, 0.171324492])
    else:
        raise ValueError('Invalid number of Gauss-Jacobi points.')
    
    if 2.0 * radd[nreg - 1] > a:
        raise ValueError('A radius is greater than half a side.')
    
    # MATLAB: track_w=zeros(1,9+nreg+4*na2*(2+ngauss*nreg*(5+2*(2*nreg-1))))
    track_w = np.zeros(9 + nreg + 4 * na2 * (2 + ngauss * nreg * (5 + 2 * (2 * nreg - 1))))
    # MATLAB: kstart=9+nreg+8*na2 ; track_w(1:3)=[4, nreg, 4*na2] ; zn1=0.
    kstart = 9 + nreg + 8 * na2
    track_w[0:3] = [4, nreg, 4 * na2]  # MATLAB track_w(1:3) = Python track_w[0:3]
    zn1 = 0.
    
    # MATLAB: ao2=a/2. ; wa=2./real(na2) ; track_w(6:9)=a ; vol=a*a ;
    ao2 = a / 2.
    wa = 2. / na2
    track_w[5:9] = a  # MATLAB track_w(6:9) = Python track_w[5:9]
    vol = a * a
    
    # MATLAB: for jjj=nreg:-1:1
    for jjj in range(nreg, 0, -1):
        # MATLAB: r2=pi*radd(jjj)^2 ; track_w(9+jjj)=vol-r2 ; vol=r2 ;
        r2 = np.pi * radd[jjj - 1]**2  # MATLAB radd(jjj) = Python radd[jjj-1]
        track_w[9 + jjj - 1] = vol - r2  # MATLAB track_w(9+jjj) = Python track_w[9+jjj-1]
        vol = r2
    
    # Track generation
    k = kstart
    # MATLAB: for ia=1:na2
    for ia in range(1, na2 + 1):
        # MATLAB: za=(2.0*real(ia)-1.)/real(na2)-1. ; phi=0.25*pi*(za+1.) ;
        za = (2.0 * ia - 1.) / na2 - 1.
        phi = 0.25 * np.pi * (za + 1.)
        # MATLAB: zn1=zn1+sin(phi)*wa ; si=sin(phi) ; co=cos(phi) ; ta=si/co ;
        zn1 += np.sin(phi) * wa
        si = np.sin(phi)
        co = np.cos(phi)
        ta = si / co
        
        # MATLAB: track_w(9+nreg+ia)=si ; track_w(9+nreg+4*na2+ia)=co ;
        track_w[9 + nreg + ia - 1] = si
        track_w[9 + nreg + 4 * na2 + ia - 1] = co
        # MATLAB: track_w(9+nreg+na2+ia)=co ; track_w(9+nreg+5*na2+ia)=-si ;
        track_w[9 + nreg + na2 + ia - 1] = co
        track_w[9 + nreg + 5 * na2 + ia - 1] = -si
        
        if phi <= 0.25 * np.pi:
            # MATLAB: track_w(9+nreg+2*na2+ia)=co ; track_w(9+nreg+6*na2+ia)=-si ;
            track_w[9 + nreg + 2 * na2 + ia - 1] = co
            track_w[9 + nreg + 6 * na2 + ia - 1] = -si
            # MATLAB: track_w(9+nreg+3*na2+ia)=si ; track_w(9+nreg+7*na2+ia)=co ;
            track_w[9 + nreg + 3 * na2 + ia - 1] = si
            track_w[9 + nreg + 7 * na2 + ia - 1] = co
            # MATLAB: jsu=4 ; x1=0. ; xlim=a ; dlim=ao2*co+(ao2-xlim)*si ;
            jsu = 4
            x1 = 0.
            xlim = a
            dlim = ao2 * co + (ao2 - xlim) * si
        else:
            # MATLAB: track_w(9+nreg+2*na2+ia)=co ; track_w(9+nreg+6*na2+ia)=si ;
            track_w[9 + nreg + 2 * na2 + ia - 1] = co
            track_w[9 + nreg + 6 * na2 + ia - 1] = si
            # MATLAB: track_w(9+nreg+3*na2+ia)=si ; track_w(9+nreg+7*na2+ia)=-co ;
            track_w[9 + nreg + 3 * na2 + ia - 1] = si
            track_w[9 + nreg + 7 * na2 + ia - 1] = -co
            # MATLAB: jsu=2 ; x1=a/ta ; xlim=0.5*(a+a/ta) ; dlim=0. ;
            jsu = 2
            x1 = a / ta
            xlim = 0.5 * (a + a / ta)
            dlim = 0.
        
        # MATLAB: for k0=nreg:-1:1
        for k0 in range(nreg, 0, -1):
            # MATLAB: km=nreg-k0+1 ; x2=min(xlim,xlim-(radd(k0)-dlim)/si) ;
            km = nreg - k0 + 1
            x2 = min(xlim, xlim - (radd[k0 - 1] - dlim) / si)  # radd(k0) = radd[k0-1]
            
            # MATLAB: if ((x1 < xlim) && (phi <= 0.25*pi)) || ((x1 < x2) && (phi > 0.25*pi))
            if ((x1 < xlim) and (phi <= 0.25 * np.pi)) or ((x1 < x2) and (phi > 0.25 * np.pi)):
                l3 = k
                # MATLAB: vap=zeros(1,nreg) ;
                vap = np.zeros(nreg)
                
                # MATLAB: for ix=1:ngauss
                for ix in range(1, ngauss + 1):
                    # MATLAB: if k0 == nreg
                    if k0 == nreg:
                        # MATLAB: s=0.5*(x2-x1)*si*wx(ix) ;
                        s = 0.5 * (x2 - x1) * si * wx[ix - 1]  # wx(ix) = wx[ix-1]
                        # MATLAB: x=x1+0.5*(x2-x1)*(1.0+zx(ix)) ;
                        x = x1 + 0.5 * (x2 - x1) * (1.0 + zx[ix - 1])
                    else:
                        # Flurig change of variable.
                        # MATLAB: s=2.*(x2-x1)*si*pwr(ix) ;
                        s = 2. * (x2 - x1) * si * pwr[ix - 1]  # pwr(ix) = pwr[ix-1]
                        # MATLAB: x=x1+(x2-x1)*alp(ix)^2 ;
                        x = x1 + (x2 - x1) * alp[ix - 1]**2
                    
                    # MATLAB: track_w(k+1:k+5)=[ia, 1, jsu, s*wa/4., 2*km-1] ;
                    track_w[k:k + 5] = [ia, 1, jsu, s * wa / 4., 2 * km - 1]  # track_w(k+1) = track_w[k]
                    # MATLAB: track_w(k+6:k+2*km+4)=abs(km-1:-1:1-km)+1+nreg-km ;
                    track_w[k + 5:k + 2 * km + 4] = np.abs(np.arange(km - 1, -km, -1)) + 1 + nreg - km
                    # MATLAB: k=k+5+(2*km-1) ;
                    k += 5 + (2 * km - 1)
                    
                    c = ao2 * si - (ao2 - x) * co
                    d = (ao2 * co + (ao2 - x) * si)**2
                    sumtrk = 0.
                    
                    # MATLAB: for kk=nreg:-1:k0+1
                    for kk in range(nreg, k0, -1):
                        corde = np.sqrt(rad[kk - 2]**2 - d)  # rad(kk-1) in MATLAB; kk-2 in Python
                        delta = c - corde
                        sumtrk += delta
                        # MATLAB: track_w(k+nreg-kk+1)=del
                        track_w[k + nreg - kk] = delta
                        vap[kk - 1] += delta * s  # vap(kk) in MATLAB; kk-1 in Python
                        c = corde
                    
                    if km != 1:
                        delta = 2.0 * corde
                        # MATLAB: track_w(k+km)=del
                        track_w[k + km - 1] = delta
                        vap[k0 - 1] += delta * s  # vap(k0) in MATLAB; k0-1 in Python
                        # MATLAB: sum(track_w(k+km-1:-1:k+2))
                        if km > 2:
                            sumtrk += delta + np.sum(track_w[k + 1:k + km - 1])
                            # MATLAB: track_w(k+km+1:k+2*km-2)=track_w(k+km-1:-1:k+2)
                            track_w[k + km:k + 2 * km - 2] = track_w[k + km - 2:k:-1]
                            # MATLAB: vap(k0+1:k0+km-2)=vap(k0+1:k0+km-2)+track_w(k+km-1:-1:k+2).*s
                            vap[k0:k0 + km - 2] += track_w[k + km - 2:k:-1] * s
                        else:
                            sumtrk += delta
                    
                    k += 2 * km - 1
                    
                    if phi <= 0.25 * np.pi:
                        delta = x / co - sumtrk
                    else:
                        delta = a / si - sumtrk
                    
                    # MATLAB: track_w(k)=del
                    track_w[k - 1] = delta
                    vap[nreg - 1] += delta * s  # vap(nreg) in MATLAB; nreg-1 in Python
                
                # Volume normalization
                if k0 < nreg:
                    dlim1 = ao2 * co + (ao2 - x2) * si
                    dlim2 = ao2 * co + (ao2 - x1) * si
                    vw1 = 0.
                    sumvap = 0.
                    
                    # MATLAB: for i=k0:nreg-1
                    for i in range(k0, nreg):  # k0 is 1-based MATLAB, range gives k0-1 to nreg-1 in Python
                        sumvap += vap[i - 1]
                        rw = rad[i - 1]  # rad(i) in MATLAB; i-1 in Python
                        vex1 = rw * rw * np.arccos(dlim1 / rw) - dlim1 * np.sqrt(rw * rw - dlim1 * dlim1)
                        if rw > dlim2:
                            vex1 -= (rw * rw * np.arccos(dlim2 / rw) - dlim2 * np.sqrt(rw * rw - dlim2 * dlim2))
                        vap[i - 1] = (vex1 - vw1) / vap[i - 1]
                        vw1 = vex1
                    
                    vex1 = 0.5 * (a * si - (a - x1 - x2) * co) * (x2 - x1) * si
                    if phi <= 0.25 * np.pi:
                        vex2 = 0.5 * ta * (x2 * x2 - x1 * x1) - vex1
                    else:
                        vex2 = (x2 - x1) * a - vex1
                    
                    vex1 = (vex1 - 0.5 * vw1) / (vex1 - 0.5 * sumvap)
                    vex2 = (vex2 - 0.5 * vw1) / (vex2 - 0.5 * sumvap)
                    
                    # MATLAB: for ix=1:ngauss
                    # l3 currently points to the start of the tracks for this k0/ia
                    for ix in range(1, ngauss + 1):
                        # MATLAB: l3=l3+5 ; km=(track_w(l3)+1)/2
                        l3 += 5  # In MATLAB, l3 now points to position containing nseg
                        # In Python, l3 now = k+5, but nseg is at k+4, so use l3-1
                        km = int((track_w[l3 - 1] + 1) / 2)  # track_w(l3) = track_w[l3-1]
                        # MATLAB: l3=l3+2*km-1
                        l3 += 2 * km - 1  # Now l3 points to last region index
                        # MATLAB: track_w(l3+km)=track_w(l3+km)*vap(k0)
                        track_w[l3 + km - 1] *= vap[k0 - 1]  # track_w(l3+km) = track_w[l3+km-1]
                        # MATLAB: fact=vap(k0+1:k0+km-2)
                        if km > 2:
                            fact = vap[k0:k0 + km - 2]  # vap(k0+1:k0+km-2) = vap[k0:k0+km-2]
                            # MATLAB: track_w(l3+km-1:-1:l3+2)=track_w(l3+km-1:-1:l3+2).*fact
                            track_w[l3 + 1:l3 + km - 1] *= fact[::-1]  # track_w(l3+2:l3+km-1) = track_w[l3+1:l3+km-1]
                            # MATLAB: track_w(l3+km+1:l3+2*km-2)=track_w(l3+km+1:l3+2*km-2).*fact
                            track_w[l3 + km:l3 + 2 * km - 2] *= fact  # track_w(l3+km+1:l3+2*km-2) = track_w[l3+km:l3+2*km-2]
                        # MATLAB: track_w(l3+1)=track_w(l3+1)*vex1
                        track_w[l3] *= vex1  # track_w(l3+1) = track_w[l3]
                        # MATLAB: track_w(l3+2*km-1)=track_w(l3+2*km-1)*vex2
                        track_w[l3 + 2 * km - 2] *= vex2  # track_w(l3+2*km-1) = track_w[l3+2*km-2]
                        # MATLAB: l3=l3+2*km-1
                        l3 += 2 * km - 1
                
                track_w[3] += ngauss
                x1 = x2
    
    track_w[4] = 1. / np.sqrt(0.25 * np.pi * zn1)
    kend = k
    
    # Apply symmetries
    track_w[kend:2 * kend - kstart] = track_w[kstart:kend]
    k = kend
    # MATLAB: for itrk=1:track_w(4)
    for itrk in range(int(track_w[3])):
        # MATLAB: track_w(k+1)=na2+track_w(k+1)
        track_w[k] = na2 + track_w[k]
        # MATLAB: nseg=track_w(k+5)
        nseg = int(track_w[k + 4])
        # MATLAB: if track_w(k+3) == 2
        if track_w[k + 2] == 2:
            track_w[k + 1] = 4
            track_w[k + 2] = 3
        # MATLAB: elseif track_w(k+3) == 4
        elif track_w[k + 2] == 4:
            track_w[k + 2] = 3
        # MATLAB: track_w(k+6+nseg:k+5+2*nseg)=track_w(k+5+2*nseg:-1:k+6+nseg)
        track_w[k + 5 + nseg:k + 5 + 2 * nseg] = track_w[k + 5 + 2 * nseg - 1:k + 5 + nseg - 1:-1]
        # MATLAB: k=k+5+2*track_w(k+5)
        k += 5 + 2 * int(track_w[k + 4])
    
    track_w[3] = 2 * track_w[3]
    kend = k
    track_w[kend:2 * kend - kstart] = track_w[kstart:kend]
    k = kend
    # MATLAB: for itrk=1:track_w(4)
    for itrk in range(int(track_w[3])):
        # MATLAB: track_w(k+1)=2*na2+track_w(k+1)
        track_w[k] = 2 * na2 + track_w[k]
        # MATLAB: if track_w(k+3) == 4
        if track_w[k + 2] == 4:
            track_w[k + 1] = 4
            track_w[k + 2] = 2
        # MATLAB: elseif track_w(k+3) == 2
        elif track_w[k + 2] == 2:
            track_w[k + 1] = 3
            track_w[k + 2] = 4
        # MATLAB: elseif track_w(k+2) == 1 && track_w(k+3) == 3
        elif track_w[k + 1] == 1 and track_w[k + 2] == 3:
            track_w[k + 1] = 3
            track_w[k + 2] = 2
        # MATLAB: elseif track_w(k+2) == 4 && track_w(k+3) == 3
        elif track_w[k + 1] == 4 and track_w[k + 2] == 3:
            track_w[k + 1] = 1
            track_w[k + 2] = 2
        # MATLAB: k=k+5+2*track_w(k+5)
        k += 5 + 2 * int(track_w[k + 4])
    
    track_w[3] = 2 * track_w[3]
    
    if k > 9 + nreg + 4 * na2 * (2 + ngauss * nreg * (5 + 2 * (2 * nreg - 1))):
        raise RuntimeError('Tracking overflow.')
    
    track = track_w[:k]
    return track


def sybrhl(track, sigt, pij):
    """
    Stamm'ler normalization of collision, escape and transmission
    probabilities in unstructured finite 2D geometry.
    
    This function performs an iterative normalization procedure to correct
    probability matrices computed by tij2d, ensuring proper conservation
    and balance equations.
    
    Parameters:
    -----------
    track : ndarray
        Tracking array from sybt2d
    sigt : array-like
        Total macroscopic cross sections for each region (cm^-1)
    pij : ndarray
        Probability matrix from tij2d (packed symmetric format)
    
    Returns:
    --------
    pij : ndarray
        Normalized probability matrix (packed symmetric format)
    
    Algorithm:
    ----------
    Uses iterative method with acceleration by residual minimization to
    compute normalization weights that satisfy neutron balance equations.
    
    (c) 2008 Alain Hebert, Ecole Polytechnique de Montreal
    Translated to Python 2026
    """
    # Parameters
    epscon = 1.0e-6   # Convergence criterion
    nitmax = 20       # Maximum iterations
    cptlb = 3         # Cycle length before acceleration
    cptac = 3         # Acceleration cycles
    
    # Extract geometry info
    nsurf = int(track[0])
    nreg = int(track[1])
    nun = nsurf + nreg  # Total number of surfaces + regions
    
    # Initialize weights array: weig[:, 0]=old, weig[:, 1]=current, weig[:, 2]=new
    weig = np.zeros((nun, 3))
    weig[:, 1:3] = 0.5
    
    # MATLAB: g=track(6:5+nsurf+nreg)
    # Python: track[5:5+nsurf+nreg]
    g = track[5:5 + nsurf + nreg].copy()
    # MATLAB: g(1:nsurf)=g(1:nsurf)./4.
    g[0:nsurf] = g[0:nsurf] / 4.
    # MATLAB: g(nsurf+1:nun)=g(nsurf+1:nun).*sigt
    g[nsurf:nun] = g[nsurf:nun] * np.asarray(sigt)
    
    # Copy pij for working
    pij_w = pij.copy()
    
    # Reconstruct real probabilities
    # MATLAB: iij=nsurf*(nsurf+1)/2
    iij = nsurf * (nsurf + 1) // 2
    # MATLAB: for ir=1:nreg
    for ir in range(1, nreg + 1):
        # MATLAB: pij_w(iij+1:iij+nsurf)=pij_w(iij+1:iij+nsurf).*sigt(ir)
        pij_w[iij:iij + nsurf] *= sigt[ir - 1]  # Python 0-based for array access
        # MATLAB: iij=iij+nsurf+ir
        iij += nsurf + ir
    
    # MATLAB: iij=nsurf*(nsurf+1)/2
    iij = nsurf * (nsurf + 1) // 2
    # MATLAB: for jr=1:nreg
    for jr in range(1, nreg + 1):
        # MATLAB: iij=iij+nsurf
        iij += nsurf
        # MATLAB: pij_w(iij+1:iij+jr)=pij_w(iij+1:iij+jr).*(sigt(1:jr)*sigt(jr))
        pij_w[iij:iij + jr] *= np.asarray(sigt[0:jr]) * sigt[jr - 1]
        # MATLAB: iij=iij+jr
        iij += jr
    
    # Main iteration loop
    nit = 0
    totcon = 9999.
    
    print(f"[sybrhl] Starting normalization iterations...")
    print(f"[sybrhl] Initial totcon: {totcon}, epscon: {epscon}")
    
    while totcon >= epscon:
        nit += 1
        if nit > nitmax:
            raise RuntimeError('Weights not converged.')
        
        # MATLAB: for ir=1:nun
        for ir in range(1, nun + 1):
            # Compute weighted sum
            # MATLAB: wfspad=g(ir)+pij_w(indpos(ir,ir))*weig(ir,3)
            wfspad = g[ir - 1] + pij_w[indpos(ir, ir) - 1] * weig[ir - 1, 2]
            
            # MATLAB: wfspad=wfspad-sum(weig(1:nun,3).*pij_w(indpos(ir,1:nun))')
            sum_term = 0.
            for j in range(1, nun + 1):
                sum_term += weig[j - 1, 2] * pij_w[indpos(ir, j) - 1]
            wfspad -= sum_term
            
            # MATLAB: wfsp=pij_w(indpos(ir,ir))+sum(pij_w(indpos(ir,1:nun)))
            wfsp = pij_w[indpos(ir, ir) - 1]
            for j in range(1, nun + 1):
                wfsp += pij_w[indpos(ir, j) - 1]
            
            if wfsp != 0:
                weig[ir - 1, 2] = wfspad / wfsp
        
        # Acceleration by residual minimization
        zmu = 1.0
        # MATLAB: if mod(nit-1,cptac+cptlb) >= cptac
        if (nit - 1) % (cptac + cptlb) >= cptac:
            # MATLAB: nom=sum((weig(:,2)-weig(:,1)).*(weig(:,3)-2.*weig(:,2)+weig(:,1)))
            nom = np.sum((weig[:, 1] - weig[:, 0]) * (weig[:, 2] - 2. * weig[:, 1] + weig[:, 0]))
            # MATLAB: denom=sum((weig(:,3)-2.*weig(:,2)+weig(:,1)).^2)
            denom = np.sum((weig[:, 2] - 2. * weig[:, 1] + weig[:, 0])**2)
            
            if denom != 0.0:
                zmu = -nom / denom
            
            if (zmu > 5.0) or (zmu < 0.0):
                zmu = 1.0
            
            # MATLAB: weig(1:nun,3)=weig(1:nun,2)+zmu*(weig(1:nun,3)-weig(1:nun,2))
            weig[0:nun, 2] = weig[0:nun, 1] + zmu * (weig[0:nun, 2] - weig[0:nun, 1])
            # MATLAB: weig(1:nun,2)=weig(1:nun,1)+zmu*(weig(1:nun,2)-weig(1:nun,1))
            weig[0:nun, 1] = weig[0:nun, 0] + zmu * (weig[0:nun, 1] - weig[0:nun, 0])
        
        # Calculate convergence metric
        # MATLAB: totcon=max(abs(weig(:,3)-weig(:,2))./weig(:,3))
        with np.errstate(divide='ignore', invalid='ignore'):
            rel_diff = np.abs(weig[:, 2] - weig[:, 1]) / weig[:, 2]
            rel_diff = np.where(np.isfinite(rel_diff), rel_diff, 0.)
            totcon = np.max(rel_diff)
        
        # Update weights
        # MATLAB: weig(1:nun,1)=weig(1:nun,2) ; weig(1:nun,2)=weig(1:nun,3)
        weig[0:nun, 0] = weig[0:nun, 1]
        weig[0:nun, 1] = weig[0:nun, 2]
    
    print(f"[sybrhl] Converged after {nit} iterations, final totcon: {totcon}")
    print(f"[sybrhl] Final weights (first 5): {weig[0:5, 1]}")
    print(f"[sybrhl] pij sum BEFORE final renormalization: {np.sum(pij)}")
    
    # Renormalize "pij" symmetric matrix
    # MATLAB: iprb=0
    iprb = 0
    # MATLAB: for ir=1:nun
    for ir in range(1, nun + 1):
        # MATLAB: pij(iprb+1:iprb+ir)=pij(iprb+1:iprb+ir).*(weig(ir,1)+weig(1:ir,1)')
        # After loop ends: weig[:,0] has old iteration, weig[:,1] has final converged
        # MATLAB column 1 = Python column 1 (after shift), so use weig[:,1]
        # BUT: MATLAB saved column 2 to column 1, so final values are in column 1
        # Actually: before the renorm loop, we did: weig[:,1] = weig[:,2]
        # So converged values from column 2 are now in column 1
        pij[iprb:iprb + ir] *= (weig[ir - 1, 1] + weig[0:ir, 1])
        # MATLAB: iprb=iprb+ir
        iprb += ir
    
    print(f"[sybrhl] pij sum AFTER final renormalization: {np.sum(pij)}")
    print(f"[sybrhl] Sample weight sums: {weig[0,1] + weig[0,1]}, {weig[4,1] + weig[0,1]}")
    
    return pij
