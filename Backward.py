'''
Backward method to solve 1D reaction-diffusion equation:
        u_t = k * u_xx
    
with Neumann boundary conditions 
at x=0: u_x(0,t) = 0 = sin(2*np.pi)
at x=L: u_x(L,t) = 0 = sin(2*np.pi)

with L = 1 and initial conditions:
u(x,0) = (1.0/2.0)+ np.cos(2.0*np.pi*x) - (1.0/2.0)*np.cos(3*np.pi*x)

u_x(x,t) = (-4.0*(np.pi**2))np.exp(-4.0*(np.pi**2)*t)*np.cos(2.0*np.pi*x) + 
            (9.0/2.0)*(np.pi**2)*np.exp(-9.0*(np.pi**2)*t)*np.cos(3*np.pi*x))
'''
def BackwardEuler(M, lambd, T = 0.5, L = 1, k = 1):
    #Parameters needed to solve the equation within the explicit method
    # M = GRID POINTS on space interval
    N = (M**2) #GRID POINTS on time interval

    # ---- Length of the wire in x direction ---- 
    x0, xL = 0, L 

    # ----- Spatial discretization step -----
    dx = (xL - x0)/(M-1)

    # ---- Final time ---- 
    t0, tF = 0, T 

    # ----- Time step -----
    dt = (tF - t0)/(N-1)

    # k = 1.0 Diffusion coefficient
    #lambd = dt*k/dx**2
    a = 1 + 2*lambd

    xspan = np.linspace(x0, xL, M)
    tspan = np.linspace(t0, tF, N)

    main_diag = (1 + 2*lambd)*np.ones((1,M))
    off_diag = -lambd*np.ones((1, M-1))
    a = main_diag.shape[1]
    diagonals = [main_diag, off_diag, off_diag]

    #Sparse Matrix diagonals
    A = sparse.diags(diagonals, [0,-1,1], shape=(a,a)).toarray()
    A[0,1] = -2*lambd
    A[M-1,M-2] = -2*lambd

    # --- Initializes matrix U -----
    U = np.zeros((M, N))

    # --- Initial condition -----
    U[:,0] = (1.0/2.0)+ np.cos(2.0*np.pi*xspan) - (1.0/2.0)*np.cos(3*np.pi*xspan)

    # ---- Neumann boundary conditions -----
    leftBC = np.arange(1, N+1)
    f = np.sin(2*leftBC*np.pi) #f = U[0,:] 

    rightBC = np.arange(1, N+1)
    g = np.sin(2*rightBC*np.pi) #g = U[-1,:]

    for i in range(1, N):
        c = np.zeros((M-2,1)).ravel()
        b1 = np.asarray([2*lambd*dx*f[i], 2*lambd*dx*g[i]])
        b1 = np.insert(b1, 1, c)
        b2 = np.array(U[0:M, i-1])
        b = b1 + b2  # Right hand side
        U[0:M, i] = np.linalg.solve(A,b)  # Solve x=A\b
    
    return (U, tspan, xspan)

#U, tspan, xspan = BackwardEuler(M = 14, lambd = 1.0/6.0)
#Uexact, x, t = ExactSolution(M = 14)
#surfaceplot(U, Uexact, tspan, xspan)
