# ---- Surface plot ----
def surfaceplotX(U, Uexact, tspan, xspan): 
    #meshgrid : Return coordinate matrices from coordinate vectors
    X, T = np.meshgrid(tspan, xspan)
    fig = plt.figure(figsize=plt.figaspect(0.5))
    #fig = plt.figure(figsize = (7.5,5.5))
    fig.suptitle("----- Method -------------------------- Error |(uexact - u)| -----", fontsize=12)
    
    ax1 = fig.add_subplot(1, 2, 1,projection='3d')
    surf = ax1.plot_surface(X, T, U, linewidth=0, cmap=cm.jet, antialiased=True)
    ax1.set_xlabel('Time')
    ax1.set_ylabel('Space')
    ax1.set_zlabel('U')
    fig.savefig(path.join("plot_method{0}.png".format(count)))
    
    ax2 = fig.add_subplot(1, 2, 2, projection='3d')
    surf2 = ax2.plot_surface(X, T, abs(Uexact-U), linewidth=0, cmap=cm.jet, antialiased=True)
    ax2.set_xlabel('Time')
    ax2.set_ylabel('Space')
    ax2.set_zlabel('U')
    
    plt.tight_layout()
    ax1.view_init(30,230)
    ax2.view_init(30,230)
    fig.savefig(path.join("plot_exact_error{0}.png".format(count)),dpi=900)
    
    fig, ax = plt.subplots()
    CS = ax.contour(X, T, abs(Uexact-U), 12)

    plt.draw()
