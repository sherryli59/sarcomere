import actin 

if __name__ == "__main__":

    n_filaments = 100
    n_bundles = 10
    n_crosslinks = 50
    boxL = 30.
    sarcomere = actin.SarcomereModel(n_filaments, n_bundles, n_crosslinks, boxL, boxL)
    sarcomere.run_simulation(1e-3, 1., 10000)

    #actin.plot_system()
    #plt.show()