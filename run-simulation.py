import actin_gpu as actin 
import matplotlib.pyplot as plt
if __name__ == "__main__":

    # n_filaments = 100
    # n_bundles = 10
    # n_crosslinks = 50
    # boxL = 30.
    n_filaments = 10
    n_bundles = 2
    n_crosslinks = 5
    boxL = 10.
    sarcomere = actin.SarcomereModel(n_filaments, n_bundles, n_crosslinks, boxL, boxL)
    sarcomere.run_simulation(1e-3, 1., 200000, make_animation=False)
