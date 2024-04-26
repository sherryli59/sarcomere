import numpy as np
import matplotlib.pyplot as plt

class Filaments():
    def __init__(self, n_filaments, Lx, Ly):
        self.n_filaments = n_filaments
        self.Lx = Lx
        self.Ly = Ly
        self.filament_length = 3.0
        self.xs = np.zeros((n_filaments, 2))
        self.xs[:,0] = np.random.uniform(0, Lx, n_filaments)
        self.xs[:,1] = np.random.uniform(0, Ly, n_filaments)
        self.thetas = np.random.uniform(0, 2*np.pi, n_filaments)

class MysosinBundle():
    def __init__(self, n_bundles, Lx, Ly):
        self.n_bundles = n_bundles
        self.Lx = Lx
        self.Ly = Ly
        self.xs = np.zeros((n_bundles, 2))
        self.xs[:,0] = np.random.uniform(0, Lx, n_bundles)
        self.xs[:,1] = np.random.uniform(0, Ly, n_bundles)
        self.radius = 1.0

class AlphaActinin():
    def __init__(self, n_crosslinks, Lx, Ly):
        self.n_crosslinks = n_crosslinks
        self.Lx = Lx  
        self.Ly = Ly
        self.xs = np.zeros((n_crosslinks, 2))
        self.xs[:,0] = np.random.uniform(0, Lx, n_crosslinks)
        self.xs[:,1] = np.random.uniform(0, Ly, n_crosslinks)
        self.radius = 0.1


class SarcomereModel():
    def __init__(self, n_filaments, n_bundles, n_crosslinks, Lx, Ly):
        self.actin = Filaments(n_filaments, Lx, Ly)
        self.myosin = MysosinBundle(n_bundles, Lx, Ly)
        self.alpha_actinin = AlphaActinin(n_crosslinks, Lx, Ly)
        self.Lx = Lx
        self.Ly = Ly
        self.e_am = -10.0
        self.e_al = -1.0
        self.f_myosin = 1.0
        self.energy = self.compute_energy()
    
    def compute_energy(self):
        return self.myosin_self_energy()  + self.actin_myosin_energy() + self.actin_alpha_energy()

    def myosin_self_energy(self):
        n_bundles = self.myosin.n_bundles
        xs = self.myosin.xs
        radius = self.myosin.radius
        energy = 0.0
        for i in range(n_bundles):
            for j in range(i+1, n_bundles):
                rij = np.linalg.norm(xs[i] - xs[j])
                if rij < 2*radius:
                    return np.inf
                else:
                    pass
        return energy
    
    def actin_myosin_energy(self):
        segments = self.actin.filament_length*np.array([np.cos(self.actin.thetas), np.sin(self.actin.thetas)]).T
        actin_endpts = np.concatenate([self.actin.xs + 0.5*segments, self.actin.xs - 0.5*segments], axis=0)
        n_bundles = self.myosin.n_bundles
        energy = 0.0
        for i in range(n_bundles):
            bundle = self.myosin.xs[i]
            for j in range(2*self.actin.n_filaments):
                actin_segment = actin_endpts[j]
                r = bundle - actin_segment
                if np.linalg.norm(r) < self.myosin.radius:
                    energy += self.e_am
                else:
                    pass
        return energy
    
    def actin_alpha_energy(self):
        n_crosslinks = self.alpha_actinin.n_crosslinks
        n_filaments = self.actin.n_filaments
        la = self.actin.filament_length
        ra = self.alpha_actinin.radius
        energy = 0.0
        for i in range(n_crosslinks):
            crosslink = self.alpha_actinin.xs[i]
            for j in range(n_filaments):
                # compute distance between crosslink and actin filament
                thetas = self.actin.thetas[j]
                actin_endpt = self.actin.xs[j] + 0.5 * la * np.array([np.cos(thetas), np.sin(thetas)])
                x_alpha = crosslink - actin_endpt
                x_actin = -la * np.array([np.cos(thetas), np.sin(thetas)])
                dx = x_alpha - x_actin
            
                # check if crosslink is between the two ends of the actin filament
                if np.dot(x_actin, dx) > 0 and np.dot(x_alpha, x_actin) > 0:
                    if np.sqrt(1 - (np.dot(x_actin, x_alpha)/(np.linalg.norm(x_actin)*la))**2) < ra:
                        energy += self.e_al
                    else:
                        pass
                    
        return energy

    #def actin_myosin_force(self):


    def run_simulation(self, dt, D, n_steps, frame_dir="frames", make_animation=False):
        import os
        if not os.path.exists(frame_dir):
            os.mkdir(frame_dir)

        for i in range(n_steps):
            self.mala_step(dt, D)
            if i % 100 == 0:
                fig, ax = self.plot_system()
                plt.savefig("{}/frame_{:04d}.png".format(frame_dir,i))
                plt.close(fig)

        # save frames as an animation
        if make_animation:
            os.system("ffmpeg -r 1 -i {}/frame_%04d.png -vcodec mpeg4 -y movie.mp4".format(frame_dir))

    # TODO: use stokes einstein equation to compute diffusion coefficient for each component
    def mala_step(self, dt, D):
        # update actin filaments
        xs_old = self.actin.xs.copy()
        thetas_old = self.actin.thetas.copy()
        self.actin.xs = self.pbc(self.actin.xs + np.sqrt(2*D*dt)*np.random.randn(self.actin.n_filaments, 2))
        self.actin.thetas += np.sqrt(2*D*dt)*np.random.randn(self.actin.n_filaments)
        delta_energy = self.compute_energy() - self.energy
        if np.random.rand() < np.exp(-delta_energy):
            self.energy+=delta_energy
        else:
            self.actin.xs = xs_old
            self.actin.thetas = thetas_old

        # update myosin bundles
        xs_old = self.myosin.xs.copy()
        self.myosin.xs = self.pbc(self.myosin.xs + np.sqrt(2*D*dt)*np.random.randn(self.myosin.n_bundles, 2))
        delta_energy = self.compute_energy() - self.energy
        if np.random.rand() < np.exp(-delta_energy):
            self.energy += delta_energy
        else:
            self.myosin.xs = xs_old

        # update alpha-actinin
        xs_old = self.alpha_actinin.xs.copy()
        self.alpha_actinin.xs = self.pbc(self.alpha_actinin.xs + np.sqrt(2*D*dt)*np.random.randn(self.alpha_actinin.n_crosslinks, 2))
        delta_energy = self.compute_energy() - self.energy
        if np.random.rand() < np.exp(-delta_energy):
            self.energy += delta_energy
        else:
            self.alpha_actinin.xs = xs_old

    def pbc(self, xs):
        xs[:,0] = xs[:,0] + self.Lx * (xs[:,0] < 0) - self.Lx * (xs[:,0] >= self.Lx)
        xs[:,1] = xs[:,1] + self.Ly * (xs[:,1] < 0) - self.Ly * (xs[:,1] >= self.Ly)
        return xs

    def plot_system(self):
        n_filaments = self.actin.n_filaments
        xs = self.actin.xs
        thetas = self.actin.thetas
        l = self.actin.filament_length
        fig, ax = plt.subplots()
        for i in range(n_filaments):
            x = xs[i,0]
            y = xs[i,1]
            theta = thetas[i]
            ax.plot([x- 0.5*l*np.cos(theta), x + 0.5*l*np.cos(theta)], [y- 0.5*l*np.sin(theta), y + 0.5*l*np.sin(theta)], lw=2.0, color="C0")

        # plot myosin bundles
        n_bundles = self.myosin.n_bundles
        xs = self.myosin.xs
        for i in range(n_bundles):
            x = xs[i,0]
            y = xs[i,1]
            circle = plt.Circle((x, y), self.myosin.radius, color="C1", alpha=0.5)
            ax.add_artist(circle)

        # plot alpha-actinin
        n_crosslinks = self.alpha_actinin.n_crosslinks
        xs = self.alpha_actinin.xs
        for i in range(n_crosslinks):
            x = xs[i,0]
            y = xs[i,1]
            circle = plt.Circle((x, y), self.alpha_actinin.radius, color="C2", alpha=0.5)
            ax.add_artist(circle)

        ax.set_xlim(0, self.Lx)
        ax.set_ylim(0, self.Ly)
        ax.set_aspect("equal")
        return fig, ax
    