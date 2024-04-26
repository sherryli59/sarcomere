import torch
import matplotlib.pyplot as plt
import numpy as np

def pairwise_vectors(x,y=None,box=None):
    if y is None:
        y = x
    pair_vec = (x.unsqueeze(-2) - y.unsqueeze(-3))
    if box is not None:
        #unsqueeze and expand the dimensions of box to match the dimensions of pair_vec
        while len(box.shape) < len(pair_vec.shape):
            box = box.unsqueeze(0)
        box = box.expand_as(pair_vec)
        pair_vec = pair_vec - box*torch.round(pair_vec/box)
    return pair_vec

def pairwise_distances(x,y=None,box=None,remove_diag=False):
    pair_vec = pairwise_vectors(x,y,box)
    distances = torch.linalg.norm(pair_vec, axis=-1)
    if y is None and remove_diag:
        rem_dims = distances.shape[:-2]
        n = distances.shape[-1]
        distances = distances.flatten(start_dim=-2)[...,1:].view(*rem_dims,n-1, n+1)[...,:-1].reshape(*rem_dims,n, n-1)
    return distances

def point_segment_distance(x, a, b, box=None):
    ab = b - a
    ap = x - a
    ab_norm = torch.linalg.norm(ab, axis=-1)
    ab = ab/ab_norm
    ap_ab = torch.sum(ap*ab, axis=-1)
    ap_ab = torch.clamp(ap_ab, 0.0, ab_norm)
    dist = torch.linalg.norm(ap - ap_ab[...,None]*ab, axis=-1)
    return dist

class Filaments():
    def __init__(self, n_filaments, Lx, Ly,length=3.0, device="cuda"):
        self.n_filaments = n_filaments
        self.Lx = Lx
        self.Ly = Ly
        self.filament_length = length
        rand = torch.rand((n_filaments,3)).to(device)
        self.xs = torch.stack([Lx*rand[...,0], Ly*rand[...,1]], axis=-1)
        self.thetas = 2*torch.pi*rand[...,2]

class MyosinBundle():
    def __init__(self, n_bundles, Lx, Ly, radius=1.0, device="cuda"):
        self.n_bundles = n_bundles
        self.Lx = Lx
        self.Ly = Ly
        rand = torch.rand((n_bundles,2)).to(device)
        self.xs = torch.stack([Lx*rand[...,0], Ly*rand[...,1]], axis=-1)
        self.radius = radius

class AlphaActinin():
    def __init__(self, n_crosslinks, Lx, Ly, radius=0.3, device="cuda"):
        self.n_crosslinks = n_crosslinks
        self.Lx = Lx  
        self.Ly = Ly
        self.xs = torch.zeros((n_crosslinks, 2))
        rand = torch.rand((n_crosslinks,2)).to(device)
        self.xs = torch.stack([Lx*rand[...,0], Ly*rand[...,1]], axis=-1)
        self.radius = radius


class SarcomereModel():
    def __init__(self, n_filaments, n_bundles, n_crosslinks, Lx, Ly, device="cuda"):
        self.actin = Filaments(n_filaments, Lx, Ly, device=device)
        self.myosin = MyosinBundle(n_bundles, Lx, Ly, device=device)
        self.alpha_actinin = AlphaActinin(n_crosslinks, Lx, Ly, device=device)
        self.Lx = Lx
        self.Ly = Ly
        self.box = torch.tensor([Lx, Ly]).to(device)
        self.e_am = -10.0
        self.e_al = -5.0
        self.f_myosin = 1.0
        #self.sarcomeric_structure()
        fig, ax = self.plot_system()
        plt.savefig("initial.png")
        plt.close(fig)
        self.energy = self.compute_energy()
        print("initial energy: ", self.energy)
        print("number of actin filaments: ", self.actin.n_filaments)
        print("number of myosin bundles: ", self.myosin.n_bundles)
        print("number of alpha-actinin crosslinks: ", self.alpha_actinin.n_crosslinks)
    
    def sarcomeric_structure(self):
        radius = self.myosin.radius*0.8
        unit_length = 2*(radius+self.actin.filament_length)
        n_cols = self.Lx//unit_length
        n_rows = self.myosin.n_bundles//n_cols
        self.myosin.n_bundles = int(n_cols*n_rows)
        x_coord = 0.5*unit_length+unit_length*torch.arange(n_cols).to(self.myosin.xs.device)
        
        am_dist = radius + self.actin.filament_length/2
        actin_x_coord = torch.cat([x_coord + am_dist, x_coord - am_dist], axis=-1)
        n_actin_rows = self.actin.n_filaments//(n_cols*2)
        n_filaments_per_bundle = int(n_actin_rows//n_rows)
        self.actin.n_filaments = int(n_filaments_per_bundle*2*self.myosin.n_bundles)

        y_coord = self.Ly/n_rows*(torch.arange(n_rows).to(self.myosin.xs.device)+0.5)
        
        self.myosin.xs = torch.meshgrid(x_coord,y_coord)
        self.myosin.xs = torch.stack(self.myosin.xs, axis=-1).reshape(-1,2) 
        #introduce offset to every other column
        self.myosin.xs[1::2,1] += radius/n_filaments_per_bundle
        
        actin_y_coord = y_coord[:,None] + radius*torch.linspace(-0.5,0.5,n_filaments_per_bundle)[None,:].to(self.myosin.xs.device)
        actin_y_coord = actin_y_coord.flatten()
        self.actin.xs = torch.stack(torch.meshgrid(actin_x_coord, actin_y_coord), axis=-1).reshape(-1,2)
        self.actin.thetas = torch.zeros(self.actin.xs.shape[-2]).to(self.myosin.xs.device)
        alpha_x_coord = x_coord + unit_length/2
        self.alpha_actinin.xs = torch.stack(torch.meshgrid(alpha_x_coord, actin_y_coord), axis=-1).reshape(-1,2)
        self.alpha_actinin.n_crosslinks = self.alpha_actinin.xs.shape[-2]


    def compute_energy(self):
        return self.myosin_self_energy()  + self.actin_myosin_energy() + self.actin_alpha_energy()

    def myosin_self_energy(self):
        xs = self.myosin.xs
        radius = self.myosin.radius
        distances = pairwise_distances(xs,box=self.box,remove_diag=True)
        return 100 * (distances<2*radius).sum()
    
    def actin_myosin_energy(self):
        segments = self.actin.filament_length*torch.stack([torch.cos(self.actin.thetas), torch.sin(self.actin.thetas)], axis=-1)
        actin_endpts = torch.cat([self.actin.xs + 0.5*segments, self.actin.xs - 0.5*segments], axis=0)
        distances = pairwise_distances(self.myosin.xs,actin_endpts,box=self.box)
        n_bonds = (distances <= self.myosin.radius).sum()
        energy = n_bonds*self.e_am
        return energy
    
    def actin_alpha_energy(self):
        la = self.actin.filament_length
        ra = self.alpha_actinin.radius
        thetas = self.actin.thetas
        actin_endpts = self.actin.xs + 0.5*la*torch.stack([torch.cos(thetas), torch.sin(thetas)], axis=-1)
        x_alpha = pairwise_vectors(actin_endpts,self.alpha_actinin.xs,box=self.box)
        x_actin = -la * torch.stack([torch.cos(thetas), torch.sin(thetas)],axis=-1)
        x_actin = x_actin.unsqueeze(-2).expand_as(x_alpha)
        dx = x_alpha - x_actin

        def dot_along_last_dim(x,y):
            return torch.einsum("...i,...i->...",x,y)
        
        #check if crosslink is between the two ends of the actin filament
        mask = (dot_along_last_dim(-x_actin, dx) >= 0) & (dot_along_last_dim(x_alpha, x_actin) >= 0)
        x_alpha_norm = torch.linalg.norm(x_alpha,dim=-1)
        distances = x_alpha_norm* torch.sqrt(1 - (dot_along_last_dim(x_actin, x_alpha)/(x_alpha_norm*la))**2)
        mask = mask & ((distances <= ra)+(x_alpha_norm == 0))
        energy = mask.sum()*self.e_al
        return energy

    def actin_myosin_force(self):
        la = self.actin.filament_length
        radius = self.myosin.radius
        thetas = self.actin.thetas
        actin_endpts = self.actin.xs + 0.5*la*torch.stack([torch.cos(thetas), torch.sin(thetas)], axis=-1)
        x_myosin = pairwise_vectors(actin_endpts,self.myosin.xs,box=self.box)
        x_actin = -la * torch.stack([torch.cos(thetas), torch.sin(thetas)],axis=-1)
        x_actin = x_actin.unsqueeze(-2).expand_as(x_myosin)
        dx = x_myosin - x_actin

        def dot_along_last_dim(x,y):
            return torch.einsum("...i,...i->...",x,y)
        
        #check if crosslink is between the two ends of the actin filament
        mask = (dot_along_last_dim(-x_actin, dx) >= 0) & (dot_along_last_dim(x_myosin, x_actin) >= 0)
        x_myosin_norm = torch.linalg.norm(x_myosin,dim=-1)
        distances = x_myosin_norm* torch.sqrt(1 - (dot_along_last_dim(x_actin, x_myosin)/(x_myosin_norm*la))**2)
        mask = mask & ((distances <= radius)+(x_myosin_norm == 0))
        orientation = torch.stack([torch.cos(thetas), torch.sin(thetas)], axis=-1)
        orientation = orientation.unsqueeze(-2).expand_as(x_myosin)
        force = self.f_myosin * mask.float().unsqueeze(-1) * orientation
        force = force.sum(dim=-3)
        return force

    def myosin_self_force(self):
        xs = self.myosin.xs
        radius = self.myosin.radius
        pair_vec = pairwise_vectors(xs,box=self.box)
        distances = torch.linalg.norm(pair_vec,dim=-1)
        mask = distances < 2*radius
        force = self.f_myosin * mask.unsqueeze(-1) * pair_vec
        force = force.sum(dim=-2)
        return force
    
    def run_simulation(self, dt, D, n_steps, frame_dir="frames", make_animation=False):
        import os
        if not os.path.exists(frame_dir):
            os.mkdir(frame_dir)
        acc = 0
        traj = {"myosin": [],"actin_angle":[], "actin": [], "alpha_actinin": []}
        for i in range(n_steps):
            acc += self.mala_step(dt, D)
            if i % 300 == 0:
                # save frame to trajectory
                traj["myosin"].append(self.myosin.xs)
                traj["actin_angle"].append(self.actin.thetas)
                traj["actin"].append(self.actin.xs)
                traj["alpha_actinin"].append(self.alpha_actinin.xs)
                fig, ax = self.plot_system()
                plt.savefig("{}/frame_{:04d}.png".format(frame_dir,i))
                plt.close(fig)
                acc_rate = acc/300
                print("Acceptance rate: ", acc_rate)
                print("Step size: ", dt)
                print("Energy: ", self.energy)
                acc = 0
                if acc_rate < 0.2:
                    dt = dt*acc_rate/0.2
                elif acc_rate > 0.6:
                    dt = dt*acc_rate/0.6

        print("Energy: ")
        print(self.energy)
        # save frames as an animation
        if make_animation:
            os.system("ffmpeg -r 1 -i {}/frame_%04d.png -vcodec mpeg4 -y movie.mp4".format(frame_dir))
        for key in traj:
            traj[key] = torch.stack(traj[key], axis=0).cpu().numpy()
        np.save("trajectory.npy", traj)
    # TODO: use stokes einstein equation to compute diffusion coefficient for each component
    def mala_step(self, dt, D):
        acc = 0
        device = self.actin.xs.device
        D = torch.tensor(D).to(device)
        # update actin filaments
        xs_old = self.actin.xs.clone()
        thetas_old = self.actin.thetas.clone()
        noise = torch.sqrt(2*D*dt)*torch.randn(self.actin.n_filaments, 2).to(device)
        self.actin.xs = self.pbc(self.actin.xs + noise)
        self.actin.thetas += torch.sqrt(2*D*dt)*torch.randn(self.actin.n_filaments).to(device)
        delta_energy = self.compute_energy() - self.energy
        rand = torch.rand_like(delta_energy).to(device)
        if rand < torch.exp(-delta_energy):
            self.energy+=delta_energy
            acc+=1
        else:
            self.actin.xs = xs_old
            self.actin.thetas = thetas_old

        # update myosin bundles
        xs_old = self.myosin.xs.clone()
        #force = self.actin_myosin_force()+self.myosin_self_force()
        force = self.myosin_self_force()
        self.myosin.xs = self.pbc(self.myosin.xs + force*dt
                                  + torch.sqrt(2*D*dt)*torch.randn(self.myosin.n_bundles, 2).to(device))
        delta_energy = self.compute_energy() - self.energy
        rand = torch.rand_like(delta_energy).to(device)
        if rand < torch.exp(-delta_energy):
            self.energy+=delta_energy
            acc+=1
        else:
            self.myosin.xs = xs_old

        # update alpha-actinin
        xs_old = self.alpha_actinin.xs.clone()
        self.alpha_actinin.xs = self.pbc(self.alpha_actinin.xs + torch.sqrt(2*D*dt)*torch.randn(
            self.alpha_actinin.n_crosslinks, 2).to(device))
        delta_energy = self.compute_energy() - self.energy
        rand = torch.rand_like(delta_energy).to(device)
        if rand < torch.exp(-delta_energy):
            self.energy+=delta_energy
            acc+=1
        else:
            self.alpha_actinin.xs = xs_old
        
        return acc/3

    def pbc(self, xs):
        xs[:,0] = xs[:,0] + self.Lx * (xs[:,0] < 0) - self.Lx * (xs[:,0] >= self.Lx)
        xs[:,1] = xs[:,1] + self.Ly * (xs[:,1] < 0) - self.Ly * (xs[:,1] >= self.Ly)
        return xs

    def plot_system(self):
        n_filaments = self.actin.n_filaments
        xs = self.actin.xs.cpu()
        thetas = self.actin.thetas.cpu()
        l = self.actin.filament_length
        fig, ax = plt.subplots()
        for i in range(n_filaments):
            x = xs[i,0]
            y = xs[i,1]
            theta = thetas[i]
            ax.plot([x- 0.5*l*torch.cos(theta), x + 0.5*l*torch.cos(theta)], [y- 0.5*l*torch.sin(theta), y + 0.5*l*torch.sin(theta)], lw=1.0, color="C0")

        # plot myosin bundles
        n_bundles = self.myosin.n_bundles
        xs = self.myosin.xs.cpu()
        #force = self.actin_myosin_force()+self.myosin_self_force()
        force = self.myosin_self_force()
        force = force.cpu().numpy()
        for i in range(n_bundles):
            x = xs[i,0]
            y = xs[i,1]
            circle = plt.Circle((x, y), self.myosin.radius, color="C1", alpha=0.5)
            ax.add_artist(circle)
            ax.arrow(x, y, force[i,0], force[i,1], head_width=0.1, head_length=0.1, fc='k', ec='k')

        # plot alpha-actinin
        n_crosslinks = self.alpha_actinin.n_crosslinks
        xs = self.alpha_actinin.xs.cpu()
        for i in range(n_crosslinks):
            x = xs[i,0]
            y = xs[i,1]
            circle = plt.Circle((x, y), self.alpha_actinin.radius, color="C2", alpha=0.5)
            ax.add_artist(circle)

        ax.set_xlim(0, self.Lx)
        ax.set_ylim(0, self.Ly)
        ax.set_aspect("equal")
        return fig, ax
    