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

def dot_along_last_dim(x,y):
            return torch.einsum("...i,...i->...",x,y)

def point_segment_distance(x, a, b, box=None):
    # compute the distances between all points in x and all segments defined by b-a
    # x: (n_points, dim)
    # a,b : (n_segments, dim)
    ab = b - a
    ab_norm = torch.linalg.norm(ab, axis=-1)
    #ap: (n_points, n_segments, dim)
    ap = pairwise_vectors(x,a,box)
    ab = ab/ab_norm.unsqueeze(-1)
    ab = ab.unsqueeze(-3).expand(ap.shape)
    ap_ab = dot_along_last_dim(ap,ab)
    ap_ab = torch.clamp(ap_ab, torch.zeros_like(ap_ab), ab_norm.unsqueeze(-2))
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
        self.energy = self.compute_energy()
        self.trajectory = None
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

    def load_traj(self,traj):
        for key in traj.keys():
            traj[key] = list(torch.tensor(traj[key]).to(self.myosin.xs.device))
        self.myosin.xs = traj["myosin"][-1]
        self.actin.xs = traj["actin"][-1]
        self.actin.thetas = traj["actin_angle"][-1]
        self.alpha_actinin.xs = traj["alpha_actinin"][-1]
        self.energy = self.compute_energy()
        self.trajectory = traj


    def compute_energy(self):
        return self.repulsive_energy()  + self.actin_myosin_energy() + self.actin_alpha_energy()

    def repulsive_energy(self):
        #myosin-myosin repulsion
        distances = pairwise_distances(self.myosin.xs,box=self.box,remove_diag=True)
        energy = 100 * (distances<2*self.myosin.radius).sum()
        #aa-aa repulsion
        # distances = pairwise_distances(self.alpha_actinin.xs,box=self.box,remove_diag=True)
        # energy += 100 * (distances<self.alpha_actinin.radius).sum()
        #aa-myosin repulsion
        distances = pairwise_distances(self.myosin.xs,self.alpha_actinin.xs,box=self.box)
        energy += 100 * (distances<self.myosin.radius+self.alpha_actinin.radius).sum()
        return energy
    
    def actin_myosin_energy(self):
        segments = self.actin.filament_length*torch.stack([torch.cos(self.actin.thetas), torch.sin(self.actin.thetas)], axis=-1)
        #actin_endpts = torch.cat([self.actin.xs + 0.5*segments, self.actin.xs - 0.5*segments], axis=0)
        actin_endpts = self.actin.xs + 0.5*segments
        distances = pairwise_distances(self.myosin.xs,actin_endpts,box=self.box)
        n_bonds = (distances <= self.myosin.radius).sum()
        energy = n_bonds*self.e_am
        return energy
    
    def actin_alpha_energy(self):
        la = self.actin.filament_length
        ra = self.alpha_actinin.radius
        thetas = self.actin.thetas
        actin_endpts = self.actin.xs + 0.5*la*torch.stack([torch.cos(thetas), torch.sin(thetas)], axis=-1)
        actin_endpts_2 = self.actin.xs - 0.5*la*torch.stack([torch.cos(thetas), torch.sin(thetas)], axis=-1)
        distances = point_segment_distance(self.alpha_actinin.xs,actin_endpts,actin_endpts_2,box=self.box)
        mask = (distances <= ra)
        energy = mask.sum()*self.e_al
        return energy

    def actin_myosin_force(self):
        la = self.actin.filament_length
        radius = self.myosin.radius
        thetas = self.actin.thetas
        actin_endpts = self.actin.xs + 0.5*la*torch.stack([torch.cos(thetas), torch.sin(thetas)], axis=-1)
        actin_endpts_2 = self.actin.xs - 0.5*la*torch.stack([torch.cos(thetas), torch.sin(thetas)], axis=-1)
        distances = point_segment_distance(self.myosin.xs,actin_endpts,actin_endpts_2,box=self.box)
        mask = (distances <= radius)
        orientation = torch.stack([torch.cos(thetas), torch.sin(thetas)], axis=-1)
        orientation = orientation.unsqueeze(-3)
        force = self.f_myosin * mask.float().unsqueeze(-1) * orientation
        force = force.sum(dim=-2)
        return force

    def myosin_self_force(self):
        xs = self.myosin.xs
        radius = self.myosin.radius
        pair_vec = pairwise_vectors(xs,box=self.box)
        distances = torch.linalg.norm(pair_vec,dim=-1)
        distances = distances + (distances == 0)
        mask = distances < 2*radius
        force = self.f_myosin * mask.unsqueeze(-1) * pair_vec/distances.unsqueeze(-1)
        force = force.sum(dim=-2)
        return force
    
    def run_simulation(self, dt, D, n_steps, frame_dir="frames", save_every=300, make_animation=True):
        import os
        if not os.path.exists(frame_dir):
            os.mkdir(frame_dir)
        acc = 0
        if self.trajectory is None:
            traj = {"myosin": [],"actin_angle":[], "actin": [], "alpha_actinin": []}
            start = 0
        else:
            traj = self.trajectory
            start = len(traj["myosin"])*save_every
        n_saved_frames = n_steps//save_every
        #count the number of digits in n_saved_frames
        n_digits = len(str(n_saved_frames))
        file_format = "{}/frame_{:0"+str(n_digits)+"d}.png"
        for i in range(start, n_steps):
            acc += self.mala_step(dt, D)
            if i % save_every == 0:
                # save frame to trajectory
                traj["myosin"].append(self.myosin.xs)
                traj["actin_angle"].append(self.actin.thetas)
                traj["actin"].append(self.actin.xs)
                traj["alpha_actinin"].append(self.alpha_actinin.xs)
                fig, ax = self.plot_system()
                plt.savefig(file_format.format(frame_dir,i//save_every))
                plt.close(fig)
                if i>0:
                    acc_rate = acc/save_every
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
            os.system("ffmpeg -r 10 -i "+file_format+" -vcodec mpeg4 -y movie.mp4".format(frame_dir))
        for key in traj:
            traj[key] = torch.stack(traj[key], axis=0).cpu().numpy()
        np.save(frame_dir+"/trajectory.npy", traj)

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
        force = self.actin_myosin_force()+self.myosin_self_force()
        #force = self.myosin_self_force()
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
        fig, ax = plt.subplots(dpi=500)
        for i in range(n_filaments):
            x = xs[i,0]
            y = xs[i,1]
            theta = thetas[i]
            delta_x = l*torch.cos(theta)
            delta_y = l*torch.sin(theta)
            ax.arrow(x-0.5*delta_x, y-0.5*delta_y, delta_x, delta_y,
                      head_width=0.1, head_length=0.1, fc='k', ec='k')
        # plot myosin bundles
        n_bundles = self.myosin.n_bundles
        xs = self.myosin.xs.cpu()
        force = self.actin_myosin_force()+self.myosin_self_force()
        force = force.cpu().numpy()
        for i in range(n_bundles):
            x = xs[i,0]
            y = xs[i,1]
            circle = plt.Circle((x, y), self.myosin.radius, color="C1", alpha=0.5)
            ax.add_artist(circle)
            if x+self.myosin.radius>self.Lx:
                circle = plt.Circle((x-self.Lx, y), self.myosin.radius, color="C1", alpha=0.5)
                ax.add_artist(circle)
            elif x-self.myosin.radius<0:
                circle = plt.Circle((x+self.Lx, y), self.myosin.radius, color="C1", alpha=0.5)
                ax.add_artist(circle)
            if y+self.myosin.radius>self.Ly:
                circle = plt.Circle((x, y-self.Ly), self.myosin.radius, color="C1", alpha=0.5)
                ax.add_artist(circle)
            elif y-self.myosin.radius<0:
                circle = plt.Circle((x, y+self.Ly), self.myosin.radius, color="C1", alpha=0.5)
                ax.add_artist(circle)
            ax.arrow(x, y, force[i,0], force[i,1],alpha=0.4, head_width=0.2, head_length=0.2,
                      fc='C1', ec='C1')

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
    