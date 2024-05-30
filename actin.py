import torch
import matplotlib.pyplot as plt
import numpy as np
from utils import plot_filaments,plot_myosin, point_point_distance, pairwise_vectors, segment_segment_distance,  point_segment_distance, dot_along_last_dim



class Filament():
    def __init__(self, n, Lx, Ly,length=3.0, device="cuda"):
        self.n = n
        self.Lx = Lx
        self.Ly = Ly
        self.length = length
        rand = torch.rand((n,3)).to(device)
        self.xs = torch.stack([Lx*rand[...,0], Ly*rand[...,1]], axis=-1)
        self.thetas = 2*torch.pi*rand[...,2]
    
    def get_endpoints(self):
        segments = self.length*torch.stack([torch.cos(self.thetas), torch.sin(self.thetas)], axis=-1)
        endpts = torch.stack([self.xs + 0.5*segments, self.xs - 0.5*segments], axis=0)
        return endpts


class MyosinBundle(Filament):
    def __init__(self, n, Lx, Ly, length=3.0, radius=1.0, device="cuda"):
        super().__init__(n, Lx, Ly, length=length, device=device)
        self.radius = radius

class AlphaActinin():
    def __init__(self, n, Lx, Ly, radius=0.2, device="cuda"):
        self.n = n
        self.Lx = Lx  
        self.Ly = Ly
        self.xs = torch.zeros((n, 2))
        rand = torch.rand((n,2)).to(device)
        self.xs = torch.stack([Lx*rand[...,0], Ly*rand[...,1]], axis=-1)
        self.radius = radius


class SarcomereModel():
    def __init__(self, n_filaments, n_bundles, n_crosslinks, Lx, Ly,
                 e_am = -1.0, e_al = -1.0, f_myosin = 20.0,
                 len_actin = 3.0, len_myosin=3.0, r_myosin = 0.5, r_alpha_actinin = 0.2,
                 catch_bond = True,
                  device="cuda"):
        self.actin = Filament(n_filaments, Lx, Ly, length=len_actin,device=device)
        #self.actin.thetas = torch.zeros_like(self.actin.thetas)
        #self.actin.thetas[1::2] = torch.pi
        self.myosin = MyosinBundle(n_bundles, Lx, Ly,length=len_myosin,radius=r_myosin, device=device)
        self.alpha_actinin = AlphaActinin(n_crosslinks, Lx, Ly,radius=r_alpha_actinin, device=device)
        self.sarcomeric_structure()
        self.Lx = Lx
        self.Ly = Ly
        self.catch_bond = catch_bond
        self.box = torch.tensor([Lx, Ly]).to(device)
        self.e_am = e_am
        self.e_al = e_al
        self.f_myosin = f_myosin
        self.energy = self.compute_energy()
        self.trajectory = None
        fig, ax = self.plot_system()
        plt.savefig("initial.png")
        plt.close(fig)
        print("initial energy: ", self.energy)
        print("number of actin filaments: ", self.actin.n)
        print("number of myosin bundles: ", self.myosin.n)
        print("number of alpha-actinin crosslinks: ", self.alpha_actinin.n)

    def sarcomeric_structure(self):
        x_coord = torch.tensor([3.,9.,15.]).to(self.myosin.xs.device)
        y_coord = torch.tensor([5.]).to(self.myosin.xs.device)
        x,y = torch.meshgrid(x_coord,y_coord)
        self.myosin.xs = torch.stack([x.flatten(),y.flatten()],axis=-1)
        self.myosin.thetas = torch.zeros_like(self.myosin.thetas)

    def load_traj(self,traj):
        print("loading trajectory")
        for key in traj.keys():
            traj[key] = list(torch.tensor(traj[key]).to(self.myosin.xs.device))
        self.myosin.xs = traj["myosin"][-1]
        self.actin.xs = traj["actin"][-1]
        self.actin.thetas = traj["actin_angle"][-1]
        self.alpha_actinin.xs = traj["alpha_actinin"][-1]
        self.energy = self.compute_energy()
        self.trajectory = traj
        fig, ax = self.plot_system()
        plt.savefig("initial.png")
        plt.close(fig)


    def compute_energy(self):
        return self.actin_myosin_energy() + self.repulsive_energy()  + self.actin_alpha_energy()

    def repulsive_energy(self):
        energy = 0
        #myosin-myosin repulsion
        myosin_endpts = self.myosin.get_endpoints()
        distances = segment_segment_distance(myosin_endpts[0],myosin_endpts[1],myosin_endpts[0],myosin_endpts[1],box=self.box,remove_diag=True)
        energy = 100 * (distances<2*self.myosin.radius).sum()
        #print("myosin-myosin repulsion: ", energy)
        #aa-aa repulsion
        distances = point_point_distance(self.alpha_actinin.xs,box=self.box,remove_diag=True)
        energy += 10 * (distances<2*self.alpha_actinin.radius).sum()
        #aa-myosin repulsion
        distances = point_segment_distance(self.alpha_actinin.xs,myosin_endpts[0],myosin_endpts[1],box=self.box)
        energy += 100 * (distances<self.myosin.radius+self.alpha_actinin.radius).sum()
        #print(energy)
        return energy
    
    def actin_myosin_energy(self):
        # segments = self.actin.length*torch.stack([torch.cos(self.actin.thetas), torch.sin(self.actin.thetas)], axis=-1)
        # actin_plus = self.actin.xs + 0.5*segments
       
        # distances = point_segment_distance(actin_plus,myosin_endpts[0],myosin_endpts[1],box=self.box)
        # #n_bonds = (distances <= self.myosin.radius).sum()
        # angles = self.actin.thetas.unsqueeze(-1)- self.myosin.thetas.unsqueeze(-2)
        # angles = angles - 2*torch.pi*torch.round(angles/(2*torch.pi))
        # strength = torch.abs(torch.cos(angles))-0.707
        # strength = strength * (strength>0)
        # energy = self.e_am*(strength * (distances <= self.myosin.radius)).sum()
        myosin_endpts = self.myosin.get_endpoints()
        actin_endpts = self.actin.get_endpoints()
        distances = segment_segment_distance(myosin_endpts[0],myosin_endpts[1],actin_endpts[0],actin_endpts[1],
                                          box=self.box)
        n_bonds = (distances<self.myosin.radius).sum()
        energy = n_bonds*self.e_am
        return energy
    
    def actin_alpha_energy(self):
        ra = self.alpha_actinin.radius
        actin_endpts = self.actin.get_endpoints()
        distances = point_segment_distance(self.alpha_actinin.xs,actin_endpts[0],actin_endpts[1],box=self.box)
        mask = (distances <= ra)
        #energy = mask.sum()*self.e_al
        energy = 0
        if self.catch_bond:
            self.aa_strength = self.directional_catch_bond(mask)
            energy += 2*self.e_al*self.aa_strength.sum()
        else:
            energy += self.e_al*mask.sum()
        return energy
    
    def directional_catch_bond(self, aa_actin_connectivity):
        #size of aa_actin_connectivity: (n_aa, n_actin)
        force_on_actin = self.actin_myosin_force()[0]
        force_mag = torch.norm(force_on_actin,dim=-1)
        norm_factor = force_mag.unsqueeze(-1)*force_mag.unsqueeze(-2)
        pairwise_cosine = dot_along_last_dim(force_on_actin.unsqueeze(-2),force_on_actin.unsqueeze(-3))/(norm_factor+1e-5)
        strength = norm_factor*torch.abs(pairwise_cosine.clamp(max=0))/self.f_myosin**2#size: (n_actin, n_actin)
        connectivity_expanded = aa_actin_connectivity.unsqueeze(-1)*aa_actin_connectivity.unsqueeze(-2)
        strength_expanded = strength.unsqueeze(-3)
        catch_bond_strength = (connectivity_expanded*strength_expanded).sum(dim=[-2,-1])
        catch_bond_strength = torch.clamp(catch_bond_strength,max=3)
        return catch_bond_strength
    
    def actin_myosin_force(self):
        actin_endpts = self.actin.get_endpoints()
        myosin_endpts = self.myosin.get_endpoints()
        distances = segment_segment_distance(actin_endpts[0],actin_endpts[1],myosin_endpts[0],myosin_endpts[1],
                                           box=self.box)
        mask = (distances <= self.myosin.radius)
        orientation = torch.stack([torch.cos(self.actin.thetas), torch.sin(self.actin.thetas)], axis=-1)
        orientation = orientation.unsqueeze(-2)
        force = self.f_myosin * (mask.float()).unsqueeze(-1) * orientation
        force_on_actin = force.sum(dim=-2)   
        self.force_on_actin = force_on_actin
        force_on_myosin = -force.sum(dim=-3)
        return force_on_actin, force_on_myosin
    
    def run_simulation(self, dt, D, n_steps, kT=1, frame_dir="frames", save_every=200,update_dt_every=200, update_myosin_every=100, make_animation=True):
        D = torch.tensor(D).to(self.myosin.xs.device)
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
            update_myosin = i % update_myosin_every == 0
            acc += self.mc_step(dt, D, kT, update_myosin=update_myosin)
            # for _ in range(200):
            #     acc += self.mc_step(dt, D, kT, update_myosin=update_myosin)
            #     print(acc)
            # exit()
            if i % save_every == 0:
                # save frame to trajectory
                traj["myosin"].append(self.myosin.xs)
                traj["actin_angle"].append(self.actin.thetas)
                traj["actin"].append(self.actin.xs)
                traj["alpha_actinin"].append(self.alpha_actinin.xs)
                fig, ax = self.plot_system()
                plt.savefig(file_format.format(frame_dir,i//save_every))
                plt.close(fig)
                traj_np = traj.copy()
                for key in traj:
                    traj_np[key] = torch.stack(traj[key], axis=0).cpu().numpy()
                print("Saving trajectory...")
                np.save(frame_dir+"/trajectory.npy", traj_np)
            if (i-start)>0 and update_dt_every>0 and (i-start)% update_dt_every == 0:
                acc_rate = acc/update_dt_every
                print(acc_rate)
                print("Step: ", i)
                print("Acceptance rate: ", acc_rate)
                print("Step size: ", dt)
                print("Energy: ", self.energy)
                acc = 0
                if acc_rate < 0.2:
                    dt = dt*max(acc_rate/0.2,0.1)
                elif acc_rate > 0.6:
                    dt = min(dt*acc_rate/0.6,0.01)
        for key in traj:
            traj[key] = torch.stack(traj[key], axis=0).cpu().numpy()
        print("Saving trajectory...")
        np.save(frame_dir+"/trajectory.npy", traj)
        # save frames as an animation
        if make_animation:
            os.system("ffmpeg -r 10 -i {}/frame_{%0".format(frame_dir)+str(n_digits)+"d}.png -vcodec mpeg4 -y movie.mp4")
        
    def mc_step(self, dt, D, kT, update_myosin=False):
        acc=0
        device = self.myosin.xs.device
        #update actin filaments
        index = torch.randint(0,self.actin.n,(1,)).item()
        xs_old = self.actin.xs[index].clone()
        thetas_old = self.actin.thetas[index].clone()
        force_on_actin, force_on_myosin = self.actin_myosin_force()
        noise = torch.sqrt(2*D*dt)*torch.randn(3).to(device)
        self.actin.xs[index] = self.pbc(self.actin.xs[index] + noise[:2])
        self.actin.thetas[index] += noise[-1]
        self.actin.thetas[index] = torch.fmod(self.actin.thetas[index],2*np.pi)

        force_on_actin_new, force_on_myosin = self.actin_myosin_force()
        force = 0.5*(force_on_actin[index] + force_on_actin_new[index])
        #print("computing actin energy")
        delta_energy = self.compute_energy() - self.energy
        work = torch.dot(force,noise[:2])
        rand = torch.rand_like(delta_energy).to(device)
        #print("a:",delta_energy,work)
        if rand < torch.exp((-delta_energy+work)/kT):
            acc+=1
            self.energy += delta_energy
        else:
            self.actin.xs[index] = xs_old
            self.actin.thetas[index] = thetas_old
        #update alpha-actinin
        index = torch.randint(0,self.alpha_actinin.n,(1,)).item()
        xs_old = self.alpha_actinin.xs[index].clone()
        noise = torch.sqrt(2*D*dt)*torch.randn(2).to(device)
        self.alpha_actinin.xs[index] = self.pbc(self.alpha_actinin.xs[index] + noise)
        delta_energy = self.compute_energy() - self.energy
        rand = torch.rand_like(delta_energy).to(device)
        if rand < torch.exp(-delta_energy/kT):
            acc+=1
            self.energy += delta_energy
        else:
            self.alpha_actinin.xs[index] = xs_old
        #update myosin
        D_myosin = D*0.2
        index = torch.randint(0,self.myosin.n,(1,)).item()
        xs_old = self.myosin.xs[index].clone()
        noise = torch.sqrt(2*D_myosin*dt)*torch.randn(3).to(device)
        self.myosin.xs[index] = self.pbc(self.myosin.xs[index] + noise[:2])
        self.myosin.thetas[index] += noise[-1]
        self.myosin.thetas[index] = torch.fmod(self.myosin.thetas[index],2*np.pi)
        force_on_actin, force_on_myosin_new = self.actin_myosin_force()
        force = 0.5*(force_on_myosin[index] + force_on_myosin_new[index])
        work = torch.dot(force,noise[:2])
        delta_energy = self.compute_energy() - self.energy
        rand = torch.rand_like(delta_energy).to(device)
        if rand < torch.exp((-delta_energy+work)/kT):
            acc+=1
            self.energy += delta_energy
        else:
            self.myosin.xs[index] = xs_old
        #print("Energy: ", self.compute_energy())
        # if update_myosin:
        #     self.myosin.xs = self.pbc(self.myosin.xs + force_on_myosin*dt)
        return acc/3

    def pbc(self, xs):
        xs[...,0] = xs[...,0] + self.Lx * (xs[...,0] < 0) - self.Lx * (xs[...,0] >= self.Lx)
        xs[...,1] = xs[...,1] + self.Ly * (xs[...,1] < 0) - self.Ly * (xs[...,1] >= self.Ly)
        return xs

    def plot_system(self):

        fig, ax = plt.subplots(dpi=500)
        plot_filaments(self.actin,ax, Lx=self.Lx, Ly=self.Ly, color='C0')
        plot_myosin(self.myosin,ax, Lx=self.Lx, Ly=self.Ly, radius=self.myosin.radius, color='C1')
      

        # plot alpha-actinin
        n_crosslinks = self.alpha_actinin.n
        xs = self.alpha_actinin.xs.cpu()
        if self.catch_bond:
            strength = torch.abs(self.aa_strength)
        for i in range(n_crosslinks):
            x = xs[i,0]
            y = xs[i,1]    
            if self.catch_bond:
                binding_affinity = float(strength[i])/2
                binding_affinity = min(1.,float(binding_affinity))
                circle = plt.Circle((x, y), self.alpha_actinin.radius, facecolor=(1-binding_affinity,1-binding_affinity,1-binding_affinity),
                                   alpha=0.5,linewidth=0.7,edgecolor="green")
            else:
                circle = plt.Circle((x, y), self.alpha_actinin.radius, color="C2", alpha=0.5)
            ax.add_artist(circle)

        ax.set_xlim(0, self.Lx)
        ax.set_ylim(0, self.Ly)
        ax.set_aspect("equal")

        return fig, ax
    