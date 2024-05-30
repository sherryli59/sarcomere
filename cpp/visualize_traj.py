import h5py
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
import math
import joblib
from joblib import Parallel, delayed
import os

def plot_filaments(xs,thetas,l,Lx,Ly,ax,color='k'):
    for i in range(xs.shape[0]):
        x = xs[i,0]
        y = xs[i,1]
        theta = thetas[i,0]
        delta_x = l*np.cos(theta)
        delta_y = l*np.sin(theta)
        ax.arrow(x-delta_x, y-delta_y, delta_x, delta_y,
                    head_width=0.1, head_length=0.1, fc=color, ec=color)
        #plot the periodic images
        if x+delta_x>Lx/2:
            ax.arrow(x-delta_x-Lx, y-delta_y, delta_x, delta_y,
                        head_width=0.1, head_length=0.1, fc=color, ec=color)
        elif x-delta_x<-Lx/2:
            ax.arrow(x-delta_x+Lx, y-delta_y, delta_x, delta_y,
                        head_width=0.1, head_length=0.1, fc=color, ec=color)
        if y+delta_y>Ly/2:
            ax.arrow(x-delta_x, y-delta_y-Ly, delta_x, delta_y,
                        head_width=0.1, head_length=0.1, fc=color, ec=color)
        elif y-delta_y<-Ly/2:
            ax.arrow(x-delta_x, y-delta_y+Ly, delta_x, delta_y,
                        head_width=0.1, head_length=0.1, fc=color, ec=color)

def plot_myosin(xs,thetas,l,Lx,Ly,radius,ax,color='k'):
    for i in range(xs.shape[0]):
        x = xs[i,0]
        y = xs[i,1]
        theta = thetas[i,0]
        delta_x = l/2 * np.cos(theta) - radius * np.sin(theta)
        delta_y = l/2 * np.sin(theta) + radius * np.cos(theta)
        #draw a dot at x,y
        ax.plot(x,y,'o',color='k')
        #plot a rectangle to represent the myosin
        ax.add_patch(Rectangle((x-delta_x,y-delta_y),width=l,height=2*radius,angle=theta*180/math.pi,fill=True,facecolor=color,alpha=0.5))
        #plot the periodic images
        if x+delta_x>Lx:
            ax.add_patch(Rectangle((x-delta_x-Lx,y-delta_y),width=l,height=2*radius,angle=theta*180/math.pi,fill=True,facecolor=color,alpha=0.5))
        elif x-delta_x<0:
            ax.add_patch(Rectangle((x-delta_x+Lx,y-delta_y),width=l,height=2*radius,angle=theta*180/math.pi,fill=True,facecolor=color,alpha=0.5))
        if y+delta_y>Ly:
            ax.add_patch(Rectangle((x-delta_x,y-delta_y-Ly),width=l,height=2*radius,angle=theta*180/math.pi,fill=True,facecolor=color,alpha=0.5))
        elif y-delta_y<0:
            ax.add_patch(Rectangle((x-delta_x,y-delta_y+Ly),width=l,height=2*radius,angle=theta*180/math.pi,fill=True,facecolor=color,alpha=0.5))

def plot_system(frame,data,myosin_radius,alpha_actinin_radius,actin_length,myosin_length, Lx,Ly):

        fig, ax = plt.subplots(dpi=500)
        plot_filaments(data["/actin/xs"][frame], data["/actin/thetas"][frame],actin_length,Lx=Lx, Ly=Ly,ax=ax, color='C0')
        plot_myosin(data["/myosin/xs"][frame], data["/myosin/thetas"][frame],myosin_length,Lx=Lx, Ly=Ly,radius=myosin_radius,ax=ax,color='C1')
      
        # plot alpha-actinin
        xs = data["/alpha_actinin/xs"][frame]
        for i in range(xs.shape[0]):
            x = xs[i,0]
            y = xs[i,1]    

            circle = plt.Circle((x, y), alpha_actinin_radius, color="C2", alpha=0.5)
            ax.add_artist(circle)

        ax.set_xlim(-Lx/2, Lx/2)
        ax.set_ylim(-Ly/2, Ly/2)
        ax.set_aspect("equal")
        return fig,ax

def plot(ind, nframes, nworkers,**kwargs):
    frame_start = int(ind*(nframes-1)/nworkers)
    frame_end = int((ind+1)*(nframes-1)/nworkers)     
    for frame in range(frame_start,frame_end):
        fig, ax = plot_system(frame=frame,**kwargs)
        plt.savefig(file_format.format(frame_dir,frame))
        plt.close(fig)


def hdf5_to_dict(hdf5_file):
    data_dict = {}
    def _hdf5_to_dict_recursive(group, group_name):
        for key in group.keys():
            if isinstance(group[key], h5py.Dataset):
                data_dict[f"{group_name}/{key}"] = np.array(group[key])
            elif isinstance(group[key], h5py.Group):
                _hdf5_to_dict_recursive(group[key], f"{group_name}/{key}")
    
    _hdf5_to_dict_recursive(hdf5_file, "")
    return data_dict

if __name__ == "__main__":
    filename = 'build/traj.h5'
    traj = h5py.File(filename, 'r')
    data = hdf5_to_dict(traj)
    energy = data["/energy/total_energy"]
    traj = h5py.File(filename, 'r')
    data = hdf5_to_dict(traj)
    nframes = data["/actin/xs"].shape[0]
    actin_length = 3
    myosin_length = 3
    Lx = 10
    Ly = 10
    myosin_radius = 0.5
    alpha_actinin_radius = 0.2
    cpu_workers = joblib.cpu_count()

    n_digits = len(str(nframes))
    file_format = "{}/frame_{:0"+str(n_digits)+"d}.png"
    frame_dir = "frames"
    if not os.path.exists(frame_dir):
        os.mkdir(frame_dir)
    Parallel(n_jobs=cpu_workers)(delayed(plot)(i, nframes, cpu_workers,
        data=data, myosin_radius=myosin_radius, alpha_actinin_radius=alpha_actinin_radius,
        actin_length=actin_length, myosin_length=myosin_length, Lx=Lx, Ly=Ly) for i in range(cpu_workers))