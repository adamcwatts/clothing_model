import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import finite_dif_pract_1
import os
import imageio as io


def gif_creator(gif_name, list_of_plots):
    # USE IMAGE IO TO CREATE A GIF
    with io.get_writer(f'{gif_name}.gif', mode='I', duration=0.1) as writer:
        # ITERATE OVER FILENAMES
        for plot in list_of_plots:
            # READ IN FILE
            image = io.imread(plot)
            # APPEND FILE TO GIF
            writer.append_data(image)
    writer.close()


analytical_df = pd.read_csv('analytical_solution.csv')
analytical_df = analytical_df.drop(analytical_df.columns[0], axis=1)

finite_df = pd.read_csv('finite_solution.csv')
finite_df = finite_df.drop(finite_df.columns[0], axis=1)

difference = analytical_df.iloc[0:10000, :] - finite_df.iloc[0:10000, :]
dx = finite_dif_pract_1.dx
dt = finite_dif_pract_1.dt
x_array_meters = np.arange(0, 0.01, dx)
x_array_mm = x_array_meters * 1000

for i in range(100 + 1):
    fig, ax = plt.subplots(ncols=1, figsize=(8, 4))
    ax.set_xlim((0, 10))
    ax.set_ylim((0, 40))
    ax.set_xlabel('Fabric Length [mm]')
    ax.set_ylabel(r'Temperature [$^{\circ}$C]')
    plt.plot(x_array_mm, analytical_df.iloc[i * 100, :], label='Fourier Analytical Solution: 500 Terms')  # analytical
    plt.plot(x_array_mm, finite_df.iloc[i * 100, :], label='Forward Finite Difference')  # FD
    plt.title(
        f'Analytical vs Numerical Forward Finite Difference Results\nHomogeneous Dirichlet BC\nTime = {round(dt * 100 * i, 3)} Seconds',
        pad=35)
    plt.legend(bbox_to_anchor=(0, 1.02, 1, 0.2), loc="lower left",
               mode="expand", borderaxespad=0, ncol=3)
    title = './plots/' + str(i).zfill(3) + '_figure.png'
    plt.savefig(title, bbox_inches='tight')
    plt.show()

file_names = sorted((fn for fn in os.listdir('./plots')))
prefix = './plots/'
for i, file in enumerate(file_names):
    file_names[i] = prefix + file

gif_creator('numerical_scheme_comparisons', file_names)
