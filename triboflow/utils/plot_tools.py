#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct  7 15:46:24 2020

Utility tools to plot the High Simmetry (HS) points for slab and interface

The module contains the following functions:

    - plot_slab_hs
    - plot_pes
    - plot_uniform_grid

    Author: Gabriele Losi (glosi000)
    Copyright 2021, Prof. M.C. Righi, TribChem, ERC-SLIDE, University of Bologna

"""

__author__ = "Gabriele Losi"
__copyright__ = "Copyright 2021, Prof. M.C. Righi, TribChem, ERC-SLIDE, University of Bologna"
__contact__ = "clelia.righi@unibo.it"
__date__ = "February 8th, 2021"

import matplotlib.pyplot as plt
import numpy as np
from pymatgen.analysis.adsorption import plot_slab

from triboflow.phys.high_symmetry import hs_dict_converter


# =============================================================================
# PLOTTING TOOLS
# =============================================================================


def export_legend(legend, filename="legend.png", expand=[-5, -5, 5, 5]):
    fig = legend.figure
    fig.canvas.draw()
    bbox = legend.get_window_extent()
    bbox = bbox.from_extents(*(bbox.extents + np.array(expand)))
    bbox = bbox.transformed(fig.dpi_scale_trans.inverted())
    fig.savefig(filename, dpi="figure", bbox_inches=bbox)


def plot_slab_hs(
    hs, slab, to_fig=None, hs_type="all", leg_size=10, in_frac=False
):
    """
    Plot the slab, displaying the atoms and the HS sites of the surface

    Parameters
    ----------
    slab : pymatgen.core.surface.Slab
        The slab to be displayed

    hs : dict
        HS sites of the slab.

    name_fig : string, optional
        Name of the image that you want to save, it will be: 'name_fig'+'.pdf'
        Suggested name: 'Name of the material' + 'Miller index'.
        The default is None and no image is saved.

    Returns
    -------
    None.

    """

    # Check the type of the hs points
    typ = list(set(type(k) for k in hs.values()))[0]
    if typ == list:
        hs = hs_dict_converter(hs, to_array=True)

    # Extract the lattice vector of the basal plane
    a = slab.lattice.matrix[0, :]
    b = slab.lattice.matrix[1, :]

    # plot the atoms and the lattice cell
    fig = plt.figure()
    ax = fig.add_subplot(111)

    plot_slab(
        slab,
        ax,
        scale=0.8,
        repeat=3,
        window=2.25,
        draw_unit_cell=True,
        decay=0.2,
        adsorption_sites=False,
    )
    xlower = min(slab.lattice.matrix[:, 0])
    xupper = max([a[0] + b[0], a[0], b[0]])
    ylower = min(slab.lattice.matrix[:, 1])
    yupper = max([a[1] + b[1], a[1], b[1]])
    ax.set(
        xlim=(xlower - 0.5, xupper + 0.5), ylim=(ylower - 0.5, yupper + 0.5)
    )

    xrange = abs(xupper - xlower)
    yrange = abs(yupper - ylower)

    ratio = min(xrange / yrange, yrange / xrange)

    # Add the HS sites with the proper labels
    for k in hs.keys():
        data = hs[k]
        if in_frac:
            data = np.dot(data, slab.lattice.matrix[:2, :2])
        if len(data.shape) == 1:
            plt.plot(
                data[0],
                data[1],
                marker="o",
                markersize=12 * (ratio),
                mew=0.5,
                linestyle="",
                zorder=10000,
                label=k,
                markeredgecolor="black",
            )
        else:
            if hs_type == "all":
                plt.plot(
                    data[:, 0],
                    data[:, 1],
                    marker="o",
                    markersize=18 * (ratio),
                    mew=0.5,
                    linestyle="",
                    zorder=10000,
                    label=k,
                    mec="black",
                )
            elif hs_type == "unique":
                plt.plot(
                    data[:, 0],
                    data[:, 1],
                    marker="o",
                    markersize=12,
                    mew=0.5,
                    linestyle="",
                    zorder=10000,
                    label=k,
                    mec="black",
                )

    plt.legend(
        bbox_to_anchor=(1.025, 1), loc="upper left", prop={"size": leg_size}
    )
    # legend = plt.legend(bbox_to_anchor=(1.025, 1), loc='upper left', prop={'size': leg_size})
    # export_legend(legend, f'{to_fig}_legend.png')
    # legend.remove()

    plt.rcParams.update({"font.size": 18})
    plt.tick_params(
        axis="both",  # changes apply to the x-axis
        which="both",  # both major and minor ticks are affected
        bottom=False,  # ticks along the bottom edge are off
        top=False,  # ticks along the top edge are off
        left=False,
        right=False,
        labelbottom=False,
        labelleft=False,
    )  # labels along the bottom edge are off

    if to_fig != None:
        plt.savefig(to_fig + ".png", dpi=300, bbox_inches="tight")

    plt.show()


def plot_pes(data, lattice, to_fig=None, vmin=None, vmax=None, plot_hs=None):
    """
    Plot the PES and eventually save it

    """

    a = lattice[0]
    b = lattice[1]
    x = data[0]
    y = data[1]
    E = data[2]

    n = 51
    if vmin and vmax:
        levels = np.linspace(vmin, vmax, n)
    else:
        levels = np.linspace(np.amin(E), np.amax(E), n)
    fig = plt.figure(figsize=(7, 7), dpi=150)
    ax = fig.add_subplot(111)
    ax.set_aspect("equal")
    anglerot = "vertical"
    # zt1=plt.contourf(x, y, E, level, extent=(-fact*a, fact*a, -fact*b, fact*b), cmap=plt.cm.RdYlBu_r)
    zt1 = plt.contourf(x, y, E, levels, cmap=plt.cm.RdYlBu_r)
    xrange, yrange = sum(abs(lattice[:2, 0])), sum(abs(lattice[:2, 1]))
    ratio = xrange / yrange
    xlim, ylim = ax.get_xlim(), ax.get_ylim()
    if ratio < 1:
        ratio = round(ratio**-1)
        ax.set_ylim((ylim[0], ylim[1] / ratio))
    else:
        ratio = round(ratio)
        ax.set_xlim((xlim[0], xlim[1] / ratio))

    cbar1 = plt.colorbar(
        zt1,
        ax=ax,
        orientation=anglerot,
        shrink=1.1 * min(yrange / xrange, 1.0),
    )
    cbar1.set_label(
        r"$E_{adh} (J/m^2)$",
        rotation=270,
        labelpad=20,
        fontsize=15,
        family="serif",
    )
    # ax.quiver(0. , 0., 1., 0.,scale=1.,scale_units='inches',width=0.01,color='white')
    # ax.quiver(0. , 0., 0., 1.,scale=1.,scale_units='inches',width=0.01,color='white')
    # ax.plot(0.,0.,'w.',ms=7)
    # ax.text(0.5,-0.5,'[1 0 1]',rotation='horizontal',color='white', fontsize=14)
    # ax.text(-0.5,1.,'[1 2 1]',rotation='vertical',color='white', fontsize=14)
    # ax.axis([-fact*min(a), fact*max(a), -fact*min(b), fact*max(b)])
    plt.xlabel("distance (Å)", fontsize=12, family="serif")
    plt.ylabel("distance (Å)", fontsize=12, family="serif")

    # plot the cell and the combined unique high symmetry points over the PES plot
    if plot_hs:
        import matplotlib.patches as patches

        x = [0, a[0], a[0] + b[0], b[0]]
        y = [0, 0, b[1], a[1] + b[1], b[1]]
        ax.add_patch(patches.Polygon(xy=list(zip(x, y)), fill=False, lw=2))

        for hsp, shift in plot_hs.items():
            pt = shift[0]
            plt.scatter(pt[0], pt[1], edgecolors="black", s=60)

    for zt1 in zt1.collections:
        zt1.set_edgecolor("face")
        zt1.set_linewidth(0.000000000001)

    if to_fig != None:
        plt.title("PES for " + str(to_fig), fontsize=15, family="serif")
        plt.savefig("PES_" + str(to_fig) + ".png", dpi=300)


def plot_uniform_grid(grid, cell, n_a, n_b):
    """
    Plot a uniform grid of n_aXn_b points on the planar base of a lattice

    """

    a = cell[0, :]
    b = cell[1, :]
    v = np.cross(a, b)

    mod_a = np.sqrt(a[0] ** 2.0 + a[1] ** 2.0 + a[2] ** 2.0)
    mod_b = np.sqrt(b[0] ** 2.0 + b[1] ** 2.0 + b[2] ** 2.0)
    A = np.sqrt(v[0] ** 2.0 + v[1] ** 2.0 + v[2] ** 2.0)

    N = n_a * n_b
    density = N / A

    # Print information
    print("1st vector:  {:} -> norm: {:.3f}".format(a, mod_a))
    print("2nd vector:  {:} -> norm: {:.3f}".format(b, mod_b))
    print("N pts: {:}   Area: {:.3f}   Density: {:.3f}".format(N, A, density))
    print("\nUniform {0}x{1} grid\n".format(n_a, n_b))
    print(grid)

    # Projection on the plane, top view
    plt.title("Projection on xy plane")
    plt.plot(grid[:, 0], grid[:, 1], "o")

    # 3D plot
    fig = plt.figure()
    ax = fig.add_subplot(111, projection="3d")
    ax.scatter(
        grid[:, 0], grid[:, 1], grid[:, 2], c="r", marker="o", label="3D grid"
    )

    # Plot the lattice edge of the plane
    x = [0, a[0], a[0] + b[0], b[0], 0]
    y = [0, a[1], a[1] + b[1], b[1], 0]
    z = [0, a[2], a[2] + b[2], b[2], 0]

    ax.plot(x, y, z)
    plt.show()


def surfen_graph(surfen_list, conv_thr=0.01, plot_title=None, to_fig=None):
    import numpy as np
    from matplotlib import pyplot as plt

    data = np.asarray(surfen_list)
    d = data.copy()
    d[:, 1] = (d[:, 1] / d[-1, 1] - 1) * 100
    fig, ax = plt.subplots()
    ax.plot(d[:, 0], d[:, 1], "ro:")
    ax.hlines(
        [-100 * conv_thr, 100 * conv_thr],
        min(d[:, 0]),
        max(d[:, 0]),
        linestyles="dashed",
        colors="k",
    )
    if plot_title:
        ax.set_title(plot_title)
    if to_fig:
        fig.savefig(to_fig + ".png", dpi=300, bbox_inches="tight")
