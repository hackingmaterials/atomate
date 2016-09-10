# coding: utf-8

from __future__ import division, print_function, unicode_literals, absolute_import

import matplotlib.pyplot as plt

__author__ = 'Kiran Mathew'
__email__ = 'kmathew@lbl.gov'


def plot_wf(wf, depth_factor=1.0, breadth_factor=2.0, labels_on=True, numerical_label=False,
            text_loc_factor=1.0, save_as=None, style='rD--', markersize=10, markerfacecolor='blue',
            fontsize=12, ):
    """
    Generate a visual representation of the workflow. Useful for checking whether the firework
    connections are in order before launching the workflow.

    Args:
        wf (Workflow): workflow object.
        depth_factor (float): adjust this to stretch the plot in y direction.
        breadth_factor (float): adjust this to stretch the plot in x direction.
        labels_on (bool): whether to label the nodes or not. The default is to lable the nodes
            using the firework names.
        numerical_label (bool): set this to label the nodes using the firework ids.
        text_loc_factor (float): adjust the label location.
        save_as (str): save the figure to the given name.
        style (str): marker style.
        markersize (int): marker size.
        markerfacecolor (str): marker face color.
        fontsize (int): font size for the node label.
    """
    keys = sorted(wf.links.keys(), reverse=True)

    # set (x,y) coordinates for each node in the workflow links
    points_map = {}
    points_map.update({keys[0]: (-0.5*breadth_factor, (keys[0]+1)*depth_factor)})
    for k in keys:
        if wf.links[k]:
            for i, j in enumerate(wf.links[k]):
                if not points_map.get(j, None):
                    points_map[j] = ((i-len(wf.links[k])/2.0)*breadth_factor, k*depth_factor)

    # connect the dots
    for k in keys:
        for i in wf.links[k]:
            plt.plot([points_map[k][0], points_map[i][0]], [points_map[k][1], points_map[i][1]],
                     style, markersize=markersize, markerfacecolor=markerfacecolor)
            if labels_on:
                label1 = wf.id_fw[k].name
                label2 = wf.id_fw[i].name
                if numerical_label:
                    label1 = str(k)
                    label2 = str(i)
                plt.text(points_map[k][0] * text_loc_factor, points_map[k][1] * text_loc_factor,
                         label1, fontsize=fontsize)
                plt.text(points_map[i][0] * text_loc_factor, points_map[i][1] * text_loc_factor,
                         label2, fontsize=fontsize)

    plt.axis('scaled')
    plt.axis('off')

    if save_as:
        plt.savefig(save_as)
