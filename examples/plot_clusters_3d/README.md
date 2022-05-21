### Plot 3d clusters

This example creates figures to visualize the prototype cluster from each orbit. It reads the ``fcc_binary/basis_sets/bset.default/clust.json`` file generated by ``ccasm bset -uf``, plots each of the prototype clusters, and writes the figure as ``clusters_3d/orbit.i.pdf``, where ``i`` is an orbit index. Parameters at the top of the script control the project being read, the orbits being read, the figure size, and the output format.

To run this example:

    python plot_clusters_3d.py