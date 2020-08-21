# -*- coding: utf-8 -*-
"""
@Created at 2020/8/20 14:59
@Author: Kurt
@file:utils.py
@Desc:
"""
import hashlib

import numpy as np
import math
import random
from pprint import pprint
import networkx as nx
import matplotlib.pyplot as plt

random.seed(1994)


def random_color():
    r = lambda: random.randint(0, 255)
    return '#{:02x}{:02x}{:02x}'.format(r(), r(), r())


def viz(info: dict):
    G = nx.DiGraph()
    pos = {}

    # add edges
    for key, val in info['A_s'].items():
        G.add_edge(key[0], key[1], weight=val)

    for key, val in info['A_t'].items():
        G.add_edge(key[0], key[1], weight=val)
        G.add_edge(key[1], key[0], weight=val)

    # set edge color
    for road in info["A_s"].keys():
        G.edges[road[0], road[1]]['color'] = 'orange'

    for rail in info["A_t"].keys():
        G.edges[rail[0], rail[1]]['color'] = 'black'
        G.edges[rail[1], rail[0]]['color'] = 'black'

    edge_color_list = [G[e[0]][e[1]]["color"] for e in G.edges()]

    for od in info['A_t'].keys():
        o, d = od
        G.nodes[o]['color'] = 'b'
        G.nodes[d]['color'] = 'b'

    # set node color
    for od in info["requests"].keys():
        o, d = od
        color = random_color()
        print(color)
        G.nodes[o]['color'] = color
        G.nodes[d]['color'] = color

    node_color_list = [node[1]["color"] for node in G.nodes(data=True)]

    pos.update(info["location"])

    nx.draw(G, pos=pos, with_labels=True, node_size=200, node_color=node_color_list,
            edge_color=edge_color_list, font_size=8, font_color='white')
    plt.show()


def dist(a, b):
    return math.sqrt((a[0] - b[0]) ** 2 + (a[1] - b[1]) ** 2)


def gen_network(num_hub, num_fixed_hub):
    # generate hub in plane[0,100] * [0,50]
    info = {"requests": {}, "A_t": {}, "A_s": {}, "hubs": [], "fixed_hubs": [], "pot_hubs": {}}
    loc_index = {}
    length = 100
    height = 50
    hub_location = [[random.randint(0, length), random.randint(0, height)] for x in range(num_hub)]
    for i in range(num_hub):
        info["hubs"].append("H" + str(i))
        loc_index["H" + str(i)] = hub_location[i]
        for j in range(i + 1, num_hub):
            D = dist(hub_location[i], hub_location[j])
            info["A_t"][("H" + str(i), "H" + str(j))] = D

    info["fixed_hubs"] = np.random.choice(info['hubs'], num_fixed_hub, replace=False).tolist()
    info["pot_hubs"] = list(set(info["fixed_hubs"]) - set(info["pot_hubs"]))

    # create requests
    origins = [[25, 0], [0, 0], [0, 25], [0, 50], [25, 50]]
    dests = [[75, 50], [100, 50], [100, 25], [100, 0], [75, 0]]

    for ind, od in enumerate(zip(origins, dests)):
        D = dist(od[0], od[1])
        info["requests"][("o" + str(ind), "d" + str(ind))] = {"starts": {}, "ends": {}}
        info["A_s"][("o" + str(ind), "d" + str(ind))] = D
        loc_index["o" + str(ind)] = od[0]
        loc_index["d" + str(ind)] = od[1]

    # update allowed start points and end points for each request.
    for od, content in info["requests"].items():
        o, d = od
        dist_o = []
        dist_d = []
        choice_o = np.random.randint(1, 3)  # random choose one or two hubs as its started points
        choice_d = np.random.randint(1, 3)
        for hub in info["hubs"]:
            dist_o.append((hub, dist(loc_index[o], loc_index[hub])))
            dist_d.append((hub, dist(loc_index[d], loc_index[hub])))

        sort_choice_o = sorted(dist_o, key=lambda x: x[1])[0: choice_o]
        sort_choice_d = sorted(dist_d, key=lambda x: x[1])[0: choice_d]

        # update allowed starts and ends
        info["requests"][od]["starts"] = sort_choice_o
        info["requests"][od]["ends"] = sort_choice_d

        # update road arcs
        for hub in sort_choice_o:
            info["A_s"][(o, hub[0])] = hub[1]
        for hub in sort_choice_d:
            info["A_s"][(hub[0]), d] = hub[1]
    info['location'] = loc_index

    return info


if __name__ == "__main__":
    Info = gen_network(10, 5)
    pprint(Info)
    viz(Info)
