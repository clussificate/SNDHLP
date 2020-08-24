# -*- coding: utf-8 -*-
"""
@Created at 2020/8/23 17:48
@Author: Kurt
@file:SubprbSovler.py
@Desc:
the subproplem is a shortest path problem with resource constraints, SPPRC.
"""
from test_data import DATA
from pprint import pprint
from copy import deepcopy
from collections import namedtuple
import math
from queue import PriorityQueue

# label
#     v: from the origin to a node v;     l_bar: previous label;     z: vector of resource consumptions
Label = namedtuple("label", ["v", "l_bar", "z"])

# vector of resource consumptions;
# we ignore transportation times in this implementation； hop can also be ignored by bidirectional_search
# cost: sum of the modified arc cost;
Consumption = namedtuple("consumptions", ['cost'])


def order_comparison(label_x, label_y):
    pass


def dominance_comparison(label_x, label_y):
    pass


def parse_node(node):
    if node is None:
        return node
    else:
        return [x.NAME for x in node]


class SpSovler:
    def __init__(self, request, Info, edges):
        self.info = deepcopy(Info)  # defensive copy
        self.request = request
        self.edges = deepcopy(edges)  # defensive copy
        self.c_s = self.info['c_s']
        self.deltas = self.info['delta_at']
        self.epsilons = {hub_req[0]: epsilons
                         for hub_req, epsilons in self.info['epsilon_h_r'].items() if request == hub_req[1]}
        self.gamma_r = self.info["gamma_r"][request]
        self.l_a = self.info["l_as"]
        self.m_r = self.info["m_r"][request]
        self.starts_r = self.info["starts"][request]
        self.ends_r = self.info["ends"][request]
        self.u = self.info["max_hop"]
        self.hubs = self.info['hubs']

        self.shortest_path = []

    def monodirectional_search(self):
        pass

    def bidirectional_search(self):
        print("---------------Generate bidirectional search network----------------")
        origin = self.request[0]
        destination = self.request[1]
        starts = self.starts_r
        ends = self.ends_r
        u = self.u
        forward_network = self.gen_network(origin, starts, u, direction="forward")
        print("Generate forward_network: \n{}".format({str(node.NAME) + "_" + str(layer_id): parse_node(node.NEXT)
                                                       for layer_id, layer in enumerate(forward_network) for node in
                                                       layer}))

        backward_network = self.gen_network(destination, ends, u, direction="backward")
        print("Generate backward_network: \n{}".format({str(node.NAME) + "_" + str(layer_id): parse_node(node.NEXT)
                                                        for layer_id, layer in enumerate(backward_network) for node in
                                                        layer}))

        self.gen_label(network=forward_network)

    def gen_network(self, source, first_layers, max_hop, direction):
        print("---------------Direction: {}----------------".format(direction))
        hubs = self.hubs
        layers = [[source], first_layers]
        number_of_layers = 0
        if direction == "forward":
            number_of_layers = math.floor(max_hop / 2)  # number of copy layers
            layers.extend([hubs] * number_of_layers)
        if direction == "backward":
            number_of_layers = (max_hop + 1) - (math.floor(max_hop / 2) + 2)  # number of copy layers
            layers.extend([hubs] * number_of_layers)

        network = [[Node(name=node) for node in layer] for layer in layers]
        total_num_layer = len(network)
        for layer_ind, layer in enumerate(network):
            for node_id, node in enumerate(layer):
                if layer_ind < total_num_layer - 1:
                    next_node_list = network[layer_ind + 1][:]  # defensive copy
                    for next_node in next_node_list:
                        if node.NAME == next_node.NAME:
                            next_node_list.remove(next_node)

                    network[layer_ind][node_id].NEXT = next_node_list
        return network

    def gen_label(self, network):
        print("---------------Generate labels----------------")
        labels = PriorityQueue()
        initial_root = network[0][0]
        initial_label = Label(initial_root.NAME, None, Consumption(0))

        def dfs(root, label):
            if not root.NEXT:
                return
            else:
                for next_node in root.NEXT:
                    new_cost = label.z.cost + 1
                    new_label = Label(next_node.NAME, label, Consumption(new_cost))
                    labels.put([new_label.z.cost, new_label])   # put PriorityQueue by cumulative cost.
                    dfs(next_node, new_label)

        dfs(initial_root, initial_label)
        print(labels)

        return labels


class Node:
    def __init__(self, name, next_node=None):
        self.NAME = name
        self.NEXT = next_node


if __name__ == "__main__":
    info = DATA.sp_info
    # pprint(info)
    r = ('o0', 'd0')
    edges = DATA.sp_edge
    sp = SpSovler(r, info, edges)
    sp.bidirectional_search()
