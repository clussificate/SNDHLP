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
from collections import namedtuple, defaultdict
from queue import PriorityQueue
import math

# label
#     v: from the origin to a node v;     l_bar: previous label;     z: vector of resource consumptions
Label = namedtuple("label", ["v", "l_bar", "z"])

# vector of resource consumptions;
# we ignore transportation times in this implementation； hop can also be ignored by bidirectional_search
# cost: sum of the modified arc cost;
Consumption = namedtuple("consumptions", ['cost', 'hop'])


class Path:
    def __init__(self, path, cost):
        self.path = path
        self.cost = cost


class Node:
    def __init__(self, name, next_node=None):
        self.NAME = name
        self.NEXT = next_node


def order_comparison(label_x, label_y):
    pass


def dominance_comparison(label_x, label_y):
    pass


def parse_node(node):
    if node is None:
        return node
    else:
        return [x.NAME for x in node]


def parse_path(merged_label, deltas, gamma, m, epsilons):
    # print("-------------------Parse path information----------------------")
    paths = defaultdict(list)
    feasible_path = []
    for hop, pairs in merged_label.items():
        for pair in pairs:
            forward_label, backward_label = pair
            cost = forward_label.z.cost + backward_label.z.cost - gamma
            try:
                cost -= m * deltas[forward_label.v, backward_label.v]
            except KeyError:
                cost -= m * deltas[backward_label.v, forward_label.v]

            forward_path_reverse = []
            backward_path = []
            while forward_label:
                forward_path_reverse.append(forward_label.v)
                if forward_label.v in epsilons:
                    cost -= epsilons[forward_label.v]
                forward_label = forward_label.l_bar

            while backward_label:
                backward_path.append(backward_label.v)
                if backward_label.v in epsilons:
                    cost -= epsilons[backward_label.v]
                backward_label = backward_label.l_bar

            path = Path(forward_path_reverse[::-1] + backward_path, cost)
            paths[hop].append(path)
            # print("path:{}, cost:{}".format(path.path, cost))
            if cost < -0.001:
                feasible_path.append([cost, path])

    # for hop, paths in paths.items():
    #     for path in paths:
    #         print("Current hop: {}, current path: {}, current cost: {}".format(hop, path.path, path.cost))
    feasible_path.sort(key=lambda x: x[0], reverse=True)
    return paths, feasible_path


class SpSovler:
    def __init__(self, request, Info):
        self.info = deepcopy(Info)  # defensive copy
        self.request = request
        # self.edges = deepcopy(edges)  # defensive copy
        self.c_s = self.info['c_s']
        self.deltas = self.info['delta_at']
        self.epsilons = {hub_req[0]: epsilons
                         for hub_req, epsilons in self.info['epsilon_h_r'].items() if request == hub_req[1]}
        self.gamma_r = self.info["gamma_r"][request]
        self.l_a = self.info["l_as"]
        self.m_r = self.info["m_r"][request]
        self.starts_r = self.info["starts"][request]
        self.ends_r = self.info["ends"][request]
        # hop means the transmission between nodes, including origin to hub, and hub to destination.
        self.u = self.info["max_hop"]
        self.hubs = self.info['hubs']

        self.feasible_path = []
        self.path = []
        self.bidirectional_search(self.u)

    def monodirectional_search(self):
        pass

    def bidirectional_search(self, max_hop):
        # print("---------------Generate bidirectional search network----------------")
        origin = self.request[0]
        destination = self.request[1]
        starts = self.starts_r
        ends = self.ends_r
        forward_network = self.gen_network(origin, starts, max_hop, direction="forward")
        # print("Generate forward_network: \n{}".format({str(node.NAME) + "_" + str(layer_id): parse_node(node.NEXT)
        #                                                for layer_id, layer in enumerate(forward_network) for node in
        #                                                layer}))

        backward_network = self.gen_network(destination, ends, max_hop, direction="backward")
        # print("Generate backward_network: \n{}".format({str(node.NAME) + "_" + str(layer_id): parse_node(node.NEXT)
        #                                                 for layer_id, layer in enumerate(backward_network) for node in
        #                                                 layer}))
        forward_labels = self.gen_label(network=forward_network)
        backward_labels = self.gen_label(network=backward_network)
        paths, feasible_path = self.merge_label(forward_labels, backward_labels, max_hop)
        self.path = paths
        self.feasible_path = feasible_path

    def gen_network(self, source, first_layers, max_hop, direction):
        # print("---------------Direction: {}----------------".format(direction))
        hubs = self.hubs[:]
        layers = [[source], first_layers]
        number_of_layers = 0
        if direction == "forward":
            number_of_layers = math.floor(max_hop / 2) - 1  # number of copy layers
            layers.extend([hubs] * number_of_layers)
        if direction == "backward":
            number_of_layers = max_hop - math.floor(max_hop / 2) - 2  # number of copy layers
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
        # print("---------------Generate labels----------------")
        labels = []
        initial_root = network[0][0]
        initial_label = Label(initial_root.NAME, None, Consumption(0, 0))

        def dfs(root, label):
            if not root.NEXT:
                return
            else:
                for next_node in root.NEXT:
                    new_cost = self.get_cost(root.NAME, next_node.NAME) + label.z.cost
                    hop = label.z.hop
                    new_label = Label(next_node.NAME, label, Consumption(new_cost, hop + 1))
                    labels.append(new_label)
                    dfs(next_node, new_label)

        dfs(initial_root, initial_label)
        # print("example labels:")
        # label(v='H9',
        # l_bar=label(v='H4', l_bar=label(v='o0', l_bar=None, z=consumptions(cost=0)), z=consumptions(cost=764.9)),
        # z=consumptions(cost=1405.2))
        # od --> H4  --->， cost: 0-->764.9--->1405
        # print(labels[0])
        # print(labels[-1])
        return labels

    def get_cost(self, source, sink):
        if (source, sink) in self.l_a or (sink, source) in self.l_a:
            try:
                arc_cost = self.c_s * self.m_r * self.l_a[source, sink]
            except KeyError:
                arc_cost = self.c_s * self.m_r * self.l_a[sink, source]
        else:
            try:
                arc_cost = -self.m_r * self.deltas[source, sink]
            except KeyError:
                arc_cost = -self.m_r * self.deltas[sink, source]
        cost = arc_cost - self.epsilons[sink]
        return cost

    def merge_label(self, forward_labels, backward_labels, max_hop):
        # print("----------------------------Merge labels----------------------------")
        split_max_hop = max_hop - 1  # forward_network and backward network generate one hop.
        max_forward_hop = math.floor(max_hop / 2)
        max_backward_hop = split_max_hop - max_forward_hop
        minimal_hop = 3  # due to different start hub and end hub.
        split_minimal_hop = minimal_hop - 1

        merged_label = defaultdict(list)
        for allowed_hop in range(split_minimal_hop, split_max_hop + 1):
            for forward_hop in range(1, min(max_forward_hop + 1, allowed_hop)):
                backward_hop = allowed_hop - forward_hop
                if backward_hop > max_backward_hop:
                    continue
                # print(allowed_hop, forward_hop, backward_hop)

                sel_forward_labels = [label for label in forward_labels if label.z.hop == forward_hop]
                sel_backward_labels = [label for label in backward_labels if label.z.hop == backward_hop]
                # print("selected labels 1:\n{},\n{}".format(sel_forward_labels, len(sel_forward_labels)))
                # print("selected labels 2:\n{},\n{}".format(sel_backward_labels, len(sel_backward_labels)))

                for sel_forward_label in sel_forward_labels:
                    for sel_backward_label in sel_backward_labels:
                        if sel_forward_label.v != sel_backward_label.v:  # merge different hubs
                            # print("for \n{},\n length:{}".format(sel_forward_label, len(sel_forward_label)))
                            # print("back \n{}, \n length:{}".format(sel_backward_label, len(sel_backward_label)))
                            merged_label[allowed_hop + 1].append([sel_forward_label, sel_backward_label])

        # print("Get merged labels: \n")
        # pprint(dict(merged_label))
        # print("\n")
        # pprint(dict(merged_label)[3])
        # print("\n")
        # pprint(dict(merged_label)[5])

        paths = parse_path(merged_label, deltas=self.deltas, gamma=self.gamma_r, m=self.m_r, epsilons=self.epsilons)
        return paths


if __name__ == "__main__":
    info = DATA.sp_info
    pprint(info)
    r = ('o4', 'd4')
    sp = SpSovler(r, info)
    for cost, path_ in sp.feasible_path:
        print("feasible path {}, cost:{}".format(path_.path, path_.cost))
    print("----------------------------------------------------")
    for hop, paths_ in sp.path.items():
        for path in paths_:
            print("hop:{} all paths {}, cost:{}".format(hop, path.path, path.cost))
