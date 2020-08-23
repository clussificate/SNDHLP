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

# label
#     v: from the origin to a node v;     l_bar: previous label;     z: vector of resource consumptions
Label = namedtuple("label", ["v", "l_bar", "z"])

# vector of resource consumptions;
# we ignore transportation times in this implementation.
# cost: sum of the modified arc cost;    hop: number of hops
Consumption = namedtuple("resource consumptions", ['cost', 'hop'])


def order_comparison(label_x, label_y):
    pass


def dominance_comparison(label_x, label_y):
    pass


class SpSovler:
    def __init__(self, request, Info, edges):
        self.info = deepcopy(Info)  # defensive copy
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
        self.shortest_path = []

    def monodirectional_search(self):
        pass

    def bidirectional_search(self):
        pass


if __name__ == "__main__":
    info = DATA.sp_info
    # pprint(info)
    r = ('o0', 'd0')
    edges = DATA.sp_edge
    sp = SpSovler(r, info, edges)
