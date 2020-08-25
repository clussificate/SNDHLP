# -*- coding: utf-8 -*-
"""
@Created at 2020/8/21 13:06
@Author: Kurt
@file: BPSolver.py
@Desc:
A branch-and-price-and-cut algorithm for service network design and hub location problem.
"""
import numpy as np
from pprint import pprint
import gurobipy as gp
from gurobipy import GRB
from queue import LifoQueue, PriorityQueue
from copy import deepcopy
import pickle
from collections import namedtuple
from SubprbSovler import SpSovler


class BPSolver:
    BEST_VAL = float("inf")
    _SOL = None

    def __init__(self, info):
        self.info = deepcopy(info)
        self.pending_nodes = LifoQueue()
        self.pending_nodes_best_bound = PriorityQueue()
        self.R = self.info['requests']

    def add_node(self, lpbound, node):
        self.pending_nodes.put(node)
        self.pending_nodes_best_bound.put((lpbound, node))

    def solve(self):
        # create root node
        initial_path = {}
        for r in self.R:
            initial_path[r] = []
            path_ = list(r)
            initial_path[r].append(path_)

        root = Node(self.info, initial_path, self)
        self.pending_nodes.put(root)
        self.pending_nodes_best_bound.put((1000000, root))

        iter_times = 0
        while not (self.pending_nodes.empty()):
            processing_node = self.pending_nodes.get()

            processing_node.process()
            iter_times += 1


class Node:
    def __init__(self, info, path, bpsolver):
        self.info = deepcopy(info)  # defensive copy
        self.k_t = self.info['capacity']
        self.H = self.info['hubs']
        self.H_to_label = {}
        self.label_to_H = {}
        self.path = path
        self.bpsolver = bpsolver

        # get all edges. A_t: bidirectional edge, A_s: normal monodirectional edge
        # self.edges = list((lambda x, y: {**x, **y})(self.info["A_t"], self.info["A_s"]))
        # self.edges.extend([(arc[1], arc[0]) for arc in self.info["A_t"]])

        for hub in self.H:
            ind = len(self.H_to_label)
            self.H_to_label[hub] = ind  # {'H0': 0}
            self.label_to_H[ind] = hub  # {0: 'H0'}

        self.r_to_label = {}
        self.label_to_r = {}
        self.R = self.info['requests']
        for req in self.R:
            ind = len(self.r_to_label)
            self.r_to_label[req] = ind
            self.label_to_r[ind] = req

        self.A_t = self.info["A_t"]
        self.at_to_label = {}
        self.label_to_at = {}
        for a_t in self.A_t:
            ind = len(self.at_to_label)
            self.at_to_label[a_t] = ind
            self.label_to_at[ind] = a_t

        self.A_s = self.info["A_s"]
        self.as_to_label = {}
        self.label_to_as = {}
        for a_s in self.A_s:
            ind = len(self.as_to_label)
            self.as_to_label[a_s] = ind
            self.label_to_as[ind] = a_s

        self.H_f = self.info['fixed_hubs']
        self.n_H = len(self.info['hubs'])

        self.c_t = self.info["rail_fee"]
        self.c_s = self.info["road_fee"]
        self.m = {key: val["demand"] for key, val in self.R.items()}  # {request: demand}

        # label for paths
        # path encoding ---->  {r_id:[0,1...]}: the p_th path of r_th request
        # path_to_label ----> {"0,1,3":(0, 1)}
        self.path_encoding = {}
        self.path_to_label = {}
        self.label_to_path = {}
        for r, paths in self.path.items():
            r_id = self.r_to_label[r]
            self.path_encoding[r_id] = []
            for path_ in paths:
                p_id = len(self.path_encoding[r_id])
                self.path_encoding[r_id].append(p_id)
                self.path_to_label[str(path_)] = (r_id, p_id)
                self.label_to_path[(r_id, p_id)] = path_

        # paths that contain hub i in each request: {h_id:{r_id:[p_id, p_id], r_id:[...]...}, hub2:{...}....}
        self.path_i = {hub_id: self.getPathByHub(hub_id) for hub_id in self.label_to_H}
        # pprint(self.path_i)

        # paths that contain arc: {arc_id:{r_id:[p_id, p_id], r_id:[...]...}, arc_id:{...}....}
        self.path_a_t = {arc_t_id: self.getPathByArc(arc_t_id, arc_type="t") for arc_t_id in
                         self.label_to_at}  # rail arc

        self.path_a_s = {arc_s_id: self.getPathByArc(arc_s_id, arc_type="s") for arc_s_id in
                         self.label_to_as}  # road arc
        # # alternative methods.
        # self.path_a_t = {self.label_to_at[arc_t_id]: self.getPathByArc(self.path, arc_t_id, arc_type="t") for
        # arc_t_id in self.label_to_at.keys()}  # rail arc self.path_a_s = {self.label_to_as[arc_s_id]:
        # self.getPathByArc(self.path, arc_s_id, arc_type="s") for arc_s_id in self.label_to_as.keys()}  # road arc
        # self.path_arc = (lambda x, y: {**x, **y})(self.path_a_t, self.path_a_s)  # merge path information

        self.master_prb = None
        self.price_prb = None
        self.solved = False

    def process(self):
        if not self.solved:
            print("Enter node")
            self._solve_lp()

    def _solve_lp(self):

        # self.price_prb = self.build_price_prb()

        # column generation loop
        max_iter_times = 1
        iter_time = 0
        while iter_time < max_iter_times:
            iter_time += 1
            master_prb, vars_tuple, constrs_tuple = self.build_master_prb()
            master_prb.setParam('outputflag', True)
            master_prb.optimize()

            # get dual variable value
            dual_c1 = master_prb.getAttr("Pi", constrs_tuple.c1)
            dual_c2 = constrs_tuple.c2.pi
            dual_c3 = master_prb.getAttr("Pi", constrs_tuple.c3)
            dual_c4 = master_prb.getAttr("Pi", constrs_tuple.c4)
            dual_c5 = master_prb.getAttr("Pi", constrs_tuple.c5)

            new_path = self._price(dual_c3, dual_c4, dual_c5)  # {request1: [path1, path2...]},  not id
            self.update_paras(new_path)

    def build_master_prb(self):
        # initial RLPM
        model = gp.Model("RLPM")

        # num_path = [(od, x) for od, paths in self.path.items() for x in range(paths)]  # ind for path variables
        y = model.addVars(list(self.label_to_path), lb=0.0, ub=1.0, vtype=GRB.CONTINUOUS, name="y")
        h = model.addVars(len(self.H), lb=0.0, ub=1, vtype=GRB.CONTINUOUS, name='h')
        t = model.addVars(len(self.A_t), lb=0.0, ub=GRB.INFINITY, vtype=GRB.CONTINUOUS, name="t")
        model.update()

        # add constraints
        c1 = model.addConstrs(
            (h[h_id] == 1 for h_id in self.label_to_H if self.label_to_H[h_id] in self.H_f))
        c2 = model.addConstr(sum([x for _, x in h.items()]), GRB.LESS_EQUAL,
                             self.n_H)  # note that the single constraint
        c3 = model.addConstrs((y.sum(r_id, "*") == 1 for r_id in self.label_to_r))

        c4 = model.addConstrs((sum([self.m[self.label_to_r[r_id]] * y[r_id, p_id]
                                    for r_id in self.label_to_r for p_id in
                                    self.path_a_t[arc_t_id][r_id]])
                               - self.k_t * t[arc_t_id] <= 0
                               for arc_t_id in self.label_to_at))

        c5 = model.addConstrs((sum([y[p_id] for p_id in self.path_i[h_id][r_id]]) <= h[h_id]
                               for h_id in self.label_to_H for r_id in self.label_to_r))  # key: (h_id, r_id)

        obj = gp.LinExpr()
        for arc_t_id, arc_t in self.label_to_at.items():
            obj += self.c_t * self.A_t[arc_t] * t[arc_t_id]

        # path_a_s ={arc_id:{r_id:[p_id]}}
        for arc_s_id, r_ids in self.path_a_s.items():
            arc_s = self.label_to_as[arc_s_id]
            l_a = self.A_s[arc_s]
            for r_id, p_id_list in r_ids.items():
                r = self.label_to_r[r_id]
                m_r = self.m[r]
                for p_id in p_id_list:
                    obj += self.c_s * l_a * m_r * y[r_id, p_id]

        model.setObjective(obj, GRB.MINIMIZE)

        Var = namedtuple('Variables', ['y', 'h', 't'])
        var = Var(y, h, t)

        Constr = namedtuple('Constraints', ['c1', 'c2', 'c3', 'c4', 'c5'])
        constr = Constr(c1, c2, c3, c4, c5)

        return model, var, constr

    def _price(self, *params):
        gammas, deltas, epsilons = params
        new_path = {}
        info = {"c_s": self.info["road_fee"], "m_r": {}, "starts": {}, "ends": {},
                "max_hop": self.info['max_hop'], "hubs": self.info['hubs'],
                "gamma_r": {}, "delta_at": {}, "l_as": {}, "epsilon_h_r": {}}

        for r, content in self.R.items():
            info["m_r"][r] = content["demand"]
            info["starts"][r] = content["starts"]
            info["ends"][r] = content["ends"]
            info["gamma_r"][r] = gammas[self.r_to_label[r]]

        for arc_s, length in self.info["A_s"].items():
            info["l_as"][arc_s] = length

        for arc_t, arc_t_id in self.at_to_label.items():
            info["delta_at"][arc_t] = deltas[arc_t_id]

        for hub, h_id in self.H_to_label.items():
            for r, r_id in self.r_to_label.items():
                info["epsilon_h_r"][hub, r] = epsilons[h_id, r_id]

        print("Current info for price problem:")
        pprint(info)

        # solve subproblem to find feasible path for each request.
        for r, _ in self.r_to_label.items():
            new_path[r] = SpSovler(r, info)
        return new_path

    def getPathByHub(self, hid):
        res = {}
        hub = self.label_to_H[hid]
        for r_id, p_ids in self.path_encoding.items():
            sel_p_id = []
            for p_id in p_ids:
                # print(r_id, p_id)
                p = self.label_to_path[r_id, p_id]
                if hub in p:
                    sel_p_id.append(p_id)
            res[r_id] = sel_p_id

        return res

    def getPathByArc(self, a_id, arc_type=None):
        res = {}
        if arc_type == "t":
            arc_map = self.label_to_at
        elif arc_type == "s":
            arc_map = self.label_to_as
        else:
            raise Exception()
        arc = arc_map[a_id]
        start, end = arc
        for r_id, p_ids in self.path_encoding.items():
            sel_pid = []
            for p_id in p_ids:
                p = self.label_to_path[r_id, p_id]
                for ind, node in enumerate(p[:-1]):
                    if p[ind] == start and p[ind + 1] == end:
                        sel_pid.append(p_id)
            res[r_id] = sel_pid
        return res

    def update_paras(self, new_path):
        """
        update:
        · self.label_to_path,  self.path_to_label, self.path_i, self.path_a_s, self.path_a_t
        · self.path , self.path_encoding
        """
        pass


if __name__ == "__main__":
    with open("info", "rb") as f:
        Info = pickle.load(f)

    bp = BPSolver(Info)
    bp.solve()
