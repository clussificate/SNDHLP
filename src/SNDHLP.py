# -*- coding: utf-8 -*-
"""
@Created at 2020/8/21 13:06
@Author: Kurt
@file:SNDHLP.py
@Desc:
"""
import numpy as np
import gurobipy as gp
from gurobipy import GRB
import pickle


class SNDHLP:
    def __init__(self, info):
        self.k_t = info['capacity']
        self.H = info['hubs']
        self.H_to_label = {}
        self.label_to_H = {}
        for hub in self.H:
            ind = len(self.H_to_label)
            self.H_to_label[hub] = ind  # {'H0': 0}
            self.label_to_H[ind] = hub  # {0: 'H0'}

        self.r_to_label = {}
        self.label_to_r = {}
        self.R = info['requests']
        for req in self.R.keys():
            ind = len(self.r_to_label)
            self.r_to_label[req] = ind
            self.label_to_r[ind] = req

        self.A_t = info["A_t"]
        self.at_to_label = {}
        self.label_to_at = {}
        for a_t in self.A_t.keys():
            ind = len(self.at_to_label)
            self.at_to_label[a_t] = ind
            self.label_to_at[ind] = a_t

        self.A_s = info["A_s"]
        self.as_to_label = {}
        self.label_to_as = {}
        for a_s in self.A_s:
            ind = len(self.as_to_label)
            self.as_to_label[a_s] = ind
            self.label_to_as[ind] = [a_s]

        self.H_f = info['fixed_hubs']
        self.n_H = len(info['hubs'])

        self.c_t = info["rail_fee"]
        self.c_s = info["road_fee"]

        self.m = {key: val["demand"] for key, val in self.R.items()}

        # initial RLPM
        self.model = gp.Model("RLPM")

        # label for paths (r_id, p_id): the p_th path of t_th request
        self.path = {}
        self.path_to_label = {}
        self.label_to_path = {}
        for r in self.R.keys():
            r_id = self.r_to_label[r]
            self.path[r_id] = []
            p_id = len(self.path[r_id])
            self.path[r_id].append(p_id)
            path = list(r)
            self.path_to_label[str(path)] = (r_id, p_id)
            self.label_to_path[(r_id, p_id)] = path

        # paths that contain hub i in each request: {h_id:{r_id:[p_id, p_id], r_id:[...]...}, hub2:{...}....}
        self.path_i = {hub_id: self.getPathByHub(self.path, hub_id) for hub_id in self.label_to_H.keys()}

        # paths that contain arc: {arc_id:{r_id:[p_id, p_id], r_id:[...]...}, arc_id:{...}....}

        self.path_a_t = {arc_t_id: self.getPathByArc(self.path, arc_t_id, type= "t") for arc_t_id in self.label_to_at.keys()}  # rail arc
        self.path_a_s = {arc_s_id: self.getPathByArc(self.path, arc_s_id, type= "s") for arc_s_id in self.label_to_as.keys()}  # road arc
        self.path_arc = (lambda x, y: {**x, **y})(self.path_a_t, self.path_a_s)  # merge path information

        # num_path = [(od, x) for od, paths in self.path.items() for x in range(paths)]  # ind for path variables
        self.y = self.model.addvars(list(self.label_to_path.keys()), lb=0.0, ub=1.0, vtype=GRB.CONTINUOUS, name="y")
        self.h = self.model.addVars(len(self.H), lb=0.0, ub=1, vtype=GRB.CONTINUOUS, name='h')
        self.t = self.model.addVars(len(self.A_t), lb=0.0, ub=GRB.INFINITY, vtype=GRB.CONTINUOUS, name="t")
        self.model.update()

        # add constraints
        self.c1 = self.model.addConstrs((self.h[self.H_to_label[hf]] == 1 for hf in self.H_f))
        self.c2 = self.model.addConstr(sum([x for _, x in self.h.items()]), GRB.LESS_EQUAL, self.n_H)
        self.c3 = self.model.addConstrs((self.y.sum(self.r_to_label[r], "*") == 1 for r in self.R.keys()))

        self.c4 = self.model.addConstrs((sum([self.m[r_id]*self.y[p_id]
                                              for r_id in self.label_to_r.keys() for p_id in self.path_a_t[arc_t_id][r_id]])
                                         <= self.k_t*self.t[arc_t_id]
                                         for arc_t_id in self.label_to_at))

        self.c5 = self.model.addConstrs((sum([self.y[p_id] for p_id in self.path_i[h_id][r_id]]) <= self.h[h_id] for h_id in self.label_to_H for r_id in self.label_to_r))

    def getPathByHub(self, path, hid):
        res = {hid: {}}
        hub = self.label_to_H[hid]
        for r_id, p_ids in path.items():
            sel_p_id = []
            for p_id in p_ids:
                path = self.label_to_path[p_id]
                if hub in path:
                    sel_p_id.append(p_id)
            res[hid][r_id] = sel_p_id

        return res

    def getPathByArc(self, path, a_id, type=None):
        res = {a_id: {}}
        if type== "t":
            arc_map = self.label_to_at
        elif type == "s":
            arc_map = self.label_to_as
        else:
            raise Exception()
        arc = arc_map[a_id]
        start, end = arc
        for r_id, p_ids in path.items():
            sel_pid = []
            for p_id in p_ids:
                path = self.label_to_path[p_id]
                for ind, node in enumerate(path[:-1]):
                    if path[ind] == start and path[ind+1] ==end:
                        sel_pid.append(p_id)
            res[a_id][r_id] = sel_pid
        return res


if __name__ == "__main__":
    with open("info", "rb") as f:
        Info = pickle.load("info")
