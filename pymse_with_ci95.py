#!/usr/bin/env python
#coding:utf-8
"""
  Author:  Yeison Cardona --<yeison.eng@gmail.com>
  Modified by John Roberts --<jfr10@leicester.ac.uk>
  Modified by Heiko Balzter --<hb91@le.ac.uk>
  Purpose:
  Created: 12/02/16
"""

import os
import numpy as np
from scipy.stats import t
import collections


########################################################################
class PyMSE:
    """"""

    #----------------------------------------------------------------------
    def __init__(self, dataset, scale=20, m=2, r=0.15, q=0.975):
        """Constructor"""
        if isinstance(dataset, collections.Iterable):
            self.data = np.array(dataset)
        elif os.path.exists(dataset):
            self.data = read_data(dataset)
        self.data_cg = self.data.copy()
        self.get(scale, m, r, q)  # Everything in get() should be in init()


    #----------------------------------------------------------------------
    def get(self, scale=20, m=2, r=0.15, q=0.975):

        assert isinstance(scale, (int, list, tuple, map, range, np.ndarray)), "scale must be int, list or tuple."
        assert isinstance(m, (int, list, tuple, map, range, np.ndarray)), "m must be int, list or tuple."
        assert isinstance(r, (int, float, list, tuple, map, range, np.ndarray)), "r must be int, float, list or tuple."
        assert isinstance(q, (float)), "q must be float."

        #Scales
        if isinstance(scale, int):
            self.SCALE = [scale]
        elif isinstance(scale, (list, tuple, range, np.ndarray)):
            self.SCALE = scale
        else:
            self.SCALE = range(1, 21, 1)
        self.scale_max = self.SCALE[-1]
        if len(self.SCALE) > 1:
            self.scale_step = self.SCALE[1] - self.SCALE[0]
        else:
            self.scale_step = 1


        # So, I _think_ that M is the range of resolutions for MSE calcs.
        if isinstance(m, int):
            self.M = [m]
        elif isinstance(m, (list, tuple, range, np.ndarray)):
            self.M = m
        else:
            self.M = range(2, 3, 1)  # This is a very roundabout way of saying '2'.
        self.m_min = self.M[0]
        self.m_max = self.M[-1]
        if len(self.M) > 1:
            self.m_step = self.M[1] - self.M[0]
        else:
            self.m_step = 1

        if isinstance(r, (int, float)):
            self.R = np.arange(r, r*1.0000000001, 0.05)
        elif isinstance(r, list, tuple, range, np.ndarray):
            self.R = r
        else:
            self.R = np.arange(r, r*1.0000000001, 0.05)
        self.r_min = self.R[0]
        self.r_max = self.R[-1]
        if len(self.R) > 1:
            self.r_step = self.R[1] - self.R[0]
        else:
            self.r_step = 0.05

        standard_deviation = self.data.std()

        SE = {}
        for sc in self.SCALE:
            self.__coarse_graining__(sc)

            for r in self.R:
                se = self.sample_entropy(r, standard_deviation, sc);
                if not r in SE:
                    SE[r] = []

                SE[r].append(se)


        self.DATA = []
        for se_r in SE.keys():

            ent = list(map(lambda *arg:arg, *SE[se_r]))

            for m_ in range(len(self.M)):
                entr = dict(map(lambda *arg:arg, self.SCALE, ent[m_]))
                self.DATA.append({"m": self.M[m_], "mse": entr, "r": se_r, "q": q})


        if len(self.DATA) == 1:
            return self.DATA[0]
        else:
            return self.DATA






    ##----------------------------------------------------------------------
    #def __read_data__(self, dataset):
        #""""""
        #if isinstance(dataset, (str, bytes)):
            #assert os.path.isfile(dataset), "Missing \"{}\" file.".format(dataset)

            #with open(dataset, "r") as file:
                #dataset = np.array(list(map(float, file.readlines())))
                #file.close()

        #return dataset



    #----------------------------------------------------------------------
    def __coarse_graining__(self, resolution):
        out_len = int(len(self.data) / resolution)
        out = np.empty([out_len, 1])
        for i, subarray in enumerate(np.split(self.data, out_len)):
            out[i] = np.mean(subarray)
        self.data_cg = out


    #----------------------------------------------------------------------
    def sample_entropy(self, r, standard_deviation, scale=1, q=0.975):
        """"""
        se = []
        ci = []

        nlin = float(len(self.data))
        nlin_j = int((nlin/scale) - self.m_max)
        r_new = r*standard_deviation

        cont = [0] * (self.m_max+2)

        for i in range(0, nlin_j):
            for l in range(i+1, nlin_j):
                k = 0
                while k < self.m_max and (np.abs(self.data_cg[i+k] - self.data_cg[l+k]) <= r_new):
                    k += 1
                    cont[k] += 1
                if k == self.m_max and (np.abs(self.data_cg[i+self.m_max] - self.data_cg[l+self.m_max]) <= r_new):
                    cont[self.m_max+1] += 1

        for i in self.M:
            if cont[i+1] == 0 or cont[i] == 0:
                if (nlin_j > 0) and ((nlin_j-1) > 0):
                    se.append(-1 * np.log(1.0/(nlin_j*(nlin_j-1))))
                    ci.append(standard_deviation * t.ppf(q, cont[i] - 1) / np.sqrt(cont[i]))
                else: 
                    se.append(0.0)
                    ci.append(0.0)
            else:
                se.append(-1 * np.log(float(cont[i+1])/cont[i]))
                ci.append(standard_deviation * t.ppf(q, cont[i] - 1) / np.sqrt(cont[i]))
        return se, ci


    def conf_int(self, r, standard_deviation, scale=1):
        """

        Parameters
        ----------
        r
        standard_deviation
        scale

        Returns
        -------

        """
        ci95 = []
        
        nlin = float(len(self.data))
        nlin_j = int((nlin/scale) - self.m_max)
        r_new = r * standard_deviation

        cont = [0] * (self.m_max+2)

        for i in range(0, nlin_j):
            for l in range(i+1, nlin_j):
                k = 0
                while k < self.m_max and (np.abs(self.data_cg[i+k] - self.data_cg[l+k]) <= r_new):
                    k += 1
                    cont[k] += 1
                if k == self.m_max and (np.abs(self.data_cg[i+self.m_max] - self.data_cg[l+self.m_max]) <= r_new):
                    cont[self.m_max+1] += 1

        for i in self.M:
            if cont[i+1] == 0 or cont[i] == 0:
                if (nlin_j > 0) and ((nlin_j-1) > 0):
                    # calculate the confidence interval of the sample entropy at scale j after Richman and Moorman (2000)
                    # We have B template matches of which A actually occur
                    # Assign 1 to the A forward matches and 0 to the B-A potential forward matches that do not occur
                    # The 95% confidence interval is then: SD * t(B-1, 0.975) / sqrt(B)
                    # where SD is the standard deviation of the time-series
                    # Here, A = cont[i]+1 and B = cont[i]
                    ci95.append(np.std(self[1:nlin_j]) * t.ppf(0.975, cont[i]-1) / np.sqrt(cont[i]))
                else: 
                    ci95.append(0.0)
            else:
                ci95.append(0.0)
                
        return ci95


    #----------------------------------------------------------------------
    def __standard_deviation__(self):
        """"""
        nlin = float(len(self.data))
        sum_ = sum(self.data)
        sum2_ = sum(self.data*self.data)

        return np.sqrt((sum2_ - sum_*(sum_/nlin))/(nlin - 1))

#----------------------------------------------------------------------
def read_data(dataset):
    """"""
    if isinstance(dataset, (str, bytes)):
        assert os.path.isfile(dataset), "Missing \"{}\" file.".format(dataset)

        with open(dataset, "r") as file:
            dataset = np.array(list(map(float, file.readlines())))
            file.close()

    return dataset