#This code was derived from the R-package 'Specs Verification' (Function: BrierDecomp) 
#by Stefan Siegert, Jonas Bhend, Igor Kroener and Matteo De Felice 
#Dept. Mathematics and Statistics, University of Exeter, UK
  #Title: Forecast Verification Routines for Ensemble Forecasts of Weather and Climate
  #CRAN Repository: https://rdrr.io/cran/SpecsVerification/src/R/BrierDecomp.R 
  #Date/Publication: 2020-02-26 15:40:06 UTC
    

import sys
import numpy as np
import math
import pandas as pd

import warnings
warnings.filterwarnings("ignore", category=RuntimeWarning)

def BrierDecomp(p, y, bins=10, bias_corrected=False):
# Brier Score decomposition
# Returns 3-components of the Brier Score: Reliability, Resolution and Uncertainty, and their standard deviations
# INPUT:
    # p: 1d np.array of predicted probabilities
    # y: 1d np.array of binary observations (0 or 1)
    # bins: to estimate the calibration function, default: 10
    # bias_corrected: logical, default=False, if false, the standard (biased) decomposition of Murphy (1973) is used, 
                                             #if true, the bias-corrected decomposition of Ferro (2012) (See References in Readme)
#OUTPUT:
    #1) dictionary with estimates of Reliability, Resolution and Uncertainty
    #2) dictionary with corresponding std. deviations
    
    
    n=len(p)

    # binning from either specified number of bins or breaks 
    if np.array(bins).size==1:
        n_bins = math.floor(bins)
        p_breaks = np.linspace(start=0, stop=1, num=n_bins+1)
    else:
        n_bins = len(bins) #- 1
        bins = np.sort(bins)
        if(min(bins)<= 0 and max(bins)>= 1): #R: stopifnot(min(bins)<= 0 & max(bins) >= 1)
            sys.exit()
        p_breaks = np.append(np.insert(bins,0,0.),1.)  #start from 0 and end with 1

    #p.binning (vector of length n with entries equal to bin-no. of p_n)
    p_binning = np.digitize(p, p_breaks, right=True)
    
    ##construct matrices and column sums
    #Define index, to set 1s. -1, because R starts at 1, and python at 0.
    m_ind = np.hstack([np.arange(n).reshape(-1,1), (p_binning-1).reshape(-1,1)])  #### Do I need to transform asmatrix?
    m_A = np.zeros(shape=(n,n_bins))
    pos = list(map(tuple, m_ind))
    rows, cols = zip(*pos)
    m_A[rows, cols] = 1
    cs_A = np.sum(m_A, axis=0)
    m_B = np.zeros(shape=(n,n_bins))
    m_B[rows, cols] = y
    cs_B = np.sum(m_B, axis=0)
    m_C = np.zeros(shape=(n,n_bins))
    m_C[rows, cols] = p
    cs_C = np.sum(m_C, axis=0)
    m_Y = y.reshape(-1,1)
    cs_Y = np.sum(m_Y, axis=0)
    
    ## sets of indices for which column-sums(A) > 0, resp. > 1
    d0 = np.ravel(np.where(cs_A<1))
    d1 = np.ravel(np.where(cs_A < 2))
    
    ## in order to avoid division by zero
    def d0_set(x):
        if len(d0)>0:
            x = np.delete(x, d0)
        else: x
        return x

    def d1_set(x):
        if len(d1)>0:
            x = np.delete(x, d1)
        else: x
        return x

    #R Source CODE (SpecsVerification) as reference:
    #rel <- 1/n * sum( d0.set((cs.B - cs.C)^2 / cs.A) )
    #res <- 1/n * sum( d0.set(cs.A * (cs.B / cs.A - cs.Y / n)^2) )
    #unc <- cs.Y * (n - cs.Y) / n / n

    rel = (1/n)*np.sum(d0_set(((cs_B-cs_C)**2)/cs_A))
    res = (1/n) * np.sum(d0_set(cs_A * (cs_B / cs_A - cs_Y / n)**2))
    unc = cs_Y * (n-cs_Y) / n / n
    
    # correction factors for bias correction
    corr_s = 1 / n * np.sum(d1_set(cs_B * (cs_A - cs_B) / (cs_A * (cs_A - 1))))
    corr_t = cs_Y * (n-cs_Y)/n/n/(n-1)

    # avoid rel<0, res<0, res>1, and unc>1/4
    alpha = np.min(np.array([rel/corr_s, 
                     np.max(np.array([res / (corr_s - corr_t), (res - 1) / (corr_s - corr_t)])), 
                             (1 - 4 * unc) / (4 * corr_t), 
                             1]))

    rel2 = rel - alpha * corr_s
    res2 = res - alpha * corr_s + alpha * corr_t
    unc2 = unc + alpha * corr_t
    
    ###################################################
    # Estimating variances of reliability, resolution, uncertainty by propagation of error
    ###################################################

    #####
    # REL
    #####
    d_rel_d_a = -1 * (cs_B - cs_C) * (cs_B - cs_C) / (n * cs_A * cs_A)
    d_rel_d_b = 2 * (cs_B - cs_C) / (n * cs_A)
    d_rel_d_c = -2 * (cs_B - cs_C) / (n * cs_A)
    if (len(d0) > 0):
        d_rel_d_a[d0] = 0
        d_rel_d_b[d0] = 0
        d_rel_d_c[d0] = 0

    jacobian_rel = np.concatenate((d_rel_d_a, d_rel_d_b, d_rel_d_c))
    m_X = np.hstack([m_A, m_B, m_C])
    #center
    col_means = np.mean(m_X, axis=0)
    centered_X = m_X - col_means
    cov_X_rel = np.dot(centered_X.T, centered_X)
    var_rel = np.matmul(np.matmul(jacobian_rel, cov_X_rel), jacobian_rel.T)

    #####
    # RES
    #####
    d_res_d_a = -1 / n * (cs_B / cs_A - cs_Y / n) * (cs_B / cs_A + cs_Y / n)
    d_res_d_b = 2 / n * (cs_B / cs_A - cs_Y / n)
    d_res_d_y = 0
    if (len(d0) > 0):
        d_res_d_a[d0] = 0
        d_res_d_b[d0] = 0

    jacobian_res = np.append(np.concatenate((d_res_d_a, d_res_d_b)), d_res_d_y)
    m_X = np.hstack([m_A, m_B, m_Y])

    #center
    col_means = np.mean(m_X, axis=0)
    centered_X = m_X - col_means
    cov_X_res = np.dot(centered_X.T, centered_X)
    var_res = np.matmul(np.matmul(jacobian_res, cov_X_res), jacobian_res.T)


    #####
    # UNC
    #####
    var_unc = np.matmul((1 - 2 * cs_Y / n) * (1  - 2 * cs_Y / n) / n / n , 
                       np.dot((m_Y - np.mean(m_Y)).T, (m_Y - np.mean(m_Y))))
    
    ######
    # REL BIAS CORRECTED
    ######
    d_rel2_d_a = -1 * ((cs_B - cs_C) * (cs_B - cs_C) + (cs_B * cs_B) / (cs_A - 1) - cs_A * cs_B * (cs_A - cs_B) / (cs_A - 1) / (cs_A - 1)) / (n * cs_A * cs_A)
    d_rel2_d_b = (2 * cs_B - 1) / (n * cs_A - n) - 2 * cs_C / (n * cs_A)
    d_rel2_d_c = -2 * (cs_B - cs_C) / (n * cs_A)
    if(len(d1) > 0):
        d_rel2_d_a[d1] = 0
        d_rel2_d_b[d1] = 0
        d_rel2_d_c[d1] = 0

    jacobian_rel = np.concatenate((d_rel2_d_a, d_rel2_d_b, d_rel2_d_c))
    var_rel2 = np.matmul(np.matmul(jacobian_rel, cov_X_rel), jacobian_rel.T)

    ######
    # RES  BIAS CORRECTED
    ######

    d_res2_d_a = -1 / n * (cs_B / cs_A - cs_Y / n) * (cs_B / cs_A + cs_Y / n) + cs_B / (n * cs_A * cs_A * (cs_A - 1) * (cs_A - 1)) * ((cs_A - cs_B) * (cs_A - cs_B) - cs_B * (cs_B - 1))
    d_res2_d_b = 2 / n * (cs_B / cs_A - cs_Y / n) - (cs_A - 2 * cs_B) / (n * cs_A * (cs_A - 1))
    d_res2_d_y = (n - 2 * cs_Y) / n / n / (n - 1)
    if (len(d1) > 0):
        d_res2_d_a[d1] = 0
        d_res2_d_b[d1] = 0

    jacobian_res2 = np.append(np.concatenate((d_res2_d_a, d_res2_d_b)), d_res2_d_y)
    var_res2 = np.matmul(np.matmul(jacobian_res2, cov_X_res), jacobian_res2.T)
    
    ######
    # UNC BIAS CORRECTED
    ######

    var_unc2 = np.matmul((n - 2 * cs_Y) * (n - 2 * cs_Y) / n / n / (n - 1) / (n - 1), 
                         np.dot((m_Y - np.mean(m_Y)).T, (m_Y - np.mean(m_Y))))
    
    #Define sqrt(0) as 0 to avoid errors.
    def sqrt0(x):
        if (x <= 0):
            return(0)
        else:
            return(np.sqrt(x))

    
    if(bias_corrected == False):
        est_dict = {'res': float(res), 'rel': rel, 'unc': unc}
        sd_dict = {'res': float(sqrt0(var_res)), 'rel': float(sqrt0(var_rel)), 'unc': float(sqrt0(var_unc))}

    else:
        est_dict = {'res': float(res2), 'rel': float(rel2), 'unc': float(unc2)}
        sd_dict = {'res': float(sqrt0(var_res2)), 'rel': float(sqrt0(var_rel2)), 'unc': float(sqrt0(var_unc2))}


    return est_dict, sd_dict
