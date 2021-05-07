#!/usr/bin/env python

import numpy
import scipy.stats

#The Savage-Dickey estimator computation is based on the implementation by Äijö et al. (2016) available at https://github.com/tare/LuxGLM (MIT lisence)

def calculate_savagedickey(prior1_mean,prior1_cov,prior2_mean,prior2_cov,samples1,samples2):
  samples1_mean = numpy.mean(samples1,0)
  samples1_cov = numpy.cov(samples1,rowvar=0)

  samples2_mean = numpy.mean(samples2,0)
  samples2_cov = numpy.cov(samples2,rowvar=0)

  numerator = scipy.stats.multivariate_normal.pdf(numpy.zeros(prior1_mean.shape),mean=prior1_mean-prior2_mean,cov=prior1_cov+prior2_cov)
  denominator = scipy.stats.multivariate_normal.pdf(numpy.zeros(prior1_mean.shape),mean=samples1_mean-samples2_mean,cov=samples1_cov+samples2_cov)

  return numerator/denominator 

def calculate_savagedickey_kde(prior1_mean,prior1_cov,prior2_mean,prior2_cov,samples1,samples2):
  Delta = samples1-samples2
  density = scipy.stats.kde.gaussian_kde(Delta,bw_method='scott')

  numerator = scipy.stats.multivariate_normal.pdf(numpy.zeros(prior1_mean.shape),mean=prior1_mean-prior2_mean,cov=prior1_cov+prior2_cov)
  denominator = density.evaluate([0])[0]

  return numerator/denominator

def calculate_savagedickey_kde_1d(prior_mean,prior_cov,samples):
  density = scipy.stats.kde.gaussian_kde(samples.T,bw_method='scott')

  numerator = scipy.stats.multivariate_normal.pdf(numpy.zeros(prior_mean.shape),mean=prior_mean,cov=prior_cov)
  denominator = density.evaluate([0,0])[0]

  return numerator/denominator

def calculate_savagedickey_kde_window(prior1_mean,prior1_cov,prior2_mean,prior2_cov,samples1,samples2):
  #samples1 and samples2 have shape (# of dims, # of samples). prior1_mean and prior2_mean have shape (#dim,1) and prior1_cov and prior2_cov have shape (#dim,#dim)
  Delta = samples1-samples2
  density = scipy.stats.kde.gaussian_kde(Delta,bw_method='scott')

  numerator = scipy.stats.multivariate_normal.pdf(numpy.zeros(Delta.shape[0]),mean=prior1_mean-prior2_mean,cov=prior1_cov+prior2_cov)
  denominator = density.evaluate(numpy.zeros(Delta.shape[0]))

  return numerator/denominator, numerator, denominator

