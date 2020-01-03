import argparse
import sys
import os
import commands
import pystan
import numpy
import scipy
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot
from pylab import plot, show, savefig, xlim, figure, hold, ylim, legend, boxplot, setp, axes, scatter

import pickle
from hashlib import md5
import logging

import time
from timeit import default_timer as timer
import subprocess
import csv

from savagedickey_spatial2_wholeWindow import calculate_savagedickey, calculate_savagedickey_kde, calculate_savagedickey_kde_window, calculate_savagedickey_kde_1d


#This script is been used for the real data experiments.
#In this version the sigmaB values have been fixed for v1, v3 and v3 sep 17.9.2018



def stan_cache(model_name, **kwargs):
  f=open(model_name, 'rb')
  model_code=f.read()
  f.close()
  code_hash = md5(model_code.encode('ascii')).hexdigest()
  cache_fn = 'cached-{}-{}.pkl'.format(model_name, code_hash)
  try:
    sm = pickle.load(open(cache_fn, 'rb'))
  except:
    sm = pystan.StanModel(file=model_name)
    with open(cache_fn, 'wb') as f:
      pickle.dump(sm, f)
  else:
    logging.info("Using cached StanModel")
  return sm.sampling(**kwargs)

def csvIntoExtractDict(csvFile,isDiagnostic):
        #Inputs:
        #csvfile: the name of the file from which data is extracted
        #isDiagnostic: logical, use when retrieving ELBO values from diagnostic files
        outputdict={}
        with open(csvFile) as cmdfile:
                reader=csv.reader(cmdfile)
                isFirstrow=0
                for row in reader:
                        #print('Current row')
                        #print(row)
                        if row[0].startswith("#"):
                                #print('Row started with #')
                                continue
                        N_rowitems=len(row)
                        #print('Printing isFirstrow')
                        #print(isFirstrow)
                        if isFirstrow==0:
                                variableNames=[]
                                #print('Initializing dictionary')
                                if isDiagnostic==1: #if the input file is a diagnostic file, the dictionary keys can't be found$
                                        outputdict['iter']=[] #Makes a dummy row of zeros as the first row
                                        outputdict['time_in_seconds']=[]
                                        outputdict['ELBO']=[]
                                        variableNames=['iter','time_in_seconds','ELBO']
                                else:
                                     	#If the file type is output file, the dictionary keys can be extracted from the file it$
                                        for i in range(0,N_rowitems):
                                                outputdict[row[i]]=[] #Makes a dummy row of zeros as the first row
                                                variableNames.append(row[i])
                                                #print("Insert dictkey %s"%(row[i]))
                                isFirstrow=1 #all dictionary keys have been inserted to the dict
                        else:
                             	#print('Inserting values to dictionary')
                                for i in range(0,N_rowitems):
                                        #print("Insert value %s to dictionary key %s"%(row[i],variableNames[i]))
                                        outputdict[variableNames[i]].append(float(row[i]))
                                        #print("Now the entries for %s are %s"%(variableNames[i],outputdict[variableNames[i]]))

        return outputdict;

#Parsing samples of theta from extracted dictionary (from variational analysis)
def extractVariable2dim(outputdict,varName,N,M,N_samples):
        #Input arguments:
        #outputdict= the dictionary object from which we should parse the output
        #varName= the name of the variable in the dictionary, a string object
        #N= size in dimension 1, integer
        #M= size in dimension 2, integer, variable is originally a NxM matrix
        #N_samples= number of samples per variable in outputdict
        #This function is for extracting matrices (e.g. 2-dimensional arrays)!!
        extractVar=numpy.zeros((N_samples+1,M,N)) #Make a 3D-matrix where all samples of element (i,j) are stored into the thir$
        for i in range(0,N): #Go through all variable indices, for which the samples have been stored under separate dictionary$
                for j in range(0,M):
                        t_i_j=varName+'.'+str(j+1)+'.'+str(i+1) #Forming the key
                        extractVar[:,j,i]=outputdict[t_i_j] #Store the samples under the key into the matrix

        return extractVar;

#Parsing samples of theta from extracted dictionary (from variational analysis)
def extractVariable1dim(outputdict,varName,N,N_samples):
        #Input arguments:
        #outputdict= the dictionary object from which we should parse the output
        #varName= the name of the variable in the dictionary, a string object
        #N= size in dimension 1, integer, variable is originally a Nx1 vector
        #N_samples= number of samples per variable in outputdict
        #This function is for extracting vectors (e.g. 1-dimensional arrays)!!
        extractVar=numpy.zeros((N_samples+1,N)) #Make a 2D-matrix where all samples are stored
        for i in range(0,N): #Go through all variable indices, for which the samples have been stored under separate dictionary$
                if N==1:
                        t_i=varName #If the original variable was a scalar, this is used as a key
                else:
                     	t_i=varName+'.'+str(i+1) #Forming the key
                extractVar[:,i]=outputdict[t_i] #Store the samples under the key into the matrix

        return extractVar;




def run_luxus_HMC(sigmaB2, lux_data, stan_model_file, plot_file_name, covariate_to_be_tested, N_outputsamples):

	time_start=time.time()
	fit = stan_cache(stan_model_file, data=lux_data, iter=N_outputsamples, chains=4)
	samples=fit.extract(permuted=True)
	bf=calculate_savagedickey_kde_1d(numpy.zeros((1,)),sigmaB2*numpy.eye(1),samples['B'][:,covariate_to_be_tested])
	time_end_full=time.time()
	print(fit)

	vars_to_be_printed=['l','sigmaE2','sigmaC2','sigmaR2','B']
	fit.plot(vars_to_be_printed)
	pyplot.tight_layout()

	savefig(plot_file_name)

	sigmaR2_mean=numpy.mean(samples['sigmaR2'])
	sigmaE2_mean=numpy.mean(samples['sigmaE2'])
	sigmaC2_mean=numpy.mean(samples['sigmaC2'])


	runtime=time_end_full-time_start
	
	return bf, runtime, sigmaR2_mean, sigmaC2_mean, sigmaE2_mean;	

def run_luxus_VI(lux_data, model_name, N_gradsamples, N_elbosamples, N_outputsamples,temp_input_data_file_name, temp_output_file_name, sigmaB2):

	time_start=time.time()
	pystan.misc.stan_rdump(lux_data,temp_input_data_file_name)
	subprocess.call("./%s variational grad_samples=%s elbo_samples=%s output_samples=%s data file=%s output file=%s"%(model_name,N_gradsamples,N_elbosamples,N_outputsamples,temp_input_data_file_name,temp_output_file_name),shell=True)
	samples=csvIntoExtractDict(temp_output_file_name,0)

	samples_B=extractVariable1dim(samples,'B',N_predictors,N_outputsamples)
	bf=calculate_savagedickey_kde_1d(numpy.zeros((1,)),sigmaB2*numpy.eye(1),samples_B[:,1])

	time_end_full=time.time()

	samples_sigmaR2=extractVariable1dim(samples,'sigmaR2',1,N_outputsamples)
	samples_sigmaC2=extractVariable1dim(samples,'sigmaC2',1,N_outputsamples)
	samples_sigmaE2=extractVariable1dim(samples,'sigmaE2',1,N_outputsamples)
	
	sigmaR2_mean=numpy.mean(samples_sigmaR2)
	sigmaC2_mean=numpy.mean(samples_sigmaC2)
	sigmaE2_mean=numpy.mean(samples_sigmaE2)

	subprocess.call("rm %s"%(temp_input_data_file_name),shell=True)
	subprocess.call("rm %s"%(temp_output_file_name),shell=True)

	runtime=time_end_full-time_start

	return bf, runtime, sigmaR2_mean, sigmaC2_mean, sigmaE2_mean;		




if __name__ == '__main__':
	parser = argparse.ArgumentParser(description='Runs LuxUS model for the given input data and returns a BF for the whole window.')
	parser.add_argument('-d','--input_data',action='store',dest='input_data',type=str,required=True,help='Name of the data file.')
	parser.add_argument('-o','--outputFolder',action='store',dest='outputFolder',type=str,required=True,help='Folder where to store the results. Format /level1/level2')	
	parser.add_argument('-i','--inputFolder',action='store',dest='inputFolder',type=str,required=True,help='Folder where the input data is stored. Format /level1/level2')
	parser.add_argument('-j','--outputFile',action='store',dest='outputFile',type=str,required=True,help='File into which the BFs are written. Will be located in folder specified in -o.')
	parser.add_argument('-c','--test_covariate',action='store',dest='test_covariate',type=int,required=True,help='Covariate to be tested. Give index (in design matrix) starting from 0.')
	parser.add_argument('-b','--sigmaB2',action='store',dest='sigmaB2',type=float,required=False,help='Variance for B. Default value is used if not specified.')
	parser.add_argument('-a','--algorithm',action='store',dest='algorithm',type=int,required=True,help='Give value 0 (use HMC) or 1 (use VI).')
	parser.add_argument('-p','--diagnostic_plots',action='store',dest='diagnostic_plots',type=int,required=False,help='Give value 0 (do not plot sample diagnostics for HMC) or 1 (plot sample diagnostics for HMC). Default value is 0.')
	parser.add_argument('-g','--N_gradsamples',action='store',dest='N_gradsamples',type=int,required=False,help='Number of gradient samples used in VI. Default value is used if not specified.')
	parser.add_argument('-e','--N_elbosamples',action='store',dest='N_elbosamples',type=int,required=False,help='Number of gradient samples used in VI. Default value is used if not specified.')
	parser.add_argument('-v','--N_outputsamples_VI',action='store',dest='N_outputsamples_VI',type=int,required=False,help='Number of posterior samples used in VI. Default value is used if not specified.')
	parser.add_argument('-m','--N_outputsamples_HMC',action='store',dest='N_outputsamples_HMC',type=int,required=False,help='Number of posterior samples per chain used in HMC. Default value is used if not specified.')
	parser.add_argument('-w','--window_index',action='store',dest='window_index',type=int,required=False,help='The index of the window being analysed. If value is not given the BF is saved without window index into the defined output file.')
	parser.add_argument('-x','--sigmaR2_file',action='store',dest='sigmaR2_file',type=str,required=True,help='File for storing sigmaR2 posterior means. REMOVE FROM FINAL VERSION.')
	parser.add_argument('-y','--sigmaC2_file',action='store',dest='sigmaC2_file',type=str,required=True,help='File for storing sigmaC2 posterior means. REMOVE FROM FINAL VERSION.')
	parser.add_argument('-z','--sigmaE2_file',action='store',dest='sigmaE2_file',type=str,required=True,help='File for storing sigmaE2 posterior means. REMOVE FROM FINAL VERSION.')
	parser.add_argument('-t','--timeFile',action='store',dest='timeFile',type=str,required=False,help='File name (and path) for storing computation time. If no file name is given the computation times will not be stored into a file.')

	options = parser.parse_args()


	currenttime=time.localtime()

	default_sigmaB2 = 5
	default_diagnostic_plots = 0 
	
	default_N_gradsamples = 2000
	default_N_elbosamples = 2000
	default_N_outputsamples_VI = 2000
	default_N_outputsamples_HMC = 1000


	if  options.sigmaB2 is None:
		print("sigmaB2 was not specified. Default value %s is used."%(default_sigmaB2))
		sigmaB2=default_sigmaB2
	else:
		sigmaB2=options.sigmaB2


	inputf=open("%s/%s"%(options.inputFolder,options.input_data),'rb')
	luxus_data=pickle.load(inputf)
	inputf.close()

	print("Debug: input data object keys")
	print(luxus_data.keys())
	print(luxus_data)

	input_data_id=options.input_data.split('.')[0]


	if options.algorithm==0:

                if  options.N_outputsamples_HMC is None:
                        print("N_outputsamples_HMC was not specified. Default value %s is used."%(default_N_outputsamples_HMC))
                        N_outputsamples_HMC=default_N_outputsamples_HMC
                else:
                     	N_outputsamples_HMC=options.N_outputsamples_HMC

        	if  options.diagnostic_plots is None:
                	print("diagnostic_plots was not specified. Default value %s is used."%(default_diagnostic_plots))
                	diagnostic_plots=default_diagnostic_plots
        	else:
             		diagnostic_plots=options.diagnostic_plots

		plot_file_name="%s/%s_diagnostic_plots.png"%(options.outputFolder,input_data_id)

		print("Estimating the model parameters with HMC.")
		BF, runtime, sigmaR2_mean, sigmaC2_mean, sigmaE2_mean = run_luxus_HMC(sigmaB2, luxus_data, 'luxus_v3_final.stan', plot_file_name, options.test_covariate, N_outputsamples_HMC)



	if options.algorithm==1:

		if  options.N_gradsamples is None:
                	print("N_gradsamples was not specified. Default value %s is used."%(default_N_gradsamples))
                	N_gradsamples=default_N_gradsamples
        	else:
                	N_gradsamples=options.N_gradsamples

		if  options.N_elbosamples is None:
                        print("N_elbosamples was not specified. Default value %s is used."%(default_N_elbosamples))
                        N_elbosamples=default_N_elbosamples
               	else:
                        N_elbosamples=options.N_elbosamples

                if  options.N_outputsamples_VI is None:
                        print("N_outputsamples_VI was not specified. Default value %s is used."%(default_N_outputsamples_VI))
                        N_outputsamples_VI=default_N_outputsamples_VI
                else:
                     	N_outputsamples_VI=options.N_outputsamples_VI


		print("Estimating the model parameters with VI.")

		variational_temp_data_file="%s/TEMP_store_variational_input_%s_%s_%s_%s_%s_%s_%s.txt"%(options.outputFolder,input_data_id,currenttime[0],currenttime[1],currenttime[2],currenttime[3],currenttime[4],currenttime[5])
		variational_temp_output_file="%s/TEMP_store_variational_results_%s_%s_%s_%s_%s_%s_%s.csv"%(options.outputFolder,input_data_id,currenttime[0],currenttime[1],currenttime[2],currenttime[3],currenttime[4],currenttime[5])

		BF, runtime, sigmaR2_mean, sigmaC2_mean, sigmaE2_mean = run_luxus_VI(luxus_data,"luxus_v3_final", N_gradsamples, N_elbosamples, N_outputsamples_VI, variational_temp_data_file, variational_temp_output_file, sigmaB2)


	if options.window_index is None:
		with open("%s/%s"%(options.outputFolder,options.outputFile),'a+') as fw:
			fw.write("%f\n"%(BF))	

	else:
		with open("%s/%s"%(options.outputFolder,options.outputFile),'a+') as fw:
			fw.write("%f\t%s\n"%(BF,options.window_index))


	with open(options.sigmaR2_file,'a+') as fw:
                        fw.write("%f\t%s\n"%(sigmaR2_mean,options.window_index))
	with open(options.sigmaC2_file,'a+') as fw:
                        fw.write("%f\t%s\n"%(sigmaC2_mean,options.window_index))
	with open(options.sigmaE2_file,'a+') as fw:
                        fw.write("%f\t%s\n"%(sigmaE2_mean,options.window_index))


	print("Calculated Bayes factor is %s"%(BF))
        print("Bayes factor calculation took %s seconds."%(runtime))

		

	if  options.timeFile is None:
		print("File for storing computation time was not specified and computation time will not be stored.")
	else:
		with open(options.timeFile,'a+') as tf:
                        tf.write("%f\t%s\n"%(runtime,options.window_index))
