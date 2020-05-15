import argparse
import sys
import os
import commands
import pystan
import numpy
import scipy
import scipy.stats
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

def run_luxus_HMC(sigmaB2, lux_data, diff, stan_file_name, BF_output_file, name_to_BF_file, isWindow,N_cyt_in_analysis, cytosineIndex, d_file,B_file,plot_file_name,modelNro,N_predictors):

	#cytosineIndex is only needed when N_cyt_in_analysis=1, otherwise any value can be given.

	time_start=time.time()
	fit = stan_cache(stan_file_name, data=lux_data, iter=1000, chains=4)
        time_end=time.time()
        samples=fit.extract(permuted=True)

	if isWindow==1:
		bf=calculate_savagedickey_kde_1d(numpy.zeros((1,)),sigmaB2*numpy.eye(1),samples['B'][:,1])
        	with open(BF_output_file,'a+') as fw:
        		fw.write("%f\t%f\t%s\n"%(bf,diff,name_to_BF_file))
	else:

		for c_ind in range(0,N_cyt_in_analysis):
                        bf = calculate_savagedickey_kde_1d(numpy.zeros((1,)),sigmaB2*numpy.eye(1),samples['B'][:,c_ind*N_predictors+1])
			if N_cyt_in_analysis==1:
				with open(BF_output_file,'a+') as f:
                                        f.write("%f\t%f\t%f\t%s\n"%(bf,diff,cytosineIndex,name_to_BF_file))
			else:
                        	with open(BF_output_file,'a+') as f:
                                	f.write("%f\t%f\t%f\t%s\n"%(bf,diff,c_ind,name_to_BF_file))

       
        time_end_full=time.time()

	print(fit)

	thetas=samples['theta'].copy()
        theta_mean=numpy.mean(thetas,1)

	return theta_mean, time_start, time_end, time_end_full;


def run_luxus_VI(lux_data, model_name, N_cytosines, N_predictors, N_replicates, N_gradsamples, N_elbosamples, N_outputsamples, temp_input_data_file_name, temp_output_file_name, sigmaB2, isWindow, name_to_BF_file, BF_output_file, diff,cytosineIndex,d_file,B_file):

	time_start=time.time()
	pystan.misc.stan_rdump(lux_data,temp_input_data_file_name)
	subprocess.call("./%s variational grad_samples=%s elbo_samples=%s output_samples=%s data file=%s output file=%s"%(model_name,N_gradsamples,N_elbosamples,N_outputsamples,temp_input_data_file_name,temp_output_file_name),shell=True)
        time_end=time.time()
	samples=csvIntoExtractDict(temp_output_file_name,0)	

	if isWindow==1:
		samples_B=extractVariable1dim(samples,'B',N_predictors,N_outputsamples)
        	bf=calculate_savagedickey_kde_1d(numpy.zeros((1,)),sigmaB2*numpy.eye(1),samples_B[:,1])
		with open(BF_output_file,'a+') as fw:
        		fw.write("%f\t%f\t%s\n"%(bf,diff,name_to_BF_file))
	else:
		samples_B=extractVariable1dim(samples,'B',N_predictors*N_cytosines,N_outputsamples)
		for c_ind in range(0,N_cytosines):
                        bf = calculate_savagedickey_kde_1d(numpy.zeros((1,)),sigmaB2*numpy.eye(1),samples_B[:,c_ind*N_predictors+1])
                        with open(BF_output_file,'a+') as f:
				if N_cytosines==1:
					f.write("%f\t%f\t%f\t%s\n"%(bf,diff,cytosineIndex,name_to_BF_file))
                                else:
					f.write("%f\t%f\t%f\t%s\n"%(bf,diff,c_ind,name_to_BF_file))

	time_end_full=time.time()


	samples_theta=extractVariable1dim(samples,'theta',2*N_replicates*N_cytosines,N_outputsamples)
        theta_mean=numpy.mean(samples_theta,1)

	subprocess.call("rm %s"%(temp_input_data_file_name),shell=True)
        subprocess.call("rm %s"%(temp_output_file_name),shell=True)

	return theta_mean, time_start, time_end, time_end_full;


if __name__ == '__main__':
	parser = argparse.ArgumentParser(description='Runs LuxUS model on a data set for all the cytosines at the same time.')
	parser.add_argument('-d','--data_id',action='store',dest='DATA_ID',type=int,required=True,help='Data set ID.')
	parser.add_argument('-l','--lengthscale',action='store',dest='l',type=int,required=True,help='Parameter for covariance matrix')
	parser.add_argument('-w','--width',action='store',dest='width',type=int,required=True,help='width of the genomic region where the N_cytosines lie in bp')
	parser.add_argument('-r','--replicates',action='store',dest='replicates',type=int,required=True,help='number of replicates')
	parser.add_argument('-c','--cytosines',action='store',dest='N_cytosines',type=int,required=True,help='number of cytosines')
	parser.add_argument('-p','--predictors',action='store',dest='N_predictors',type=int,required=True,help='number of predictors')
	parser.add_argument('-q','--reads',action='store',dest='reads',type=int,required=True,help='number of BS calls')
	parser.add_argument('-f','--folder',action='store',dest='folder',type=str,required=True,help='Folder where to store the generated data files. Format /level1/level2')	
	parser.add_argument('-i','--inputFolder',action='store',dest='inputFolder',type=str,required=True,help='Folder where the input data is stored. Format /level1/level2')
	parser.add_argument('-b','--sigmaB2',action='store',dest='sigmaB2',type=float,required=True,help='Variance for B (set value)')
	parser.add_argument('-a','--diff',action='store',dest='diff',type=int,required=True,help='Differential or not. Give value 0 or 1.')
	parser.add_argument('-x','--run_old_HMC',action='store',dest='run_old_HMC',type=int,required=True,help='Give value 0 (do not run HMC for old model) or 1. (run HMC for old model)')
	parser.add_argument('-y','--run_old_VI',action='store',dest='run_old_VI',type=int,required=True,help='Give value 0 (do not run VI for old model) or 1. (run VI for old model)')
	parser.add_argument('-z','--run_v2',action='store',dest='run_v2',type=int,required=True,help='Give value 0 (do not run v2 model) or 1. (run v2 model)')
	parser.add_argument('-t','--run_v3',action='store',dest='run_v3',type=int,required=False,help='Give value 0 (do not run v3 model) or 1. (run v3 model)')
	parser.add_argument('-s','--run_v3_single',action='store',dest='run_v3_single',type=int,required=False,help='Run v3 for every cytosine separately. Give value 0 (do not run v3 model) or 1. (run v3 model)')
	
	options = parser.parse_args()

	N_iter=1 #how many times the calculation of BF is repeated for this dataset
	N_cytosines=options.N_cytosines
	N_predictors=options.N_predictors
	reads=options.reads
	replicates=options.replicates
	replicates_2=replicates*2
	l=options.l
	width=options.width
	outputfolder=options.folder
	inputfolder=options.inputFolder

	for diff in [options.diff]:
		print("Differential:%s"%(diff))

		id="C%s_Q%s_R%s_D%s_l%s_w%s"%(N_cytosines,reads,replicates,diff,l,width)
		idv2="C%s_Q%s_R%s_W%s"%(N_cytosines,reads,replicates,width)
		


		#Input files template
		currenttime=time.localtime()
		inputfAllv3="%s/luxusv3_all_%s_set%s_diff%s.pickle"%(inputfolder,idv2,options.DATA_ID,diff)
		inputfAllv2="%s/luxusv2_all_%s_set%s_diff%s.pickle"%(inputfolder,idv2,options.DATA_ID,diff)
		inputfAll="%s/v2generated_all_%s_set%s_diff%s.pickle"%(inputfolder,idv2,options.DATA_ID,diff)
		outputAll="/scratch/work/hallav1/LuxUS/ROC/generated_data/thetas_from_partial_all_%s_%s.txt"%(id,options.DATA_ID)
		bf_output_file="%s/sigmaEtest_bf_%s_%s.txt"%(outputfolder,id,options.DATA_ID)
		bf_window_output_file="%s/sigmaEtest_bf_window_%s_%s.txt"%(outputfolder,id,options.DATA_ID)
		B_output_file="%s/B_values_%s_%s"%(outputfolder,id,options.DATA_ID)
		variational_temp_data_file="%s/TEMP_store_variational_input_%s_%s_%s_%s_%s_%s_%s_%s.txt"%(outputfolder,id,options.DATA_ID,currenttime[0],currenttime[1],currenttime[2],currenttime[3],currenttime[4],currenttime[5])
                variational_temp_output_file="%s/TEMP_store_variational_results_%s_%s_%s_%s_%s_%s_%s_%s.csv"%(outputfolder,id,options.DATA_ID,currenttime[0],currenttime[1],currenttime[2],currenttime[3],currenttime[4],currenttime[5])
		variational_temp_data_file_v2="%s/TEMP_store_variational_input_v2_%s_%s_%s_%s_%s_%s_%s_%s.txt"%(outputfolder,id,options.DATA_ID,currenttime[0],currenttime[1],currenttime[2],currenttime[3],currenttime[4],currenttime[5])                
		variational_temp_output_file_v2="%s/TEMP_store_variational_results_v2_%s_%s_%s_%s_%s_%s_%s_%s.csv"%(outputfolder,id,options.DATA_ID,currenttime[0],currenttime[1],currenttime[2],currenttime[3],currenttime[4],currenttime[5])
		variational_temp_data_file_v3="%s/TEMP_store_variational_input_v3_%s_%s_%s_%s_%s_%s_%s_%s.txt"%(outputfolder,id,options.DATA_ID,currenttime[0],currenttime[1],currenttime[2],currenttime[3],currenttime[4],currenttime[5])
                variational_temp_output_file_v3="%s/TEMP_store_variational_results_v3_%s_%s_%s_%s_%s_%s_%s_%s.csv"%(outputfolder,id,options.DATA_ID,currenttime[0],currenttime[1],currenttime[2],currenttime[3],currenttime[4],currenttime[5])
                N_gradsamples=10
                N_elbosamples=200
                N_outputsamples=2000
		variational_temp_init_file="%s/TEMP_store_variational_init_%s_%s_%s_%s_%s_%s_%s_%s.txt"%(outputfolder,id,options.DATA_ID,currenttime[0],currenttime[1],currenttime[2],currenttime[3],currenttime[4],currenttime[5])
		variational_temp_output_file2="%s/TEMP_store_variational_results2_%s_%s_%s_%s_%s_%s_%s_%s.csv"%(outputfolder,id,options.DATA_ID,currenttime[0],currenttime[1],currenttime[2],currenttime[3],currenttime[4],currenttime[5])
		theta_output_file="%s/theta_means_%s_%s.txt"%(outputfolder,id,options.DATA_ID)
		runtime_output_file="%s/runtimes_%s_%s.txt"%(outputfolder,id,options.DATA_ID)
		diagnostic_file_id="%s/diagnostics_file_%s_%s_%s_%s_%s_%s_%s_%s"%(outputfolder,id,options.DATA_ID,currenttime[0],currenttime[1],currenttime[2],currenttime[3],currenttime[4],currenttime[5])

                

		if options.run_v3==1:

			#Load pickled data
                	inputf=open(inputfAllv3,'rb')
        	        luxv3_data=pickle.load(inputf)
	                inputf.close()

			diagnostics_file_v3="%s_%s.txt"%(diagnostic_file_id,"v3")
			diagnostics_file_v3_VI="%s_%s.txt"%(diagnostic_file_id,"v3_VI")

			v3_plot_file_name="%s/plot_samples_v3_%s_%s.png"%(outputfolder,id,options.DATA_ID)

			mean_theta_v3,v3_time_start,v3_time_end,v3_time_end_full =run_luxus_HMC(options.sigmaB2, luxv3_data, diff, 'luxus_v3_final.stan', bf_window_output_file, "v3",1,N_cytosines,0,diagnostics_file_v3,B_output_file,v3_plot_file_name,2,N_predictors)
                        mean_theta_VI_v3,v3_VI_time_start,v3_VI_time_end,v3_VI_time_end_full=run_luxus_VI(luxv3_data,"luxus_v3_final", N_cytosines, N_predictors, replicates, N_gradsamples, N_elbosamples,N_outputsamples, variational_temp_data_file_v3,variational_temp_output_file_v3,options.sigmaB2,1,"v3_VI",bf_window_output_file, diff,0,diagnostics_file_v3_VI,B_output_file)

			#Write the means of thetas into a file
                        for c_ind in range(0,N_cytosines*2*replicates):
                                with open(theta_output_file,'a+') as f:
                                        f.write("%f\t%f\t%s\n"%(mean_theta_VI_v3[c_ind],c_ind,"v3_VI"))
                                        f.write("%f\t%f\t%s\n"%(mean_theta_v3[c_ind],c_ind,"v3"))

			#Write runtimes to file
			with open(runtime_output_file,'a+') as f:
                                if options.run_v3==1:
                                        f.write("%f\t%f\t%s\n"%(v3_VI_time_end-v3_VI_time_start,v3_VI_time_end_full-v3_VI_time_start,"v3_VI"))
                                        f.write("%f\t%f\t%s\n"%(v3_time_end-v3_time_start,v3_time_end_full-v3_time_start,"v3"))




		if options.run_v3_single==1:

			total_time=0
			total_time_VI=0

			total_time_full=0
			total_time_VI_full=0

			theta_ind=0

			for ind_c in range(0,N_cytosines):
				

				inputFileSepv3="%s/v2generated_sep_%s_cyt%s_set%s_diff%s.pickle"%(inputfolder,idv2,ind_c,options.DATA_ID,diff)


				diagnostics_file_v3_sep="%s_%s_%s.txt"%(diagnostic_file_id,"v3_sep",ind_c)
				diagnostics_file_v3_sep_VI="%s_%s_%s.txt"%(diagnostic_file_id,"v3_sep_VI",ind_c)

				#Load pickled data
				inputf=open(inputFileSepv3,'rb')
				luxv3_data_s=pickle.load(inputf)
				inputf.close()

				variational_temp_data_file_v3_sep="%s/TEMP_store_variational_init_v3_sep_%s_%s_%s_%s_%s_%s_%s_%s_%s.txt"%(outputfolder,id,options.DATA_ID,currenttime[0],currenttime[1],currenttime[2],currenttime[3],currenttime[4],currenttime[5],ind_c)
				variational_temp_output_file_v3_sep="%s/TEMP_store_variational_results_v3_sep_%s_%s_%s_%s_%s_%s_%s_%s_%s.csv"%(outputfolder,id,options.DATA_ID,currenttime[0],currenttime[1],currenttime[2],currenttime[3],currenttime[4],currenttime[5],ind_c)

				v3_sep_plot_file_name="%s/plot_samples_v3_sep_cyt_%s_%s_%s.png"%(outputfolder,ind_c,id,options.DATA_ID)

				theta_mean_v3_sep, time_start_db, time_end_db, time_end_full_db = run_luxus_HMC(options.sigmaB2, luxv3_data_s, diff, 'luxus_v3_1cytosine_withRepEffect_final.stan', bf_output_file, "v3_sep_f", 0,1,ind_c, diagnostics_file_v3_sep,B_output_file,v3_sep_plot_file_name,3,N_predictors)

				theta_mean_v3_sep_VI, time_start_VI_db, time_end_VI_db, time_end_VI_full_db = run_luxus_VI(luxv3_data_s, "luxus_v3_1cytosine_withRepEffect_final", 1, N_predictors, replicates, N_gradsamples, N_elbosamples, N_outputsamples, variational_temp_data_file_v3_sep, variational_temp_output_file_v3_sep, options.sigmaB2, 0, "v3_sep_VI", bf_output_file,diff,ind_c,diagnostics_file_v3_sep_VI,B_output_file)

				total_time=total_time+(time_end_db-time_start_db)
				total_time_full=total_time_full+(time_end_full_db-time_start_db)
	
				total_time_VI=total_time_VI+(time_end_VI_db-time_start_VI_db)
                                total_time_VI_full=total_time_VI_full+(time_end_VI_full_db-time_start_VI_db)


				#Write the means of thetas into a file
				for t_ind in range(0,2*replicates):
					with open(theta_output_file,'a+') as f:
						f.write("%f\t%f\t%s\n"%(theta_mean_v3_sep[t_ind],theta_ind,"v3_sep"))
						f.write("%f\t%f\t%s\n"%(theta_mean_v3_sep_VI[t_ind],theta_ind,"v3_sep_VI"))
					theta_ind=theta_ind+1

			#Write runtimes to file
			with open(runtime_output_file,'a+') as f:
				f.write("%f\t%f\t%s\n"%(total_time,total_time_full,"v3_sep"))
				f.write("%f\t%f\t%s\n"%(total_time_VI,total_time_VI_full,"v3_sep_VI"))










