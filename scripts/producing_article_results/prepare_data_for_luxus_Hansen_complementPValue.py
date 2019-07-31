#Script for making input files for LuxUS. The idea is to make windows of fixed size 
#and analyze the cytosines in that window at a time. One input file is produced per window.
#It is assumed that the input data file has cytosines only from one chromosome.

import argparse
import numpy
import numpy.random
import scipy
import scipy.stats
import pickle
import json


def ftest(y,X,testCov):

    X_nested=numpy.delete(X,testCov,axis=1)
    
    b1,_,_,_=numpy.linalg.lstsq(X_nested,y.T)
    b2,_,_,_=numpy.linalg.lstsq(X,y.T)

    RSS1=numpy.sum((y-numpy.dot(X_nested,b1))**2)
    RSS2=numpy.sum((y-numpy.dot(X,b2))**2)

    F=((RSS1-RSS2)/(X.shape[1]-X_nested.shape[1]))/(RSS2/(len(y)-X.shape[1]))

    pvalue=1-scipy.stats.f.cdf(F,X.shape[1]-X_nested.shape[1],len(y)-X.shape[1])

    return pvalue



if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Takes in the data in text file format and then prepares it for running LuxUS. The data objects will be saved to the specified location.')
    parser.add_argument('-i','--BS_data_file_name',action='store',dest='INPUT_NAME',type=str,required=True,help='The input data file name and path.')
    parser.add_argument('-o','--output_folder',action='store',dest='OUTPUT_FOLDER',type=str,required=True,help='The output location. If output folder is not specified input folder is used instead.')
    parser.add_argument('-d','--design_matrix_file',action='store',dest='DESIGN_MATRIX',type=str,required=True,help='The design matrix file (and file path if needed), where there is one row for each sample. Include also the column of ones corresponding to the intercept.')
    parser.add_argument('-w','--window_width',action='store',dest='WIDTH',type=int,required=False,help='Width of the genomic window within which the cytosines have to lie. If not specified, 1000 bp is used as default.')
    parser.add_argument('-c','--N_cytosines',action='store',dest='N_CYTOSINES',type=int,required=False,help='The maximum number of cytosines in a window. If not specified x is used as default.')
    parser.add_argument('-r','--N_replicates',action='store',dest='N_replicates',type=int,required=True,help='Number of replicates.')
    parser.add_argument('-b','--sigmaB2',action='store',dest='sigmaB2',type=float,required=False,help='Variance sigmaB2 for the coefficients B.')
    parser.add_argument('-p','--alpha_l',action='store',dest='alpha_l',type=float,required=False,help='Hyperparameter alpha for prior distribution of l. If not specified, default value x is used.')
    parser.add_argument('-q','--beta_l',action='store',dest='beta_l',type=float,required=False,help='Hyperparameter beta for prior distribution of l. If not specified, default value y is used.')
    parser.add_argument('-m','--alpha_e',action='store',dest='alpha_e',type=float,required=False,help='Hyperparameter alpha for prior distribution of e. If not specified, default value x is used.')
    parser.add_argument('-n','--beta_e',action='store',dest='beta_e',type=float,required=False,help='Hyperparameter beta for prior distribution of e. If not specified, default value y is used.')
    parser.add_argument('-k','--alpha_r',action='store',dest='alpha_r',type=float,required=False,help='Hyperparameter alpha for prior distribution of sigmaR2. If not specified, default value is used.')
    parser.add_argument('-l','--beta_r',action='store',dest='beta_r',type=float,required=False,help='Hyperparameter beta for prior distribution of sigmaR2. If not specified, default value is used.')
    parser.add_argument('-g','--alpha_c',action='store',dest='alpha_c',type=float,required=False,help='Hyperparameter alpha for prior distribution of sigmaC2. If not specified, default value is used.')
    parser.add_argument('-f','--beta_c',action='store',dest='beta_c',type=float,required=False,help='Hyperparameter beta for prior distribution of sigmaC2. If not specified, default value is used.')
    parser.add_argument('-s','--N_required_samples',action='store',dest='N_required_samples',type=int,required=False,help='In how many samples each cytosine must be present for it to be included in the analysis.')
    parser.add_argument('-a','--bsEff',action='store',dest='bsEff',type=str,required=False,help='Bisulphite conversion efficiency in format [bsEff_1,...,bsEff_N_replicates]. If not specified, default values are used.')
    parser.add_argument('-e','--bsbEff',action='store',dest='bsbEff',type=str,required=False,help='Incorrect bisulphite conversion efficiency in in format [bsbEff_1,...,bsbEff_N_replicates]. If not specified, default values are used.')
    parser.add_argument('-j','--seqErr',action='store',dest='seqErr',type=str,required=False,help='Sequencing error in format [seqErr_1,...,seqErr_N_replicates]. If not specified, default values are used.')
    parser.add_argument('-t','--covariate_to_be_tested',action='store',dest='t_covariate',type=int,required=False,help='Give integer value indicating the column which is the covariate to be tested. Assumed to be a binary covariate. Indexing starts from 0. If not specified, default value x is used.')
    parser.add_argument('-u','--required_meth_diff',action='store',dest='req_meth_diff',type=float,required=False,help='Required p-value for the coefficient for the variable given as --covariate_to_be_tested for the cytosine window to be included in the analysis.')
    parser.add_argument('-v','--required_coverage',action='store',dest='required_coverage',type=int,required=False,help='Required average coverage over window for a replicate. If not specified, default value x is used.')
    parser.add_argument('-y','--meanCovFile',action='store',dest='meanCovFile',type=str,required=True,help='File name for storing mean coverages for each window')
    parser.add_argument('-z','--cytNFile',action='store',dest='cytNFile',type=str,required=True,help='File name for storing the number of cytosines in each window.')


    options = parser.parse_args()

    #Define some default parameter values
    default_width=1000
    default_N_cytosines=20
    default_alpha_l=38
    default_beta_l=1
    default_alpha_e=2
    default_beta_e=2
    default_sigmaB2=15
    default_N_required_samples=1 #This should be scaled according to the number of samples available.
    default_alpha_c=2
    default_beta_c=2
    default_alpha_r=2
    default_beta_r=2
    default_bsEff=1 #THESE VALUES REPRESENT "PERFECT" EXPERIMENT
    default_bsbEff=0
    default_seqErr=0
    default_t_covariate=1 #The second column in design matrix
    default_required_coverage=5


    #Import data file
    BSseq_data=numpy.genfromtxt(options.INPUT_NAME,dtype=None,missing_values="NA")
    #Import design matrix
    design_matrix=numpy.loadtxt(options.DESIGN_MATRIX)

    #Print some information of the data

    print("The length of the specified BS-sequencing data file was %s."%(len(BSseq_data)))

    #Print design matrix
    print("Printing the specified design matrix")
    print(design_matrix)

    N_replicates=options.N_replicates


    #Handling specified/unspecified parameters

    if 	options.WIDTH is None:
       	print("No window width was specified. The default width %s bp is used instead."%(default_width))
        window_width=default_width
    else:
       	window_width=options.WIDTH

    if  options.N_CYTOSINES is None:
        print("No number of cytosines was specified. The default value %s is used instead."%(default_N_cytosines))
        window_N_cytosines=default_N_cytosines
    else:
        window_N_cytosines=options.N_CYTOSINES

    if  options.alpha_l is None:
        print("No alpha for prior of l was specified. The default value %s is used instead."%(default_alpha_l))
        alpha_l=default_alpha_l
    else:
	alpha_l=options.alpha_l

    if  options.beta_l is None:
        print("No beta for prior of l was specified. The default value %s is used instead."%(default_beta_l))
        beta_l=default_beta_l
    else:
	beta_l=options.beta_l

    if  options.alpha_e is None:
        print("No alpha for prior of e was specified. The default value %s is used instead."%(default_alpha_e))
        alpha_e=default_alpha_e
    else:
	alpha_e=options.alpha_e

    if  options.beta_e is None:
        print("No beta for prior of e was specified. The default value %s is used instead."%(default_beta_e))
        beta_e=default_beta_e
    else:
	beta_e=options.beta_e

    if  options.N_required_samples is None:
        print("No number of required samples was specified. The default value %s is used instead."%(default_N_required_samples))
        N_required_samples=default_N_required_samples
    else:
        N_required_samples=options.N_required_samples
 
    if  options.alpha_c is None:
        print("No alpha for prior of sigmaC2 was specified. The default value %s is used instead."%(default_alpha_c))
        alpha_c=default_alpha_c
    else:
	alpha_c=options.alpha_c

    if  options.beta_c is None:
        print("No beta for prior of sigmaC2 was specified. The default value %s is used instead."%(default_beta_c))
        beta_c=default_beta_c
    else:
	beta_c=options.beta_c

    if  options.alpha_r is None:
        print("No alpha for prior of sigmaR2 was specified. The default value %s is used instead."%(default_alpha_r))
        alpha_r=default_alpha_r
    else:
	alpha_r=options.alpha_r

    if  options.beta_r is None:
        print("No beta for prior of sigmaR2 was specified. The default value %s is used instead."%(default_beta_r))
        beta_r=default_beta_r
    else:
	beta_r=options.beta_r

    if  options.bsEff is None:
        print("bsEff was not specified. The default value %s is used instead."%(default_bsEff))
        bsEff=numpy.ones(options.N_replicates)*default_bsEff

    else:
        bsEff=numpy.array(json.loads(options.bsEff))
	
    if  options.bsbEff is None:
        print("bsbEff was not specified. The default value %s is used instead."%(default_bsbEff))
        bsbEff=numpy.ones(options.N_replicates)*default_bsbEff

    else:
	bsbEff=numpy.array(json.loads(options.bsbEff))

    if  options.seqErr is None:
        print("seqErr was not specified. The default value %s is used instead."%(default_seqErr))
        seqErr=numpy.ones(options.N_replicates)*default_seqErr

    else:
	seqErr=numpy.array(json.loads(options.seqErr))

    if  options.t_covariate is None:
        print("The covariate to be tested was not specified. The default value %s is used instead."%(default_t_covariate))
        t_covariate=default_t_covariate

    else:
        t_covariate=options.t_covariate

    if  options.sigmaB2 is None:
        print("sigmaB2 was not specified. The default value %s is used instead."%(default_sigmaB2))
        sigmaB2=default_sigmaB2
    else:
	sigmaB2=options.sigmaB2


    if  options.required_coverage is None:
        print("required_coverage was not specified. The default value %s is used instead."%(default_required_coverage))
        required_coverage=default_required_coverage
    else:
	required_coverage=options.required_coverage


    ind_cytosine=0

    N_predictors=design_matrix.shape[1]

    #print(BSseq_data.shape)
    print(BSseq_data)

    count_data=numpy.zeros((len(BSseq_data),N_replicates))
    meth_data=numpy.zeros((len(BSseq_data),N_replicates))


    for k in range(0,N_replicates):
        count_data[:,k]=BSseq_data["f%s"%(k+3)]
        meth_data[:,k]=BSseq_data["f%s"%(k+3+N_replicates)]*BSseq_data["f%s"%(k+3)]

    #Nan values are replaced with zero in meth_data
    meth_data[numpy.isnan(meth_data)]=0


    print("Printing methylation data")
    print(meth_data)

    print("Printing count data")
    print(count_data)

    #INITIALIZATION 

    coordinates_window=numpy.zeros(window_width) #The maximum number of cytosines in a window is window_width
    counts_window=numpy.zeros((window_width,N_replicates))
    methylated_window=numpy.zeros((window_width,N_replicates))

    WINDOW_START_COORD=BSseq_data[0][1]

    RUNNING_INDEX_SEP=0
    RUNNING_INDEX=0

    NEW_WINDOW=0 #Indicator for initializing a new window

    WINDOW_COUNT=0

    CYTOSINE_INDEX_WINDOW=0

    INCLUDED_IN_ANALYSIS=numpy.zeros(len(BSseq_data))
    INCLUDED_IN_ANALYSIS_WINDOW=numpy.zeros(len(BSseq_data))
    INCLUDED_IN_ANALYSIS_WINDOW_SAVED=1

    print("len(BSSeq_data)=%s"%(len(BSseq_data)))

    number_of_cytosines_in_window=numpy.zeros(len(BSseq_data))
    mean_coverage_in_window=numpy.zeros((len(BSseq_data),N_replicates))


    #Loop over the BS-seq file, make data files to be saved
    for i in range(0,len(BSseq_data)):


        if NEW_WINDOW==1:
            print("Initializing new window.")
            coordinates_window=numpy.zeros(window_width) #The maximum number of cytosines in a window is window_width
            counts_window=numpy.zeros((window_width,N_replicates))
            methylated_window=numpy.zeros((window_width,N_replicates))
            WINDOW_START_COORD=BSseq_data[i][1]
            CYTOSINE_INDEX_WINDOW=0
            INCLUDED_IN_ANALYSIS_WINDOW=numpy.zeros(len(BSseq_data))
            NEW_WINDOW=0

        print("Cytosine %s"%(i))        
        
        if len(numpy.intersect1d(numpy.nonzero(count_data[i,:])[0],numpy.nonzero(design_matrix[:,t_covariate])[0]))>=N_required_samples and len(numpy.intersect1d(numpy.nonzero(count_data[i,:])[0],numpy.nonzero(design_matrix[:,t_covariate]==0)[0]))>=N_required_samples: 
            print("Adding cytosine %s to window."%(i))
            counts_window[CYTOSINE_INDEX_WINDOW,:]=count_data[i,:]
            methylated_window[CYTOSINE_INDEX_WINDOW,:]=meth_data[i,:] 
            coordinates_window[CYTOSINE_INDEX_WINDOW]=BSseq_data[i][1]
            INCLUDED_IN_ANALYSIS_WINDOW[i]=1
            CYTOSINE_INDEX_WINDOW=CYTOSINE_INDEX_WINDOW+1


        if (((i==(len(BSseq_data)-1)) and CYTOSINE_INDEX_WINDOW>=1) or (((BSseq_data[i+1][1]-WINDOW_START_COORD)>window_width)) and CYTOSINE_INDEX_WINDOW>=1) or CYTOSINE_INDEX_WINDOW>=window_N_cytosines:

            observed_reps=numpy.nonzero(numpy.sum(counts_window,axis=0))
            mean_methylation_level=numpy.sum(methylated_window[:,observed_reps],axis=0).astype(float)/numpy.sum(counts_window[:,observed_reps],axis=0).astype(float)
            mean_methylation_level=mean_methylation_level[0]

            print("Mean methylation levels for the replicates were")
            print(mean_methylation_level)

            #Restrict the mean methylation levels in y into range [0.00001,0.99999] to prevent infinite values after logit transformation
            mean_methylation_level[numpy.nonzero(mean_methylation_level>0.99999)[0]]=0.99999
            mean_methylation_level[numpy.nonzero(mean_methylation_level<0.00001)[0]]=0.00001

            print("After restricting the values between [0.00001,0.99999]")
            print(mean_methylation_level)

            y=scipy.special.logit(mean_methylation_level)
            
            print("After logit transformation")
            print(y)

            print("Design matrix for the window")
            print(design_matrix[observed_reps[0],:])

            pvalue=ftest(y,design_matrix[observed_reps[0],:],t_covariate)

            average_coverages=numpy.mean(counts_window[0:CYTOSINE_INDEX_WINDOW,:],axis=0)
            N_samples_with_required_coverage_group1=len(numpy.intersect1d(numpy.nonzero(average_coverages>=required_coverage)[0],numpy.nonzero(design_matrix[:,t_covariate])[0]))
            N_samples_with_required_coverage_group2=len(numpy.intersect1d(numpy.nonzero(average_coverages>=required_coverage)[0],numpy.nonzero(design_matrix[:,t_covariate]==0)[0]))

            if pvalue>options.req_meth_diff and N_samples_with_required_coverage_group1>=N_required_samples and N_samples_with_required_coverage_group2>=N_required_samples:             

                print("The F-test p-value for linear model coefficient of the variable of interest was %s and smaller than %s. Saving window %s."%(pvalue,options.req_meth_diff,WINDOW_COUNT))

                mean_coverage_in_window[INCLUDED_IN_ANALYSIS_WINDOW_SAVED-1,:]=numpy.mean(counts_window[0:CYTOSINE_INDEX_WINDOW,:],axis=0)


                print("The data is saved in a format appropriate for LuxUSv3")
                luxus_data={'n_cytosines': int(CYTOSINE_INDEX_WINDOW),'n_replicates': int(N_replicates),'n_predictors': int(N_predictors),
                                     'bsEff': numpy.tile(bsEff,(CYTOSINE_INDEX_WINDOW)), 'bsBEff': numpy.tile(bsbEff,(CYTOSINE_INDEX_WINDOW)), 'seqErr': numpy.tile(seqErr,(CYTOSINE_INDEX_WINDOW)),
                                     'bsC': methylated_window[0:CYTOSINE_INDEX_WINDOW].flatten().astype(int), 'bsTot': counts_window[0:CYTOSINE_INDEX_WINDOW].flatten().astype(int),
                                     'X': numpy.tile(design_matrix,(CYTOSINE_INDEX_WINDOW,1)), 'Z_R': numpy.tile(numpy.arange(1,N_replicates+1),(CYTOSINE_INDEX_WINDOW)).astype(int), 
                                     'Z_C': numpy.kron(numpy.arange(1,CYTOSINE_INDEX_WINDOW+1),numpy.ones((N_replicates))).astype(int), 'alpha': alpha_e, 'beta': beta_e,
                                     'alphaR': alpha_r, 'betaR': beta_r,
                                     'alphal': alpha_l, 'betal': beta_l,
                                     'sigmaB2': sigmaB2, 'alphaC': alpha_c, 'betaC': beta_c,
                                     'coordinates': coordinates_window[0:CYTOSINE_INDEX_WINDOW]}
                print(luxus_data)
                #Saving the generated datasets
                output4=open("%s/input_for_luxus_v3_%s.txt"%(options.OUTPUT_FOLDER,INCLUDED_IN_ANALYSIS_WINDOW_SAVED),'ab+')
                pickle.dump(luxus_data,output4)
                output4.close()


                number_of_cytosines_in_window[INCLUDED_IN_ANALYSIS_WINDOW_SAVED-1]=CYTOSINE_INDEX_WINDOW

                INCLUDED_IN_ANALYSIS=INCLUDED_IN_ANALYSIS+INCLUDED_IN_ANALYSIS_WINDOW_SAVED*INCLUDED_IN_ANALYSIS_WINDOW 
                INCLUDED_IN_ANALYSIS_WINDOW_SAVED+=1

            else:
                print("Window %s is not saved. Either the F-test p-value for linear model coefficient of the variable of interest was %s (bigger than %s) OR the number of samples with required coverage was too low (controls %s, cases %s)."%(WINDOW_COUNT,pvalue,options.req_meth_diff,N_samples_with_required_coverage_group2,N_samples_with_required_coverage_group1))


            #Initialize new window in next loop
            NEW_WINDOW=1
            WINDOW_COUNT=WINDOW_COUNT+1

        else:
            if (((i==(len(BSseq_data)-1)) and CYTOSINE_INDEX_WINDOW==0) or (((BSseq_data[i+1][1]-WINDOW_START_COORD)>window_width)) and CYTOSINE_INDEX_WINDOW==0):
                NEW_WINDOW=1



    outputid1=options.INPUT_NAME.split('.')[0]
    outputid=outputid1.split('/')[-1]

    print("The whole BS-seq file has been processed.")
    numpy.savetxt("%s/%s_in_analysis_indicator.txt"%(options.OUTPUT_FOLDER,outputid),INCLUDED_IN_ANALYSIS,fmt='%10.1f')
    numpy.savetxt("%s"%(options.cytNFile),number_of_cytosines_in_window[0:INCLUDED_IN_ANALYSIS_WINDOW_SAVED],fmt='%10.1f')
    numpy.savetxt("%s"%(options.meanCovFile),mean_coverage_in_window[0:INCLUDED_IN_ANALYSIS_WINDOW_SAVED,:],fmt='%10.5f')
